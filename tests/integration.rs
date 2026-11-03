//! Native CLI tests use independently constructed sequences, never caller models.
use pathallele::{
    catalogue::{Catalogue, Temporary},
    io,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufRead, Write},
    path::{Path, PathBuf},
    process::{Command, Output},
};
struct Rng(u64);
impl Rng {
    fn dna(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                b"ACGT"[(self.0 & 3) as usize] as char
            })
            .collect()
    }
}
fn temporary() -> Temporary {
    Temporary::new(&std::env::temp_dir(), "pathallele integration ").unwrap()
}
fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_pathallele")
}
fn aligner() -> PathBuf {
    std::env::var_os("PATHALLELE_ALIGNER")
        .expect("Set PATHALLELE_ALIGNER and activate the native calling environment")
        .into()
}
fn success(command: &mut Command) -> Output {
    let x = command.output().unwrap();
    assert!(
        x.status.success(),
        "command {command:?}\nstdout {}\nstderr {}",
        String::from_utf8_lossy(&x.stdout),
        String::from_utf8_lossy(&x.stderr)
    );
    x
}
fn write_fixture(
    root: &Path,
    reference: &str,
    nodes: &[(String, String)],
    paths: &[(String, Vec<String>)],
    samples: &[String],
) {
    fs::create_dir_all(root).unwrap();
    let mut text = String::from("H\tVN:Z:1.0\n");
    for (n, s) in nodes {
        text += &format!("S\t{n}\t{s}\n");
    }
    let edges: BTreeSet<_> = paths
        .iter()
        .flat_map(|(_, w)| w.windows(2).map(|p| (p[0].clone(), p[1].clone())))
        .collect();
    for (a, b) in edges {
        text += &format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            &a[..a.len() - 1],
            &a[a.len() - 1..],
            &b[..b.len() - 1],
            &b[b.len() - 1..]
        );
    }
    for (n, w) in paths {
        text += &format!("P\t{n}\t{}\t*\n", w.join(","));
    }
    fs::write(root.join("graph.gfa"), text).unwrap();
    fs::write(root.join("reference.fa"), format!(">ref\n{reference}\n")).unwrap();
    fs::write(root.join("samples.txt"), samples.join("\n") + "\n").unwrap();
}
fn command(input: &Path, output: &Path, threads: usize) -> Command {
    let mut c = Command::new(bin());
    c.arg("call")
        .arg("--graph")
        .arg(input.join("graph.gfa"))
        .arg("--reference")
        .arg(input.join("reference.fa"))
        .arg("--samples")
        .arg(input.join("samples.txt"))
        .arg("--output")
        .arg(output)
        .arg("--aligner")
        .arg(aligner())
        .args(["--threads", &threads.to_string(), "--jobs", "1"]);
    c
}
fn independent_reconstruction(
    path: &Path,
    reference: &str,
) -> (BTreeMap<String, String>, Vec<String>) {
    let mut names = vec![];
    let mut edits: Vec<Vec<(usize, usize, String)>> = vec![];
    let mut signature = vec![];
    for line in io::reader(path).unwrap().lines() {
        let line = line.unwrap();
        if line.starts_with("#CHROM") {
            names = line.split('\t').skip(9).map(str::to_string).collect();
            edits = vec![vec![]; names.len()];
        }
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        let pos = f[1].parse::<usize>().unwrap() - 1;
        let reference = f[3];
        let alts: Vec<_> = f[4].split(',').collect();
        signature.push(format!(
            "{}\t{}\t{}\t{}",
            f[1],
            f[3],
            f[4],
            f[9..].join("\t")
        ));
        for (si, cell) in f[9..].iter().enumerate() {
            let g = cell.split(':').next().unwrap();
            if g == "." {
                continue;
            }
            let g = g.parse::<usize>().unwrap();
            if g == 0 {
                continue;
            }
            let alt = alts[g - 1];
            let mut l = 0;
            let mut r = reference.len();
            let mut a = alt.len();
            while l < r.min(a) && reference.as_bytes()[l] == alt.as_bytes()[l] {
                l += 1;
            }
            while r > l && a > l && reference.as_bytes()[r - 1] == alt.as_bytes()[a - 1] {
                r -= 1;
                a -= 1;
            }
            edits[si].push((pos + l, pos + r, alt[l..a].into()));
        }
    }
    let mut result = BTreeMap::new();
    for (name, mut changes) in names.into_iter().zip(edits) {
        changes.sort();
        let mut result_sequence = String::new();
        let mut cursor = 0;
        for (s, e, alt) in changes {
            assert!(s >= cursor);
            result_sequence += &reference[cursor..s];
            result_sequence += &alt;
            cursor = e;
        }
        result_sequence += &reference[cursor..];
        result.insert(name, result_sequence);
    }
    signature.sort();
    (result, signature)
}
fn simple(root: &Path, empty: bool) -> (String, BTreeMap<String, String>) {
    let mut rng = Rng(71029);
    let r = rng.dna(6000);
    let cuts = [0, 900, 901, 2000, 3000, 3003, 4000, 6000];
    let mut nodes: Vec<_> = cuts
        .windows(2)
        .enumerate()
        .map(|(i, p)| (format!("ref_{i}"), r[p[0]..p[1]].into()))
        .collect();
    let refwalk: Vec<_> = (0..nodes.len()).map(|i| format!("ref_{i}+")).collect();
    let change = if &r[900..901] == "A" { "C" } else { "A" };
    nodes.push(("snv".into(), change.into()));
    nodes.push(("ins".into(), "GGACT".into()));
    let mut alt: Vec<_> = refwalk
        .iter()
        .enumerate()
        .filter(|&(i, _)| (i != 4))
        .map(|(i, n)| if i == 1 { "snv+".into() } else { n.clone() })
        .collect();
    alt.insert(3, "ins+".into());
    if empty {
        alt = refwalk.clone();
    }
    let reverse = alt
        .iter()
        .rev()
        .map(|n| format!("{}-", &n[..n.len() - 1]))
        .collect();
    let samples = vec!["same".into(), "changed".into(), "reverse".into()];
    let paths = vec![
        ("ref#0#chr".into(), refwalk.clone()),
        ("same#1#ctg".into(), refwalk),
        ("changed#1#ctg".into(), alt),
        ("reverse#1#ctg".into(), reverse),
    ];
    write_fixture(root, &r, &nodes, &paths, &samples);
    let changed = if empty {
        r.clone()
    } else {
        format!(
            "{}{change}{}GGACT{}{}",
            &r[..900],
            &r[901..2000],
            &r[2000..3000],
            &r[3003..]
        )
    };
    (
        r.clone(),
        BTreeMap::from([
            ("same".into(), r),
            ("changed".into(), changed.clone()),
            ("reverse".into(), changed),
        ]),
    )
}

#[test]
fn invalid_resource_settings_fail_before_io() {
    for (key, value) in [
        ("--threads", "0"),
        ("--jobs", "0"),
        ("--minimum-anchor", "18"),
        ("--mapping-k", "0"),
        ("--minimum-mapq", "61"),
        ("--max-matrix-gb", "NaN"),
        ("--mapping-batch-bases", "0"),
    ] {
        let x = Command::new(bin())
            .args([
                "call",
                "--graph",
                "absent",
                "--reference",
                "absent",
                "--samples",
                "absent",
                "--output",
                "absent",
                "--aligner",
                "absent",
                key,
                value,
            ])
            .output()
            .unwrap();
        assert_eq!(x.status.code(), Some(2));
        assert!(String::from_utf8_lossy(&x.stderr).contains("Invalid resource"));
    }
}

#[test]
#[ignore = "Requires native tools and PATHALLELE_ALIGNER"]
fn snv_indel_reverse_paths_resume_export_and_portable_validation() {
    let t = temporary();
    let input = t.path.join("input with spaces");
    let output = t.path.join("output with spaces");
    let (reference, truth) = simple(&input, false);
    success(&mut command(&input, &output, 2));
    assert_eq!(
        independent_reconstruction(&output.join("Graph.vcf.gz"), &reference).0,
        truth
    );
    let c = Catalogue::open(&output, reference.len()).unwrap();
    assert_eq!(c.loci.len(), 3);
    drop(c);
    assert_eq!(
        command(&input, &output, 2).output().unwrap().status.code(),
        Some(2)
    );
    success(command(&input, &output, 2).arg("--resume"));
    // Resume the last committed projection after an interrupted VCF stage.
    let manifest_path = output.join("run.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).unwrap()).unwrap();
    manifest["status"] = serde_json::json!("failed");
    manifest["stages"].as_object_mut().unwrap().remove("vcf");
    fs::write(&manifest_path, serde_json::to_vec(&manifest).unwrap()).unwrap();
    fs::remove_file(output.join("Graph.vcf.gz")).unwrap();
    success(command(&input, &output, 2).arg("--resume"));
    assert_eq!(
        independent_reconstruction(&output.join("Graph.vcf.gz"), &reference).0,
        truth
    );
    // A second process must fail while the output lock is held.
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(output.join(".lock"))
        .unwrap();
    fs2::FileExt::lock_exclusive(&lock).unwrap();
    let blocked = command(&input, &output, 2)
        .arg("--resume")
        .output()
        .unwrap();
    assert_eq!(blocked.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&blocked.stderr).contains("Another process owns"));
    drop(lock);
    assert_eq!(
        command(&input, &output, 2)
            .args(["--resume", "--sequence-anchor", "35"])
            .output()
            .unwrap()
            .status
            .code(),
        Some(2)
    );
    let exported = t.path.join("complete.vcf");
    success(
        Command::new(bin())
            .arg("export")
            .arg(&output)
            .arg("--reference")
            .arg(input.join("reference.fa"))
            .arg("--vcf")
            .arg(&exported),
    );
    assert_eq!(independent_reconstruction(&exported, &reference).0, truth);
    success(
        Command::new("bcftools")
            .arg("view")
            .args(["-Oz", "-o"])
            .arg(t.path.join("complete.vcf.gz"))
            .arg(&exported),
    );
    success(
        Command::new("bcftools")
            .args(["index", "-t"])
            .arg(t.path.join("complete.vcf.gz")),
    );
    let relocated = t.path.join("relocated");
    fs::rename(&input, &relocated).unwrap();
    success(
        Command::new(bin())
            .arg("validate")
            .arg(&output)
            .arg("--graph")
            .arg(relocated.join("graph.gfa"))
            .arg("--reference")
            .arg(relocated.join("reference.fa"))
            .arg("--samples")
            .arg(relocated.join("samples.txt")),
    );
    fs::remove_dir_all(output.join("work")).unwrap();
    fs::remove_dir_all(relocated).unwrap();
    let portable = success(
        Command::new(bin())
            .arg("validate")
            .arg(&output)
            .arg("--outputs-only"),
    );
    let portable: serde_json::Value = serde_json::from_slice(&portable.stdout).unwrap();
    let schema: serde_json::Value =
        serde_json::from_slice(&fs::read(output.join("schema.json")).unwrap()).unwrap();
    assert_eq!(portable["schema"], schema["version"]);
    assert_eq!(
        portable["classification_schema"],
        schema["classification"]["version"]
    );
    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(output.join("Graph.loci.tsv.gz"))
        .unwrap();
    f.write_all(b"corrupt").unwrap();
    assert_eq!(
        Command::new(bin())
            .arg("validate")
            .arg(&output)
            .arg("--outputs-only")
            .output()
            .unwrap()
            .status
            .code(),
        Some(2)
    );
}

#[test]
#[ignore = "Requires native tools and PATHALLELE_ALIGNER"]
fn all_reference_has_valid_empty_vcf_and_catalogue() {
    let t = temporary();
    let input = t.path.join("input");
    let output = t.path.join("output");
    let (reference, truth) = simple(&input, true);
    success(&mut command(&input, &output, 2));
    assert_eq!(
        independent_reconstruction(&output.join("Graph.vcf.gz"), &reference).0,
        truth
    );
    assert!(Catalogue::open(&output, reference.len())
        .unwrap()
        .loci
        .is_empty());
    success(
        Command::new(bin())
            .arg("validate")
            .arg(&output)
            .arg("--outputs-only"),
    );
}

fn structural(root: &Path, n: usize, transform: &str) -> (String, BTreeMap<String, String>) {
    let mut rng = Rng(301 + n as u64);
    let (l, m, r) = (rng.dna(4000), rng.dna(n), rng.dna(4000));
    let k = n / 2;
    let a = if &m[k..k + 1] == "A" { "C" } else { "A" };
    let snv = format!("{}{a}{}", &m[..k], &m[k + 1..]);
    let mut nodes: Vec<(String, String)> = vec![
        ("left".into(), l.clone()),
        ("middle".into(), m.clone()),
        ("right".into(), r.clone()),
        ("snv".into(), snv.clone()),
    ];
    let mut paths: Vec<(String, Vec<String>)> = vec![
        (
            "ref#0#chr".into(),
            vec!["left+".into(), "middle+".into(), "right+".into()],
        ),
        (
            "same#1#ctg".into(),
            vec!["left+".into(), "middle+".into(), "right+".into()],
        ),
        (
            "inverted#1#ctg".into(),
            vec!["left+".into(), "middle-".into(), "right+".into()],
        ),
        (
            "snv#1#ctg".into(),
            vec!["left+".into(), "snv+".into(), "right+".into()],
        ),
        (
            "reverse_ref#1#ctg".into(),
            vec!["right-".into(), "middle-".into(), "left-".into()],
        ),
        (
            "reverse_inv#1#ctg".into(),
            vec!["right-".into(), "middle+".into(), "left-".into()],
        ),
    ];
    if transform == "complex" {
        paths[2].1[1] = "snv-".into();
        let mut local = l.clone();
        local.replace_range(1000..1001, if &l[1000..1001] == "A" { "C" } else { "A" });
        nodes.push(("left_local_snv".into(), local));
        paths[2].1[0] = "left_local_snv+".into();
    }
    if transform == "repeat" {
        paths[2].1.insert(1, "middle-".into());
    }
    if transform == "unshared" || transform == "unshared_split" {
        let inverted = io::reverse_complement(&m);
        let cut = if transform == "unshared" {
            inverted.len()
        } else {
            inverted.len() / 3
        };
        let mut names = vec![];
        for (i, chunk) in inverted.as_bytes().chunks(cut).enumerate() {
            let name = format!("reverse_sequence_{i}");
            nodes.push((name.clone(), String::from_utf8(chunk.to_vec()).unwrap()));
            names.push(name);
        }
        paths[2].1 = std::iter::once("left+".into())
            .chain(names.iter().map(|n| format!("{n}+")))
            .chain(std::iter::once("right+".into()))
            .collect();
        paths[5].1 = std::iter::once("right-".into())
            .chain(names.iter().rev().map(|n| format!("{n}-")))
            .chain(std::iter::once("left-".into()))
            .collect();
    }
    if transform == "reorder" {
        nodes.reverse();
        paths.reverse();
    }
    if transform == "rename" {
        let names: BTreeMap<_, _> = nodes
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.clone(), format!("external:segment_{}", 1000 - i)))
            .collect();
        for (n, _) in &mut nodes {
            *n = names[n].clone();
        }
        for (_, walk) in &mut paths {
            for step in walk {
                *step = format!(
                    "{}{}",
                    names[&step[..step.len() - 1]],
                    &step[step.len() - 1..]
                );
            }
        }
    }
    let reference = format!("{l}{m}{r}");
    let inverted = format!("{l}{}{r}", io::reverse_complement(&m));
    let samples = vec![
        "same".into(),
        "inverted".into(),
        "snv".into(),
        "reverse_ref".into(),
        "reverse_inv".into(),
    ];
    write_fixture(root, &reference, &nodes, &paths, &samples);
    let truth = BTreeMap::from([
        ("same".into(), reference.clone()),
        ("inverted".into(), inverted.clone()),
        ("snv".into(), format!("{l}{snv}{r}")),
        ("reverse_ref".into(), reference.clone()),
        ("reverse_inv".into(), inverted),
    ]);
    (reference, truth)
}
#[test]
#[ignore = "Requires native tools and PATHALLELE_ALIGNER"]
fn bounded_inversions_and_graph_representation_invariance() {
    let t = temporary();
    let mut baseline = None;
    for (size, transform) in [
        (200, "original"),
        (1000, "original"),
        (4000, "original"),
        (10000, "original"),
        (4000, "reorder"),
        (4000, "rename"),
        (4000, "unshared"),
        (4000, "unshared_split"),
    ] {
        let root = t.path.join(format!("{size}_{transform}"));
        let input = root.join("input");
        let output = root.join("output");
        let (reference, truth) = structural(&input, size, transform);
        let threads = if transform == "original" { 1 } else { 2 };
        let mut cmd = command(&input, &output, threads);
        if threads == 2 {
            cmd.args([
                "--mapping-batch-bases",
                "1000",
                "--mapping-cache-bytes",
                "1000",
            ]);
        }
        success(&mut cmd);
        if let Ok(root) = std::env::var("PATHALLELE_TEST_ARTIFACTS") {
            let root = Path::new(&root).join(format!("{size}_{transform}"));
            fs::create_dir_all(&root).unwrap();
            for n in ["graph.gfa", "reference.fa", "samples.txt"] {
                fs::copy(input.join(n), root.join(n)).unwrap();
            }
        }
        let (observed, signature) =
            independent_reconstruction(&output.join("Graph.vcf.gz"), &reference);
        assert_eq!(observed, truth, "{size}/{transform}");
        let mut vcf = pathallele::variant::Reader::open(&output.join("Graph.vcf.gz")).unwrap();
        assert!(
            std::iter::from_fn(|| vcf.next_record().unwrap())
                .any(|r| r.info["TYPE"].split(',').any(|v| v == "INV")),
            "{size}/{transform}"
        );
        if size == 4000 && transform == "original" {
            baseline = Some(signature);
        } else if threads == 2 {
            assert_eq!(Some(signature), baseline);
        }
    }
}
#[test]
#[ignore = "Requires native tools and PATHALLELE_ALIGNER"]
fn unresolved_structures_preserve_unique_local_sequence_but_not_competing_copies() {
    let t = temporary();
    for transform in ["complex", "repeat"] {
        let input = t.path.join(transform).join("input");
        let output = t.path.join(transform).join("output");
        structural(&input, 4000, transform);
        success(&mut command(&input, &output, 1));
        let mut table = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(io::reader(&output.join("Graph.structures.tsv.gz")).unwrap());
        assert!(table
            .records()
            .map(|r| r.unwrap())
            .any(|r| &r[0] == "inverted" && r[6].starts_with("UNRESOLVED_")));
        let mut vcf = pathallele::variant::Reader::open(&output.join("Graph.vcf.gz")).unwrap();
        let si = vcf.header.sample_names().get_index_of("inverted").unwrap();
        let mut focal = 0;
        let mut local = 0;
        while let Some(r) = vcf.next_record().unwrap() {
            if transform == "complex" && r.pos == 1001 {
                local += 1;
                assert!(r.gt[si]>0,"A unique forward source block must retain its known SNV away from the inversion junction");
            }
            if r.pos - 1 <= 6000 && 6000 < r.pos - 1 + r.reference.len() {
                focal += 1;
                assert_eq!(r.gt[si], -1);
                assert_ne!(r.ft[si], "PASS");
            }
        }
        assert!(focal > 0);
        if transform == "complex" {
            assert_eq!(local, 1);
        }
    }
}

#[test]
#[ignore = "Requires native tools and PATHALLELE_ALIGNER"]
fn compound_haplotypes_keep_complete_gt_and_unknown_sequence() {
    let t = temporary();
    let input = t.path.join("input");
    let output = t.path.join("output");
    let mut rng = Rng(912);
    let (left, right) = (rng.dna(2000), rng.dna(2000));
    let middle = "A".repeat(60) + "CG";
    let alt = "CG".to_string() + &"T".repeat(59) + "C";
    let snv = format!("{}G{}", &middle[..30], &middle[31..]);
    let partial = format!("{}{}{}", &middle[..10], "N".repeat(10), &middle[20..]);
    let edge = "G".to_string() + &middle[1..];
    let reference = format!("{left}{middle}{right}");
    let nodes = vec![
        ("left".into(), left.clone()),
        ("right".into(), right.clone()),
        ("ref".into(), middle.clone()),
        ("alt".into(), alt.clone()),
        ("snv".into(), snv.clone()),
        ("partial".into(), partial),
        ("edge".into(), edge.clone()),
    ];
    let mut paths = vec![(
        "ref#0#chr".into(),
        vec!["left+".into(), "ref+".into(), "right+".into()],
    )];
    for (name, mid) in [
        ("same", "ref"),
        ("changed", "alt"),
        ("small", "snv"),
        ("unknown", "partial"),
        ("edge", "edge"),
    ] {
        paths.push((
            format!("{name}#1#ctg"),
            vec!["left+".into(), format!("{mid}+"), "right+".into()],
        ));
    }
    paths.push((
        "reverse#1#ctg".into(),
        vec!["right-".into(), "alt-".into(), "left-".into()],
    ));
    let samples: Vec<String> = ["same", "changed", "small", "unknown", "reverse", "edge"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    write_fixture(&input, &reference, &nodes, &paths, &samples);
    success(&mut command(&input, &output, 2));
    let observed = independent_reconstruction(&output.join("Graph.vcf.gz"), &reference).0;
    for (sample, sequence) in [
        ("same", &middle),
        ("small", &snv),
        ("changed", &alt),
        ("reverse", &alt),
        ("edge", &edge),
    ] {
        assert_eq!(observed[sample], format!("{left}{sequence}{right}"));
    }
    let mut vcf = pathallele::variant::Reader::open(&output.join("Graph.vcf.gz")).unwrap();
    let mut compound = None;
    while let Some(r) = vcf.next_record().unwrap() {
        if r.info["SEQTYPE"].split(',').any(|kind| kind == "MNV") {
            compound = Some(r);
        }
    }
    let compound = compound.unwrap();
    assert_eq!(compound.gt[2], 0);
    assert_eq!(compound.gt[3], -1);
    let catalogue = Catalogue::open(&output, reference.len()).unwrap();
    let li = catalogue.index[&compound.id];
    assert!(catalogue.gt(li, 2) > 0);
    assert_eq!(catalogue.projected(li, 2), 0);
    assert_eq!(catalogue.gt(li, 3), -1);
    success(
        Command::new(bin())
            .arg("validate")
            .arg(&output)
            .arg("--outputs-only"),
    );
    let exported = t.path.join("complete.vcf");
    success(
        Command::new(bin())
            .arg("export")
            .arg(&output)
            .arg("--reference")
            .arg(input.join("reference.fa"))
            .arg("--vcf")
            .arg(&exported),
    );
    let mut vcf = pathallele::variant::Reader::open(&exported).unwrap();
    let mut found = false;
    let mut previous = 0;
    while let Some(r) = vcf.next_record().unwrap() {
        assert!(r.pos >= previous);
        previous = r.pos;
        if r.id == compound.id {
            found = true;
            assert!(r.gt[2] > 0);
            assert_eq!(r.gt[3], -1);
            assert_eq!(&r.alts[r.gt[2] as usize - 1][1..], snv);
        }
    }
    assert!(found);
    // The adapter is executed even if no partial source happens to require it.
    let fake = t.path.join("fake-aligner");
    fs::write(&fake, "#!/bin/sh\nexit 0\n").unwrap();
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(&fake, fs::Permissions::from_mode(0o700)).unwrap();
    assert!(!Command::new(bin())
        .arg("doctor")
        .arg("--aligner")
        .arg(fake)
        .output()
        .unwrap()
        .status
        .success());
}

#[test]
#[ignore = "Requires native tools and PATHALLELE_ALIGNER"]
fn structural_groups_pool_sequence_subtypes_without_changing_exact_calls() {
    for (size, expected) in [(84, vec![0, 1, 1]), (120, vec![0, 1, -1])] {
        let t = temporary();
        let input = t.path.join("input");
        let output = t.path.join("output");
        let mut rng = Rng(991377);
        let left = rng.dna(2000);
        let core = rng.dna(size);
        let right = rng.dna(2000);
        let reference = format!("{left}{core}{right}");
        let inverse = pathallele::io::reverse_complement(&core);
        let mut subtype = inverse.clone();
        subtype.replace_range(
            size / 2..size / 2 + 1,
            if &inverse[size / 2..size / 2 + 1] == "A" {
                "C"
            } else {
                "A"
            },
        );
        let nodes = vec![
            ("left".into(), left),
            ("core".into(), core),
            ("right".into(), right),
            ("inverse".into(), inverse),
            ("subtype".into(), subtype),
        ];
        let walk = |middle: &str| vec!["left+".into(), format!("{middle}+"), "right+".into()];
        let paths = vec![
            ("ref#0#chr".into(), walk("core")),
            ("same#1#chr".into(), walk("core")),
            ("inverted#1#chr".into(), walk("inverse")),
            ("subtype#1#chr".into(), walk("subtype")),
        ];
        let samples = vec!["same".into(), "inverted".into(), "subtype".into()];
        write_fixture(&input, &reference, &nodes, &paths, &samples);
        success(&mut command(&input, &output, 2));
        let mut r =
            pathallele::variant::Reader::open(&output.join(pathallele::structure::VCF)).unwrap();
        let mut matching = vec![];
        while let Some(row) = r.next_record().unwrap() {
            if row.info["TYPE"].starts_with("INV") {
                matching.push(row);
            }
        }
        assert_eq!(matching.len(), 1);
        assert_eq!(
            matching[0].gt,
            expected,
            "{}",
            io::reader(&output.join("Graph.loci.vcf.gz"))
                .unwrap()
                .lines()
                .map(|l| l.unwrap())
                .filter(|l| !l.starts_with('#'))
                .collect::<Vec<_>>()
                .join("\n")
        );
        assert_eq!(
            matching[0].info["AC"],
            expected.iter().filter(|&&g| g == 1).count().to_string()
        );
        success(
            Command::new(bin())
                .arg("validate")
                .arg(&output)
                .arg("--outputs-only"),
        );
        let data = Catalogue::open(&output, reference.len()).unwrap();
        let p = pathallele::structure::Plan::new(&output, &data).unwrap();
        let members = output.join(pathallele::structure::MEMBERS);
        let original = io::reader(&members)
            .unwrap()
            .lines()
            .map(|x| x.unwrap())
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        let mut w = io::writer(&members).unwrap();
        w.write_all(original.replace("SEED", "INVALID").as_bytes())
            .unwrap();
        w.finish().unwrap();
        assert!(p.validate(&output, &data, "ref").is_err());
    }
}
