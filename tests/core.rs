use pathallele::{
    alleles,
    anchors::{self, SeedCounts, Uniqueness},
    catalogue::{self, Catalogue, Temporary},
    graph::Graph,
    io,
    model::{Coverage, Edit, Evidence, Ranges},
    project::Projection,
    source::Parent,
    variant,
};
use std::{collections::BTreeSet, fs, path::Path};
struct Rng(u64);
impl Rng {
    fn next(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 as usize) % n
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| b"ACGT"[self.next(4)] as char).collect()
    }
}
fn temporary() -> Temporary {
    Temporary::new(&std::env::temp_dir(), "pathallele-test").unwrap()
}
fn edit(s: usize, e: usize, a: &str, complex: bool) -> Edit {
    Edit {
        start: s,
        end: e,
        alt: a.into(),
        complex,
    }
}

#[test]
fn exact_chain_matches_exhaustive_optima() {
    let mut rng = Rng(19871005);
    for _ in 0..400 {
        let mut runs = BTreeSet::new();
        for _ in 0..8 {
            runs.insert((rng.next(20), rng.next(20), rng.next(6) + 1));
        }
        let runs: Vec<_> = runs.into_iter().collect();
        let mut best = 0;
        let mut common = (1usize << runs.len()) - 1;
        for mask in 0usize..1 << runs.len() {
            let selected: Vec<_> = runs
                .iter()
                .enumerate()
                .filter_map(|(i, &r)| (mask & (1 << i) != 0).then_some(r))
                .collect();
            if selected
                .windows(2)
                .any(|w| w[0].0 + w[0].2 > w[1].0 || w[0].1 + w[0].2 > w[1].1)
            {
                continue;
            }
            let score: usize = selected.iter().map(|x| x.2).sum();
            match score.cmp(&best) {
                std::cmp::Ordering::Greater => {
                    best = score;
                    common = mask;
                }
                std::cmp::Ordering::Equal => common &= mask,
                std::cmp::Ordering::Less => {}
            }
        }
        let actual = anchors::obligatory(&runs).unwrap();
        assert_eq!(
            actual,
            (0..runs.len())
                .map(|i| common & (1 << i) != 0)
                .collect::<Vec<_>>(),
            "{runs:?}"
        );
    }
}
#[test]
fn exponentially_many_ties_do_not_overflow() {
    let mut runs = vec![];
    for i in 0..180 {
        runs.push((i * 10, i * 10, 2));
        runs.push((i * 10 + 3, i * 10, 2));
    }
    runs.push((2000, 2000, 5));
    let chosen = anchors::obligatory(&runs).unwrap();
    assert_eq!(chosen.iter().filter(|&&x| x).count(), 1);
    assert!(chosen[360]);
}
#[test]
fn sparse_coverage_matches_dense_counting() {
    let mut rng = Rng(250100);
    for _ in 0..200 {
        let mut rows = vec![];
        let mut uncertain = vec![];
        let mut expected = [0usize; 97];
        let mut mask = [false; 97];
        for _ in 0..35 {
            let (a, b) = (rng.next(98), rng.next(98));
            let (s, e) = (a.min(b), a.max(b));
            rows.push((s, e));
            for x in &mut expected[s..e] {
                *x += 1;
            }
        }
        for _ in 0..5 {
            let (a, b) = (rng.next(98), rng.next(98));
            let (s, e) = (a.min(b), a.max(b));
            uncertain.push((s, e));
            mask[s..e].fill(true);
        }
        let known = Coverage::new(&[rows], &[uncertain], 97).unwrap();
        for i in 0..97 {
            assert_eq!(
                known.state(0, i),
                if mask[i] { 2 } else { expected[i].min(2) as u8 }
            );
        }
        for _ in 0..20 {
            let (a, b) = (rng.next(98), rng.next(98));
            let (s, e) = (a.min(b), a.max(b));
            assert_eq!(
                known.unique[0].covers(s, e),
                (s..e).all(|p| expected[p] == 1 && !mask[p])
            );
        }
    }
}
#[test]
fn coverage_memory_is_independent_of_reference_length() {
    let c = Coverage::new(&[vec![(1, 4), (8, 10)]], &[vec![(2, 3)]], 1_000_000_000).unwrap();
    assert_eq!(c.unique[0].rows, vec![(1, 2), (3, 4), (8, 10)]);
    assert_eq!(c.multiple[0].rows, vec![(2, 3)]);
}
#[test]
fn exact_uniqueness_index_matches_overlapping_string_oracle() {
    let mut rng = Rng(70012);
    let count = |s: &str, w: &str| {
        s.as_bytes()
            .windows(w.len())
            .filter(|x| *x == w.as_bytes())
            .count()
    };
    for _ in 0..160 {
        let mut r = rng.dna(120);
        let q = r.clone();
        if rng.next(2) == 0 {
            r.replace_range(60..80, &q[10..30]);
        }
        let rc = [SeedCounts::new(&r, 7), SeedCounts::new(&r, 19)];
        let u = Uniqueness::new(&r, &q, &rc);
        for k in [7, 19, 32] {
            for _ in 0..6 {
                let p = rng.next(q.len() - k + 1);
                let word = &q[p..p + k];
                let rev = io::reverse_complement(word);
                assert_eq!(
                    u.oriented(word),
                    word != rev
                        && count(&r, word) == 1
                        && count(&q, word) == 1
                        && count(&r, &rev) == 0
                        && count(&q, &rev) == 0
                );
            }
        }
    }
    assert!(!anchors::occurs_once("AAAAA", "AAAA"));
}
#[test]
fn sequence_anchors_are_strand_invariant_and_bounded() {
    let mut rng = Rng(1307);
    let middle = rng.dna(96);
    let r = format!("{}CG{middle}C", "A".repeat(60));
    let a = format!("CG{}C{middle}A", "T".repeat(59));
    for k in [19, 31, 51] {
        let runs = anchors::exact_runs(&r, &a, k).unwrap();
        assert!(runs.iter().any(|x| (x.r, x.q, x.n) == (62, 62, 96)));
        let reverse =
            anchors::exact_runs(&io::reverse_complement(&r), &io::reverse_complement(&a), k)
                .unwrap();
        let mut back: Vec<_> = reverse
            .iter()
            .map(|x| (r.len() - x.r - x.n, a.len() - x.q - x.n, x.n))
            .collect();
        back.sort();
        assert_eq!(
            runs.iter().map(|x| (x.r, x.q, x.n)).collect::<Vec<_>>(),
            back
        );
    }
    assert!(anchors::exact_runs("ACGTACGTACGTN", "ACGTACGTACGTA", 3)
        .unwrap()
        .is_empty());
}
#[test]
fn indexed_edit_queries_include_insertion_boundaries() {
    let mut rng = Rng(280102);
    for _ in 0..80 {
        let mut edits = vec![];
        for _ in 0..15 {
            let (a, b) = (rng.next(31), rng.next(31));
            let (s, e) = (a.min(b), a.max(b));
            edits.push((s, e, "T".into(), false));
        }
        let ev = Evidence::new(
            "A".repeat(30),
            vec![edits.clone()],
            vec![vec![(0, 30)]],
            vec![vec![]],
            vec![vec![]],
            vec![vec![]],
            vec![vec![(0, 30)]],
        )
        .unwrap();
        for s in 0..=30 {
            for e in s..=30 {
                let mut expected: Vec<_> = edits
                    .iter()
                    .filter(|(a, b, _, _)| {
                        if s == e {
                            *a <= s && s <= *b
                        } else if a == b {
                            s <= *a && *a <= e
                        } else {
                            s < *b && *a < e
                        }
                    })
                    .cloned()
                    .collect();
                expected.sort();
                expected.dedup();
                let found: Vec<_> = ev.samples[0]
                    .overlapping(&ev.pool, s, e)
                    .iter()
                    .map(|&id| ev.pool[id].input())
                    .collect();
                assert_eq!(found, expected);
            }
        }
    }
}
#[test]
fn compound_projection_requires_already_emitted_witnesses() {
    let r = "ACGTCAGT";
    let mm = vec![vec![(2, 3, "A".into(), false), (5, 6, "T".into(), false)]];
    let make = || {
        Projection::new(
            r.into(),
            vec!["s".into()],
            mm.clone(),
            vec![vec![(0, 8)]],
            vec![vec![]],
            vec![vec![]],
            vec![vec![]],
            vec![vec![(0, 8)]],
        )
        .unwrap()
    };
    let mut p = make();
    let x = p.project("compound", 1, 7, false).unwrap();
    assert!(x.gt[0] > 0);
    assert_eq!(x.projected[0], -1);
    assert_eq!(x.export[0], 9);
    let mut p = make();
    p.project("a", 2, 3, true).unwrap();
    assert_eq!(p.project("compound", 1, 7, false).unwrap().projected[0], -1);
    p.project("b", 5, 6, true).unwrap();
    let x = p.project("full", 1, 7, false).unwrap();
    assert!(x.gt[0] > 0);
    assert_eq!(x.projected[0], 0);
    assert!(x.links.contains("a,b"));
}
#[test]
fn competing_copies_stay_missing_despite_a_covering_chain() {
    let mut p = Projection::new(
        "ACGT".into(),
        vec!["s".into()],
        vec![vec![]],
        vec![vec![(0, 4)]],
        vec![vec![(1, 2)]],
        vec![vec![(1, 2)]],
        vec![vec![(1, 2)]],
        vec![vec![(0, 4)]],
    )
    .unwrap();
    let x = p.project("snp", 1, 2, true).unwrap();
    assert_eq!(x.gt, vec![-1]);
    assert_eq!(x.source, vec![2]);
}

#[test]
fn structural_context_does_not_erase_independently_bounded_local_observations() {
    for (chains, expected) in [(vec![(0, 12)], 0), (vec![(0, 6), (6, 12)], -1)] {
        let mut p = Projection::new(
            "ACGTACGTACGT".into(),
            vec!["s".into()],
            vec![vec![]],
            vec![vec![(0, 12)]],
            vec![vec![]],
            vec![vec![]],
            vec![vec![(0, 12)]],
            vec![chains],
        )
        .unwrap();
        let snp = p.project("snp", 5, 6, true).unwrap();
        let locus = p.project("edit", 4, 8, false).unwrap();
        assert_eq!(snp.gt, vec![expected]);
        assert_eq!(locus.gt, vec![expected]);
        if expected < 0 {
            assert_eq!(snp.source, vec![10]);
            assert_eq!(locus.source, vec![10]);
        }
    }
}
#[test]
fn rejected_directional_normalization_does_not_mask_distant_snv() {
    let r = "A".repeat(80);
    let f = vec![edit(1, 2, "G", false), edit(40, 50, "", false)];
    let b = vec![edit(1, 2, "G", false), edit(41, 51, "", false)];
    let pos = alleles::disagreements(&r, 0, &f, &b);
    assert!(!pos.is_empty());
    assert!(!pos.contains(&1));
    assert!(pos.iter().all(|&p| p >= 39));
}
#[test]
fn mixed_complete_units_do_not_require_gapped_alignment() {
    let p = Parent {
        id: "unit".into(),
        start: 0,
        alleles: vec!["GACGTAC".into(), "GTGGTTC".into(), "GACGTTAC".into()],
        anchors: vec![vec![]; 3],
        carriers: vec![vec![]; 3],
        kinds: vec!["REFERENCE".into(), "SEQUENCE".into(), "SEQUENCE".into()],
    };
    let e = alleles::decompose(&[p.clone()]).unwrap();
    assert_eq!(e[0][1].len(), 1);
    assert!(e[0][1][0].complex);
    assert_eq!(
        pathallele::model::replay(
            &p.alleles[0],
            0,
            p.alleles[0].len(),
            &e[0][2].iter().collect::<Vec<_>>()
        )
        .unwrap(),
        p.alleles[2]
    );
}
#[test]
fn physical_replay_detects_collapsed_repeat_copy() {
    let r = "GAAAAAC";
    let parents = vec![
        Parent {
            id: "a".into(),
            start: 0,
            alleles: vec!["GAA".into(), "GAAA".into()],
            carriers: vec![vec![], vec![0]],
            anchors: vec![vec![]; 2],
            kinds: vec!["REFERENCE".into(), "SEQUENCE".into()],
        },
        Parent {
            id: "b".into(),
            start: 3,
            alleles: vec!["AAA".into(), "AAAA".into()],
            carriers: vec![vec![], vec![0]],
            anchors: vec![vec![]; 2],
            kinds: vec!["REFERENCE".into(), "SEQUENCE".into()],
        },
    ];
    let c = Coverage::new(&[vec![(0, 7)]], &[vec![]], 7).unwrap();
    let mask = vec![Ranges::default()];
    let valid = vec![BTreeSet::from([
        edit(1, 1, "A", false),
        edit(4, 4, "A", false),
    ])];
    assert!(alleles::validate_models(&parents, &valid, r, &c, &mask).is_ok());
    assert!(alleles::validate_models(
        &parents,
        &vec![BTreeSet::from([edit(1, 1, "A", false)])],
        r,
        &c,
        &mask
    )
    .is_err());
}

fn small_graph(root: &Path, text: &str) -> Result<Graph, String> {
    fs::write(root.join("graph.gfa"), text).unwrap();
    fs::write(root.join("reference.fa"), ">ref\nACGTTGCA\n").unwrap();
    fs::write(root.join("samples.txt"), "s\n").unwrap();
    Graph::read(
        &root.join("graph.gfa"),
        &root.join("reference.fa"),
        &root.join("samples.txt"),
        None,
        None,
    )
}
const GFA:&str="H\tVN:Z:1.1\nS\talpha\tACGT\nS\tomega\tTGCA\nL\talpha\t+\tomega\t+\t0M\nP\tref#0#chr\talpha+,omega+\t*\nP\ts#1#ctg\talpha+,omega+\t*\n";
#[test]
fn malformed_and_unsupported_gfa_fail_with_line_number() {
    let t = temporary();
    let cases = [
        GFA.replace("\t0M\n", "\t*\n"),
        GFA.replace("alpha+,omega+\t*", "alpha+,omega+\t0M,0M,0M"),
        format!("{GFA}W\ts\t1\textra\t0\t4\t>alpha\n"),
        format!("{GFA}S\tbroken\n"),
        GFA.replace("S\talpha\tACGT", "S\talpha\tACGT\tLN:i:3"),
        format!("{GFA}P\talpha\talpha+\t*\n"),
        GFA.replace("VN:Z:1.1", "VN:Z:2.0"),
    ];
    for case in cases {
        let e = small_graph(&t.path, &case).err().unwrap();
        assert!(e.contains("GFA line"), "{e}");
    }
}
#[test]
fn segment_names_survive_record_order() {
    let t = temporary();
    let a = small_graph(&t.path, GFA).unwrap();
    let b = small_graph(
        &t.path,
        &GFA.replace(
            "S\talpha\tACGT\nS\tomega\tTGCA",
            "S\tomega\tTGCA\nS\talpha\tACGT",
        ),
    )
    .unwrap();
    assert_eq!(a.spell(a.sample_paths[0][0]), b.spell(b.sample_paths[0][0]));
    assert_eq!(
        a.nodes[a.paths[a.sample_paths[0][0]].nodes[0] as usize].name,
        b.nodes[b.paths[b.sample_paths[0][0]].nodes[0] as usize].name
    );
}
#[test]
fn unselected_paths_are_still_validated() {
    let t = temporary();
    assert!(small_graph(&t.path, &format!("{GFA}P\tother#1#ctg\talpha+,omega+\t*\n")).is_ok());
    assert!(small_graph(&t.path, &format!("{GFA}P\tother#1#ctg\tomega+,alpha+\t*\n")).is_err());
}
#[test]
fn classification_uses_trimmed_lengths_and_not_net_svlen() {
    assert_eq!(variant::classify("TA", "TG").unwrap(), ("SNV", "SMALL"));
    assert_eq!(variant::classify("T", "TG").unwrap(), ("INS", "SMALL"));
    assert_eq!(
        variant::classify(
            &format!("T{}", "A".repeat(50)),
            &format!("T{}", "C".repeat(50))
        )
        .unwrap(),
        ("MNV", "SV")
    );
    assert_eq!(
        variant::classify(
            &format!("T{}", "A".repeat(51)),
            &format!("T{}", "C".repeat(50))
        )
        .unwrap(),
        ("REPLACEMENT", "SV")
    );
    for bad in ["<INV>", "*", "AN", "", "A"] {
        assert!(variant::classify("A", bad).is_err());
    }
}
#[test]
fn complete_vcf_roundtrip_preserves_missingness_and_empty_dimensions() {
    let t = temporary();
    let mut p = Projection::new(
        "ACGT".into(),
        vec!["ref".into(), "alt".into(), "missing".into()],
        vec![vec![], vec![(1, 2, "T".into(), false)], vec![]],
        vec![vec![(0, 4)], vec![(0, 4)], vec![]],
        vec![vec![]; 3],
        vec![vec![]; 3],
        vec![vec![]; 3],
        vec![vec![(0, 4)]; 3],
    )
    .unwrap();
    let row = p.project("v0", 1, 2, true).unwrap();
    let mut writer = catalogue::Writer::new(
        &t.path,
        &["ref".into(), "alt".into(), "missing".into()],
        "ref",
        "ACGT",
    )
    .unwrap();
    let classifier = pathallele::classify::Classifier::new("ACGT", 31).unwrap();
    let classes = vec![classifier.annotate(1, 2, "T", "SEQUENCE").unwrap()];
    writer.add("v0", 1, 2, "SNP", &row, &classes).unwrap();
    writer.finish(2).unwrap();
    let c = Catalogue::open(&t.path, 4).unwrap();
    assert_eq!(
        (0..3).map(|si| c.gt(0, si)).collect::<Vec<_>>(),
        vec![0, 1, -1]
    );
    assert_eq!(c.status(0, 2, false), "NO_DATA");
    drop(c);
    let t = temporary();
    catalogue::Writer::new(&t.path, &["s".into()], "ref", "ACGT")
        .unwrap()
        .finish(1)
        .unwrap();
    assert!(Catalogue::open(&t.path, 4).unwrap().loci.is_empty());
}
#[test]
fn work_checkpoint_rejects_corruption_and_trailing_data() {
    let t = temporary();
    let p = t.path.join("work.bin");
    let expected = vec!["ACGT".to_string(), "GATTACA".to_string()];
    io::checkpoint(&p, &expected).unwrap();
    assert_eq!(io::load_checkpoint::<Vec<String>>(&p).unwrap(), expected);
    let original = fs::read(&p).unwrap();
    for bytes in [
        original[..5].to_vec(),
        original[..original.len() - 1].to_vec(),
        [original.as_slice(), b"trailing"].concat(),
    ] {
        fs::write(&p, bytes).unwrap();
        assert!(io::load_checkpoint::<Vec<String>>(&p).is_err());
    }
    let mut bad = original;
    bad[0] = b'X';
    fs::write(&p, bad).unwrap();
    assert!(io::load_checkpoint::<Vec<String>>(&p).is_err());
}

#[test]
fn paf_reverse_coordinates_and_exact_cigar_contract() {
    let t = temporary();
    let p = t.path.join("a.paf");
    let text = "s\t200\t20\t120\t-\tr\t500\t100\t201\t99\t101\t60\ttp:A:P\tcg:Z:50=1D50=\n";
    fs::write(&p, text).unwrap();
    let records = pathallele::context::read_alignments(&p, 20).unwrap();
    let r = &records["s"][0];
    assert_eq!((r.qs, r.qe, r.native_qs, r.native_qe), (80, 180, 20, 120));
    assert_eq!(r.index.runs, vec![(80, 130, 100), (130, 180, 151)]);
    fs::write(&p, text.replace("50=1D50=", "100M")).unwrap();
    assert!(pathallele::context::read_alignments(&p, 20).is_err());
}

#[test]
fn complete_vcf_preserves_boundary_padding_projection_and_reference_only_loci() {
    use pathallele::classify::Classifier;
    use pathallele::project::ProjectionRow;
    let t = temporary();
    let samples = vec!["ref".into(), "alt".into(), "unknown".into()];
    let reference = "ACGT";
    let classifier = Classifier::new(reference, 31).unwrap();
    let mut writer = catalogue::Writer::new(&t.path, &samples, "ref", reference).unwrap();
    let row = |alleles: Vec<String>, gt, projected| ProjectionRow {
        alleles,
        alts: vec![],
        cells: String::new(),
        links: String::new(),
        status_order: vec![0, 1],
        status_counts: vec![],
        gt,
        projected,
        source: vec![0, 0, 1],
        export: vec![0, 0, 1],
    };
    // Unsorted construction, terminal empty REF/ALT and a projection that differs from GT.
    writer
        .add(
            "right",
            4,
            4,
            "EDIT",
            &row(vec!["".into(), "GG".into()], vec![0, 1, -1], vec![0, 0, -1]),
            &[classifier.annotate(4, 4, "GG", "SEQUENCE").unwrap()],
        )
        .unwrap();
    writer
        .add(
            "left",
            0,
            1,
            "EDIT",
            &row(vec!["A".into(), "".into()], vec![0, 1, -1], vec![0, 1, -1]),
            &[classifier.annotate(0, 1, "", "SEQUENCE").unwrap()],
        )
        .unwrap();
    writer
        .add(
            "reference",
            1,
            2,
            "SNP",
            &row(vec!["C".into()], vec![0, 0, -1], vec![0, 0, -1]),
            &[],
        )
        .unwrap();
    writer.finish(2).unwrap();
    let c = Catalogue::open(&t.path, reference.len()).unwrap();
    assert_eq!(c.loci, ["left", "reference", "right"]);
    assert_eq!(c.alleles[c.index["left"]], ["A", ""]);
    assert_eq!(c.alleles[c.index["right"]], ["", "GG"]);
    assert_eq!(c.gt(c.index["right"], 1), 1);
    assert_eq!(c.projected(c.index["right"], 1), 0);
    assert_eq!(c.status(0, 2, false), "NO_DATA");
    assert!(!t.path.join("Graph.locus_genotypes.npz").exists());
    drop(c);
    use std::io::{Read, Write};
    let mut original = String::new();
    io::reader(&t.path.join(catalogue::FILE))
        .unwrap()
        .read_to_string(&mut original)
        .unwrap();
    for (before, after) in [
        ("1:PASS:0:PASS", "5:PASS:0:PASS"),
        (".:NO_DATA:.:NO_DATA", "0:NO_DATA:.:NO_DATA"),
        ("LSTART=4", "LSTART=5"),
        ("AN=2", "AN=3"),
        ("GT:FT:PG:PF", "GT:FT:PG"),
        ("complete_locus\n", "event_projection\n"),
    ] {
        assert!(original.contains(before));
        let mut out = io::writer(&t.path.join(catalogue::FILE)).unwrap();
        out.write_all(original.replace(before, after).as_bytes())
            .unwrap();
        out.finish().unwrap();
        assert!(
            Catalogue::open(&t.path, reference.len()).is_err(),
            "{before} -> {after}"
        );
    }
}
