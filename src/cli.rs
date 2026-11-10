//! One executable, explicit stages, immutable identities and resumable checkpoints.
use crate::{
    alignment::Adapter,
    alleles,
    catalogue::Catalogue,
    classify,
    graph::Graph,
    io,
    model::Result,
    source::{self, SourceConfig, Sources},
    variant,
};
use clap::{Args, Parser};
use fs2::FileExt;
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Parser, Clone)]
#[command(
    name = "pathallele",
    version,
    about = "Call source-resolved haploid variants from GFA paths",
    arg_required_else_help = true,
    help_template = "{name} {version}\n{about}\n\nUsage: {usage}\n\n{all-args}{after-help}",
    after_help = "Example: pathallele -g graph.gfa -r ref.fa -o out\nKeep out/work/ for --resume. See pathallele.1 for file formats."
)]
struct Call {
    #[arg(
        short = 'g',
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        help = "GFA1 graph with explicit S sequences and P paths"
    )]
    graph: PathBuf,
    #[arg(
        short = 'r',
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        help = "Single-sequence reference FASTA; must match a graph path"
    )]
    reference: PathBuf,
    #[arg(
        short = 's',
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        help = "Optional sample subset/order, one ID per line [all non-reference samples]"
    )]
    samples: Option<PathBuf>,
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        help_heading = "Input/output",
        help = "New output directory; use --resume for an existing run"
    )]
    output: PathBuf,
    #[arg(
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        help = "Path-to-sample TSV with path and sample columns [PanSN names]"
    )]
    path_map: Option<PathBuf>,
    #[command(flatten)]
    settings: Settings,
    #[arg(
        long,
        help_heading = "Input/output",
        help = "Resume with identical inputs, executable, tools and settings"
    )]
    resume: bool,
}
#[derive(Args, Clone, Serialize)]
struct Settings {
    #[arg(
        short = 't',
        long,
        value_name = "INT",
        help_heading = "Resources",
        default_value_t = 16,
        help = "Worker threads for mapping, source processing and compression"
    )]
    threads: usize,
    #[arg(
        short = 'j',
        long,
        value_name = "INT",
        help_heading = "Resources",
        default_value_t = 8,
        help = "Concurrent internal allele-alignment jobs"
    )]
    jobs: usize,
    #[arg(
        long,
        value_name = "INT",
        help_heading = "Mapping/anchors",
        default_value_t = 31,
        help = "Minimum source/graph homology anchor length, bp (>=19)"
    )]
    minimum_anchor: usize,
    #[arg(
        short = 'k',
        long,
        value_name = "INT",
        help_heading = "Mapping/anchors",
        default_value_t = 19,
        help = "Minimap2 and fallback exact-seed k-mer length (>=19)"
    )]
    mapping_k: usize,
    #[arg(
        long,
        value_name = "INT",
        help_heading = "Mapping/anchors",
        default_value_t = 1000,
        help = "Maximum gap for graph-context anchor chaining, bp"
    )]
    maximum_graph_gap: usize,
    #[arg(
        long,
        value_name = "INT",
        help_heading = "Mapping/anchors",
        default_value_t = 31,
        help = "Exact-sequence anchor length within source alleles, bp (>=19)"
    )]
    sequence_anchor: usize,
    #[arg(
        short = 'q',
        long,
        value_name = "INT",
        help_heading = "Mapping/anchors",
        default_value_t = 20,
        help = "Minimum mapping quality for source alignments (0..60)"
    )]
    minimum_mapq: u32,
    #[arg(
        long,
        value_name = "STR",
        help_heading = "Input/output",
        help = "Reference GFA path name [infer from reference sequence]"
    )]
    reference_path: Option<String>,
    #[arg(
        long,
        value_name = "FLOAT",
        help_heading = "Resources",
        default_value_t = 8.0,
        help = "Numeric genotype row-buffer budget, GiB; not a process RSS limit"
    )]
    max_matrix_gb: f64,
    #[arg(
        long,
        value_name = "INT",
        help_heading = "Resources",
        default_value_t = 50_000_000,
        help = "Minimap2 query mini-batch size, bases"
    )]
    mapping_batch_bases: usize,
    #[arg(
        long,
        value_name = "INT",
        help_heading = "Resources",
        default_value_t = 50_000_000,
        help = "Release minimap2 per-thread cache above this many bytes; 0 disables"
    )]
    mapping_cache_bytes: usize,
}
fn canonical(path: &Path) -> Result<PathBuf> {
    fs::canonicalize(path).map_err(|e| format!("{}: {e}", path.display()))
}
fn read_json(path: &Path) -> Result<Value> {
    serde_json::from_reader(File::open(path).map_err(|e| format!("{}: {e}", path.display()))?)
        .map_err(|e| format!("{}: {e}", path.display()))
}
fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |t| t.as_secs())
}
fn executable(name: &Path) -> Result<PathBuf> {
    executable_on_path(name, &std::env::var_os("PATH").unwrap_or_default())
}
fn executable_on_path(name: &Path, search: &std::ffi::OsStr) -> Result<PathBuf> {
    use std::os::unix::fs::PermissionsExt;
    let candidates: Vec<_> = if name.components().count() > 1 {
        vec![name.to_path_buf()]
    } else {
        std::env::split_paths(search)
            .map(|p| p.join(name))
            .collect()
    };
    for p in candidates {
        if let Ok(m) = p.metadata() {
            if m.is_file() && m.permissions().mode() & 0o111 != 0 {
                return canonical(&p);
            }
        }
    }
    Err(format!("Executable not found: {}", name.display()))
}

// The installed caller owns its alignment engine; no external backend discovery.
fn resolve_aligner() -> Result<Adapter> {
    if !cfg!(feature = "native-aligner") {
        return Err("Calling requires the native-aligner build. Install the Bioconda package or build with --features native-aligner.".into());
    }
    let caller = std::env::current_exe().map_err(|e| e.to_string())?;
    Ok(Adapter {
        executable: canonical(&caller)?,
    })
}
fn libraries(binary: &Path) -> Result<BTreeMap<String, String>> {
    let mut command = if cfg!(target_os = "macos") {
        let mut c = Command::new("otool");
        c.arg("-L");
        c
    } else {
        Command::new("ldd")
    };
    let run = command.arg(binary).output().map_err(|e| e.to_string())?;
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&run.stdout),
        String::from_utf8_lossy(&run.stderr)
    );
    if text.contains("not found") {
        return Err(format!("Unresolved dynamic dependency: {text}"));
    }
    if !run.status.success()
        && !text.contains("not a dynamic executable")
        && !text.contains("statically linked")
    {
        return Err(format!("Cannot inspect dependencies: {text}"));
    }
    let mut result = BTreeMap::new();
    for word in text.split_whitespace() {
        let word = word.trim_end_matches(':');
        if word.starts_with('/') && Path::new(word).is_file() {
            let path = canonical(Path::new(word))?;
            if path != binary {
                result.insert(path.to_string_lossy().into_owned(), io::digest(&path)?);
            }
        }
    }
    Ok(result)
}
fn runtime(aligner: &Adapter) -> Result<Value> {
    let mut result = serde_json::Map::new();
    for name in ["minimap2", "bcftools", "samtools"] {
        let path = executable(Path::new(name))?;
        let run = Command::new(&path)
            .arg("--version")
            .output()
            .map_err(|e| e.to_string())?;
        if !run.status.success() {
            return Err(format!("{name} --version failed"));
        }
        let text = format!(
            "{}{}",
            String::from_utf8_lossy(&run.stdout),
            String::from_utf8_lossy(&run.stderr)
        );
        result.insert(name.into(),json!({"path":path,"sha256":io::digest(&path)?,"version":text.lines().next().unwrap_or(""),"libraries":libraries(&path)?}));
    }
    let path = &aligner.executable;
    result.insert(
        "align_alleles".into(),
        json!({"path":path,"sha256":io::digest(path)?,"libraries":libraries(path)?}),
    );
    result.get_mut("align_alleles").unwrap()["backend"] = json!("embedded-vcflib-wfa2");
    Ok(Value::Object(result))
}
fn code() -> Result<Value> {
    Ok(
        json!({"binary_sha256":io::digest(&std::env::current_exe().map_err(|e|e.to_string())?)?,"cargo_lock_sha256":io::hash_bytes(include_bytes!("../Cargo.lock")),"architecture":std::env::consts::ARCH,"os":std::env::consts::OS}),
    )
}
fn identities(root: &Path, names: &[&str]) -> Result<Value> {
    let mut hashes = serde_json::Map::new();
    for &name in names {
        hashes.insert(name.into(), json!(io::digest(&root.join(name))?));
    }
    Ok(Value::Object(hashes))
}
fn verify_identity_maps(root: &Path, values: &[&Value]) -> Result<()> {
    let mut expected = BTreeMap::new();
    for value in values {
        for (name, hash) in value.as_object().ok_or("Invalid manifest identity map")? {
            let hash = hash.as_str().ok_or("Invalid manifest SHA256")?;
            if expected.insert(name, hash).is_some_and(|old| old != hash) {
                return Err(format!("Conflicting manifest identities: {name}"));
            }
        }
    }
    for (name, hash) in expected {
        let relative = Path::new(name);
        if relative
            .components()
            .any(|c| !matches!(c, Component::Normal(_)))
        {
            return Err("Invalid manifest relative path".into());
        }
        if io::digest(&root.join(relative))? != hash {
            return Err(format!("Recorded file identity differs: {name}"));
        }
    }
    Ok(())
}
pub fn verify_manifest(output: &Path, manifest: &Value, outputs_only: bool) -> Result<()> {
    if manifest["status"] != "complete" {
        return Err("Run is not complete".into());
    }
    if manifest["outputs"].as_object().is_none_or(|v| v.is_empty()) {
        return Err("Manifest has no completed output identities".into());
    }
    let mut identities = vec![&manifest["outputs"]];
    if !outputs_only {
        identities.extend(
            manifest["stages"]
                .as_object()
                .ok_or("Missing stages")?
                .values(),
        );
    }
    verify_identity_maps(output, &identities)
}

struct Meter {
    start: Instant,
    own: f64,
    child: f64,
}
fn usage(who: i32) -> (f64, u64) {
    let mut r = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    let ok = unsafe { libc::getrusage(who, r.as_mut_ptr()) };
    if ok != 0 {
        return (0.0, 0);
    }
    let r = unsafe { r.assume_init() };
    let cpu = r.ru_utime.tv_sec as f64
        + r.ru_utime.tv_usec as f64 / 1e6
        + r.ru_stime.tv_sec as f64
        + r.ru_stime.tv_usec as f64 / 1e6;
    let scale = if cfg!(target_os = "macos") { 1 } else { 1024 };
    (cpu, r.ru_maxrss as u64 * scale)
}
impl Meter {
    fn new() -> Self {
        Self {
            start: Instant::now(),
            own: usage(libc::RUSAGE_SELF).0,
            child: usage(libc::RUSAGE_CHILDREN).0,
        }
    }
    fn elapsed(&self) -> Value {
        let (own, peak) = usage(libc::RUSAGE_SELF);
        let (child, child_peak) = usage(libc::RUSAGE_CHILDREN);
        json!({"wall_seconds":self.start.elapsed().as_secs_f64(),"cpu_self_seconds":own-self.own,"cpu_children_seconds":child-self.child,"peak_self_bytes":peak,"peak_child_bytes":child_peak,"memory_scope":"Cumulative process/terminated-child high-water marks; not concurrent process-tree RSS"})
    }
}
struct Run {
    output: PathBuf,
    manifest: Value,
    _lock: File,
}
impl Run {
    fn done(&self, name: &str) -> Result<bool> {
        if let Some(value) = self.manifest["stages"].get(name) {
            verify_identity_maps(&self.output, &[value])?;
            Ok(true)
        } else {
            Ok(false)
        }
    }
    fn save(&self) -> Result<()> {
        io::json(&self.output.join("run.json"), &self.manifest)
    }
    fn finish_stage(&mut self, name: &str, files: &[&str], meter: &Meter) -> Result<()> {
        self.manifest["stages"][name] = identities(&self.output, files)?;
        self.manifest["metrics"][name] = meter.elapsed();
        io::progress(format!(
            "[{name}] {:.3} s",
            meter.start.elapsed().as_secs_f64()
        ));
        self.save()
    }
}

fn call(args: Call) -> Result<()> {
    let total = Meter::new();
    let s = &args.settings;
    if s.threads == 0
        || s.jobs == 0
        || s.maximum_graph_gap == 0
        || s.mapping_batch_bases == 0
        || s.minimum_anchor.min(s.mapping_k).min(s.sequence_anchor) < 19
        || s.minimum_mapq > 60
        || !s.max_matrix_gb.is_finite()
        || s.max_matrix_gb <= 0.0
    {
        return Err(
            "Invalid resource/mapping settings; counts must be positive, anchors >=19, MAPQ 0..60"
                .into(),
        );
    }
    let mut inputs = BTreeMap::new();
    for (name, path) in [("graph", &args.graph), ("reference", &args.reference)] {
        let path = canonical(path)?;
        inputs.insert(name, json!({"path":path,"sha256":io::digest(&path)?}));
    }
    for (name, path) in [("samples", &args.samples), ("path_map", &args.path_map)] {
        if let Some(path) = path {
            let path = canonical(path)?;
            inputs.insert(name, json!({"path":path,"sha256":io::digest(&path)?}));
        }
    }
    let aligner = resolve_aligner()?;
    let tools = runtime(&aligner)?;
    let signature = json!({"version":env!("CARGO_PKG_VERSION"),"inputs":inputs,"settings":s,"tools":tools,"code":code()?});
    // Validate all inputs before creating any result files.
    let graph = Graph::read(
        &args.graph,
        &args.reference,
        args.samples.as_deref(),
        s.reference_path.as_deref(),
        args.path_map.as_deref(),
    )?;
    if graph.samples.len() as f64 * 32.0 / 1073741824.0 > s.max_matrix_gb {
        return Err("Numeric genotype row buffers exceed --max-matrix-gb".into());
    }
    fs::create_dir_all(&args.output).map_err(|e| e.to_string())?;
    let output = canonical(&args.output)?;
    let lock = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(output.join(".lock"))
        .map_err(|e| e.to_string())?;
    lock.try_lock_exclusive()
        .map_err(|_| "Another process owns this output directory")?;
    let nonempty = fs::read_dir(&output)
        .map_err(|e| e.to_string())?
        .any(|e| e.map_or(true, |e| e.file_name() != ".lock"));
    let manifest = if nonempty {
        if !args.resume || !output.join("run.json").is_file() {
            return Err(
                "Output is nonempty; use a fresh directory or --resume for the same run".into(),
            );
        }
        let v = read_json(&output.join("run.json"))?;
        if v["signature"] != signature {
            return Err("Resume refused: input, binary, settings or tool identity changed".into());
        }
        v
    } else {
        json!({"signature":signature,"stages":{},"metrics":{},"started_unix_seconds":now()})
    };
    let mut run = Run {
        output: output.clone(),
        manifest,
        _lock: lock,
    };
    if run.manifest["status"] == "complete" {
        verify_manifest(&output, &run.manifest, false)?;
        println!("Verified complete run: {}", output.display());
        return Ok(());
    }
    run.manifest["status"] = json!("running");
    run.save()?;
    io::set_log(&output.join("run.log"))?;
    io::progress(format!(
        "PathAllele {}; {} samples",
        env!("CARGO_PKG_VERSION"),
        graph.samples.len()
    ));
    let result = execute(&mut run, &args, &graph, &total, &aligner);
    if let Err(e) = &result {
        run.manifest["status"] = json!("failed");
        run.manifest["error"] = json!(e);
        let _ = run.save();
        io::progress(format!("[error] {e}"));
    }
    result?;
    println!("Completed: {}", output.display());
    Ok(())
}
fn execute(
    run: &mut Run,
    args: &Call,
    graph: &Graph,
    total: &Meter,
    aligner: &Adapter,
) -> Result<()> {
    let output = run.output.clone();
    let work = output.join("work");
    fs::create_dir_all(&work).map_err(|e| e.to_string())?;
    let s = &args.settings;
    let sources = if run.done("sources")? {
        io::load_checkpoint::<Sources>(&work.join("sources.bin"))?
    } else {
        let meter = Meter::new();
        graph.write_inventory(&output)?;
        fs::write(
            work.join("reference.fa"),
            format!(">{}\n{}\n", graph.chrom, graph.reference),
        )
        .map_err(|e| e.to_string())?;
        let status = Command::new("samtools")
            .arg("faidx")
            .arg(work.join("reference.fa"))
            .status()
            .map_err(|e| e.to_string())?;
        if !status.success() {
            return Err("samtools faidx failed".into());
        }
        let config = SourceConfig {
            threads: s.threads,
            minimum_anchor: s.minimum_anchor,
            mapping_k: s.mapping_k,
            maximum_gap: s.maximum_graph_gap,
            minimum_mapq: s.minimum_mapq,
            mapping_batch_bases: s.mapping_batch_bases,
            mapping_cache_bytes: s.mapping_cache_bytes,
        };
        let value = source::build(graph, &work, &output, &config)?;
        run.finish_stage(
            "sources",
            &[
                "work/reference.fa",
                "work/reference.fa.fai",
                "work/sources.bin",
                "work/paths.fa",
                "work/context.paf",
                "Graph.paths.tsv.gz",
                "Graph.structures.tsv.gz",
                "Graph.nodes.tsv.gz",
                "Graph.input_paths.tsv.gz",
                "summary.context.json",
            ],
            &meter,
        )?;
        value
    };
    if !run.done("projection")? {
        let meter = Meter::new();
        let projection = work.join("projection");
        if projection.exists() {
            fs::remove_dir_all(&projection).map_err(|e| e.to_string())?;
        }
        let prepared = alleles::prepare(
            &sources,
            &graph.reference,
            &graph.chrom,
            &work.join("reference.fa"),
            &projection,
            s.jobs,
            aligner,
            s.sequence_anchor,
        )?;
        variant::project(
            &prepared,
            &sources,
            graph,
            &output,
            &projection,
            s.minimum_anchor,
            s.threads,
        )?;
        run.finish_stage(
            "projection",
            &[
                "Graph.alleles.tsv.gz",
                "Graph.anchors.tsv.gz",
                "Graph.alignment_ambiguity.tsv.gz",
                "Graph.loci.tsv.gz",
                classify::FILE,
                crate::catalogue::FILE,
                "Graph.projection_links.tsv.gz",
                "schema.json",
                "summary.projection.json",
                "work/projection/Graph.raw.vcf",
            ],
            &meter,
        )?;
    };
    if !run.done("vcf")? {
        let meter = Meter::new();
        variant::finish_vcf(&output, &work.join("reference.fa"), &work)?;
        run.finish_stage(
            "vcf",
            &["Graph.vcf.gz", "Graph.vcf.gz.tbi", "Graph.loci.vcf.gz.tbi"],
            &meter,
        )?;
    }
    let data = Catalogue::open(&output, graph.reference.len())?;
    if !run.done("structure")? {
        let meter = Meter::new();
        let plan = crate::structure::Plan::new(&output, &data)?;
        plan.write(
            &output,
            &data,
            &graph.chrom,
            graph.reference.len(),
            s.threads,
        )?;
        io::json(&output.join("summary.structure.json"), &plan.summary())?;
        run.finish_stage(
            "structure",
            &[
                crate::structure::VCF,
                "Graph.sv.vcf.gz.tbi",
                crate::structure::MEMBERS,
                "summary.structure.json",
            ],
            &meter,
        )?;
    }
    run.manifest["status"] = json!("complete");
    run.manifest["finished_unix_seconds"] = json!(now());
    let mut names = vec![];
    for item in fs::read_dir(&output).map_err(|e| e.to_string())? {
        let item = item.map_err(|e| e.to_string())?;
        let name = item.file_name().to_string_lossy().into_owned();
        if item.file_type().map_err(|e| e.to_string())?.is_file()
            && !matches!(name.as_str(), "run.json" | "run.log" | ".lock")
        {
            names.push(name);
        }
    }
    run.manifest["outputs"] = identities(
        &output,
        &names.iter().map(String::as_str).collect::<Vec<_>>(),
    )?;
    run.manifest["metrics"]["total"] = total.elapsed();
    run.manifest.as_object_mut().unwrap().remove("error");
    run.save()?;
    io::progress(format!(
        "[complete] {:.3} s",
        total.start.elapsed().as_secs_f64()
    ));
    Ok(())
}

pub fn run() -> Result<()> {
    // Private process-isolation entry point; never a public subcommand.
    let mut args = std::env::args_os().skip(1);
    if args.next().is_some_and(|arg| arg == "internal-align") {
        let input = args.next().ok_or("Missing internal alignment input")?;
        if args.next().is_some() {
            return Err("Unexpected internal alignment argument".into());
        }
        return crate::alignment::worker(Path::new(&input));
    }
    call(Call::parse())
}

#[cfg(test)]
mod help_tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn calling_owns_its_alignment_engine() {
        let args = ["pathallele", "-g", "g", "-r", "r", "-s", "s", "-o", "o"];
        assert!(Call::try_parse_from(args).is_ok());
        for removed in ["call", "validate", "export", "doctor", "build-aligner"] {
            assert!(Call::try_parse_from(["pathallele", removed]).is_err());
        }
        assert!(Call::try_parse_from(["pathallele", "-g", "g", "-r", "r", "-o", "o"]).is_ok());
        let mut override_args = args.to_vec();
        override_args.extend(["--aligner", "external"]);
        assert!(Call::try_parse_from(override_args).is_err());
        assert!(Call::try_parse_from(["pathallele", "doctor", "--aligner", "external"]).is_err());
    }

    #[test]
    fn engine_availability_matches_build() {
        assert_eq!(resolve_aligner().is_ok(), cfg!(feature = "native-aligner"));
    }

    #[test]
    fn every_command_and_argument_has_help() {
        fn check(command: &clap::Command) {
            assert!(command.get_about().is_some(), "{}", command.get_name());
            for arg in command.get_arguments() {
                assert!(
                    arg.get_help()
                        .is_some_and(|text| !text.to_string().is_empty()),
                    "{}: {} has no help",
                    command.get_name(),
                    arg.get_id()
                );
            }
            for subcommand in command.get_subcommands() {
                check(subcommand);
            }
        }
        Call::command().debug_assert();
        check(&Call::command());
    }

    #[test]
    fn short_options_preserve_call_defaults() {
        let parsed = Call::try_parse_from([
            "pathallele",
            "-g",
            "graph.gfa",
            "-r",
            "ref.fa",
            "-s",
            "samples.txt",
            "-o",
            "out",
        ])
        .unwrap();
        let call = parsed;
        assert_eq!(call.graph, PathBuf::from("graph.gfa"));
        assert_eq!(call.settings.threads, 16);
        assert_eq!(call.settings.jobs, 8);
        assert_eq!(call.settings.minimum_anchor, 31);
        assert_eq!(call.settings.mapping_k, 19);
        assert_eq!(call.settings.sequence_anchor, 31);
        assert_eq!(call.settings.minimum_mapq, 20);
    }
}
