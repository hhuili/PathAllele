//! One executable, explicit stages, immutable identities and resumable checkpoints.
use crate::{
    alleles::{self, Prepared},
    catalogue::{Catalogue, Temporary},
    classify,
    graph::Graph,
    io,
    model::{self, EditInput, Result},
    source::{self, SourceConfig, Sources},
    validate, variant,
};
use clap::{Args, Parser, Subcommand};
use fs2::FileExt;
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Parser)]
#[command(
    name = "pathallele",
    version,
    about = "Call source-resolved haploid variants from GFA paths",
    arg_required_else_help = true,
    disable_help_subcommand = true,
    help_template = "{name} {version}\n{about}\n\nUsage: {usage}\n\n{all-args}{after-help}",
    after_help = "Run 'pathallele <command> -h' for options. See pathallele.1 for file formats."
)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Call exact alleles and structural groups
    #[command(
        after_help = "Example: pathallele call -g graph.gfa -r ref.fa -s samples.txt -o out --aligner pathallele-align\nFull validation is included. Keep out/work/ for revalidation and resume."
    )]
    Call(Call),
    /// Validate a completed run
    Validate(Check),
    /// Export selected complete-locus haplotypes
    Export(Export),
    /// Check native tools and library dependencies
    Doctor {
        #[arg(
            long,
            value_name = "EXE",
            help = "BiWFA adapter executable (path or PATH name)"
        )]
        aligner: PathBuf,
    },
    /// Build the vcflib/BiWFA alignment adapter
    BuildAligner(Build),
}
#[derive(Args, Clone)]
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
        help = "Non-reference haploid sample IDs, one per line"
    )]
    samples: PathBuf,
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
        value_name = "EXE",
        help_heading = "Input/output",
        help = "BiWFA adapter executable (path or PATH name)"
    )]
    aligner: PathBuf,
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
        help = "Concurrent external allele-alignment jobs"
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
#[derive(Args)]
struct Check {
    #[arg(
        value_name = "DIR",
        help = "Completed run directory containing run.json"
    )]
    output: PathBuf,
    #[arg(
        long,
        help = "Check published outputs only; no original graph or work/ required"
    )]
    outputs_only: bool,
    #[arg(
        short = 'g',
        long,
        value_name = "FILE",
        help = "Relocated original GFA [recorded path; identical SHA256 required]"
    )]
    graph: Option<PathBuf>,
    #[arg(
        short = 'r',
        long,
        value_name = "FILE",
        help = "Relocated original reference FASTA [recorded path]"
    )]
    reference: Option<PathBuf>,
    #[arg(
        short = 's',
        long,
        value_name = "FILE",
        help = "Relocated original sample list [recorded path]"
    )]
    samples: Option<PathBuf>,
    #[arg(
        long,
        value_name = "FILE",
        help = "Relocated original path-to-sample TSV [recorded path]"
    )]
    path_map: Option<PathBuf>,
}
#[derive(Args)]
struct Export {
    #[arg(value_name = "DIR", help = "Completed run directory")]
    output: PathBuf,
    #[arg(
        short = 'r',
        long,
        value_name = "FILE",
        help = "Original single-sequence reference FASTA"
    )]
    reference: PathBuf,
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        help = "New uncompressed VCF; existing files are not overwritten"
    )]
    vcf: PathBuf,
    #[arg(
        short = 'l',
        long,
        value_name = "ID",
        help = "Complete locus ID; repeat to select multiple loci [all variant loci]"
    )]
    locus: Vec<String>,
}
#[derive(Args)]
struct Build {
    #[arg(
        short = 'p',
        long,
        value_name = "DIR",
        help = "Native dependency prefix containing include/ and lib/"
    )]
    prefix: PathBuf,
    #[arg(
        long,
        value_name = "DIR",
        help = "WFA2 source headers containing wavefront/wfa.hpp"
    )]
    wfa_headers: PathBuf,
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        help = "New adapter executable; existing files are not overwritten"
    )]
    output: PathBuf,
    #[arg(
        long,
        value_name = "EXE",
        default_value = "c++",
        help = "C++17 compiler executable (path or PATH name)"
    )]
    cxx: String,
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
    use std::os::unix::fs::PermissionsExt;
    let candidates: Vec<_> = if name.components().count() > 1 {
        vec![name.to_path_buf()]
    } else {
        std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default())
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
fn runtime(aligner: &Path) -> Result<Value> {
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
    let path = executable(aligner)?;
    let probe = Temporary::new(&std::env::temp_dir(), "pathallele-probe")?;
    let input = probe.path.join("probe.vcf");
    std::fs::write(&input, "##fileformat=VCFv4.2\n##contig=<ID=ref,length=8>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nref\t1\tprobe\tACGTACGT\tACCTACGT\t.\t.\t.\n").map_err(|e| e.to_string())?;
    let check = Command::new(&path)
        .arg(&input)
        .output()
        .map_err(|e| format!("Cannot execute alignment adapter: {e}"))?;
    if !check.status.success()
        || String::from_utf8_lossy(&check.stdout).trim() != "probe\t3\tG\tC\t0"
    {
        return Err(format!(
            "Alignment adapter functional check failed: {}",
            String::from_utf8_lossy(&check.stderr)
        ));
    }
    result.insert(
        "align_alleles".into(),
        json!({"path":path,"sha256":io::digest(&path)?,"libraries":libraries(&path)?}),
    );
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
fn verify_hashes(root: &Path, value: &Value) -> Result<()> {
    for (name, hash) in value.as_object().ok_or("Invalid manifest identity map")? {
        let relative = Path::new(name);
        if relative
            .components()
            .any(|c| !matches!(c, Component::Normal(_)))
        {
            return Err("Invalid manifest relative path".into());
        }
        if io::digest(&root.join(relative))?.as_str()
            != hash.as_str().ok_or("Invalid manifest SHA256")?
        {
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
    verify_hashes(output, &manifest["outputs"])?;
    if !outputs_only {
        for stage in manifest["stages"]
            .as_object()
            .ok_or("Missing stages")?
            .values()
        {
            verify_hashes(output, stage)?;
        }
    }
    Ok(())
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
            verify_hashes(&self.output, value)?;
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
    for (name, path) in [
        ("graph", &args.graph),
        ("reference", &args.reference),
        ("samples", &args.samples),
    ] {
        let path = canonical(path)?;
        inputs.insert(name, json!({"path":path,"sha256":io::digest(&path)?}));
    }
    if let Some(path) = &args.path_map {
        let path = canonical(path)?;
        inputs.insert("path_map", json!({"path":path,"sha256":io::digest(&path)?}));
    }
    let tools = runtime(&args.aligner)?;
    let signature = json!({"version":env!("CARGO_PKG_VERSION"),"inputs":inputs,"settings":s,"tools":tools,"code":code()?});
    // Validate all inputs before creating any result files.
    let graph = Graph::read(
        &args.graph,
        &args.reference,
        &args.samples,
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
        "PathAllele {}; {} samples; standalone Rust",
        env!("CARGO_PKG_VERSION"),
        graph.samples.len()
    ));
    let result = execute(&mut run, &args, &graph, &total);
    if let Err(e) = &result {
        run.manifest["status"] = json!("failed");
        run.manifest["error"] = json!(e);
        let _ = run.save();
        io::progress(format!("[error] {e}"));
    }
    result?;
    println!("Completed and validated: {}", output.display());
    Ok(())
}
fn execute(run: &mut Run, args: &Call, graph: &Graph, total: &Meter) -> Result<()> {
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
    let (prepared, positive) = if run.done("projection")? {
        (
            io::load_checkpoint::<Prepared>(&work.join("prepared.bin"))?,
            io::load_checkpoint::<Vec<Vec<EditInput>>>(&work.join("positive.bin"))?,
        )
    } else {
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
            &args.aligner,
            s.sequence_anchor,
        )?;
        let positive = variant::project(
            &prepared,
            &sources,
            graph,
            &output,
            &projection,
            s.minimum_anchor,
            s.threads,
        )?;
        io::checkpoint(&work.join("prepared.bin"), &prepared)?;
        io::checkpoint(&work.join("positive.bin"), &positive)?;
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
                "work/prepared.bin",
                "work/positive.bin",
                "work/projection/model_validation.json",
                "work/projection/Graph.raw.vcf",
            ],
            &meter,
        )?;
        (prepared, positive)
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
    let meter = Meter::new();
    let validation = validate::full(&output, graph, &sources, &prepared, &positive, &data)?;
    run.manifest["metrics"]["validation"] = meter.elapsed();
    io::json(&output.join("validation.json"), &validation)?;
    run.manifest["validation"] = validation;
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
        "[validation] {:.3} s; total {:.3} s",
        meter.start.elapsed().as_secs_f64(),
        total.start.elapsed().as_secs_f64()
    ));
    Ok(())
}

fn portable_data(output: &Path) -> Result<(Value, Value, Catalogue)> {
    let manifest = read_json(&output.join("run.json"))?;
    verify_manifest(output, &manifest, true)?;
    let schema = read_json(&output.join("schema.json"))?;
    if schema["version"] != variant::OUTPUT_SCHEMA || schema["ploidy"] != 1 {
        return Err("Unsupported output schema; use the original binary for older runs".into());
    }
    let length = schema["reference"]["length"]
        .as_u64()
        .ok_or("Invalid reference length")? as usize;
    let data = Catalogue::open(output, length)?;
    Ok((manifest, schema, data))
}
fn check(args: Check) -> Result<Value> {
    let output = canonical(&args.output)?;
    let (mut manifest, schema, data) = portable_data(&output)?;
    let chrom = schema["reference"]["contig"]
        .as_str()
        .ok_or("Invalid reference contig")?;
    let length = schema["reference"]["length"].as_u64().unwrap() as usize;
    if args.outputs_only {
        return validate::portable(&output, &data, chrom, length);
    }
    if manifest["signature"]["version"] != env!("CARGO_PKG_VERSION")
        || manifest["signature"]["code"] != code()?
    {
        return Err("Full validation requires the recorded PathAllele binary; use --outputs-only for portable files".into());
    }
    verify_manifest(&output, &manifest, false)?;
    let inputs = manifest["signature"]["inputs"]
        .as_object_mut()
        .ok_or("Missing input identities")?;
    for (key, path) in [
        ("graph", args.graph),
        ("reference", args.reference),
        ("samples", args.samples),
        ("path_map", args.path_map),
    ] {
        if let Some(path) = path {
            let entry = inputs
                .get_mut(key)
                .ok_or("Cannot rebind an input absent from the run")?;
            entry["path"] = json!(canonical(&path)?);
        }
    }
    for item in inputs.values() {
        let path = item["path"].as_str().ok_or("Invalid input path")?;
        if item["sha256"] != io::digest(Path::new(path))? {
            return Err(format!("Original input changed: {path}"));
        }
    }
    let path = |key: &str| -> Result<PathBuf> {
        Ok(PathBuf::from(
            inputs
                .get(key)
                .and_then(|x| x["path"].as_str())
                .ok_or("Missing original input")?,
        ))
    };
    let (g, r, s, m) = (
        path("graph")?,
        path("reference")?,
        path("samples")?,
        if inputs.contains_key("path_map") {
            Some(path("path_map")?)
        } else {
            None
        },
    );
    let graph = Graph::read(
        &g,
        &r,
        &s,
        manifest["signature"]["settings"]["reference_path"].as_str(),
        m.as_deref(),
    )?;
    let work = output.join("work");
    validate::full(
        &output,
        &graph,
        &io::load_checkpoint(&work.join("sources.bin"))?,
        &io::load_checkpoint(&work.join("prepared.bin"))?,
        &io::load_checkpoint::<Vec<Vec<EditInput>>>(&work.join("positive.bin"))?,
        &data,
    )
}

fn export(args: Export) -> Result<()> {
    let output = canonical(&args.output)?;
    let (manifest, schema, data) = portable_data(&output)?;
    let reference = canonical(&args.reference)?;
    if manifest["signature"]["inputs"]["reference"]["sha256"] != io::digest(&reference)? {
        return Err("Reference SHA256 differs from run".into());
    }
    let records = io::read_fasta(&reference)?;
    if records.len() != 1 {
        return Err("Expected one reference sequence".into());
    }
    let (chrom, reference) = &records[0];
    if schema["reference"]["contig"] != *chrom || schema["reference"]["length"] != reference.len() {
        return Err("Reference identity differs from schema".into());
    }
    let wanted: HashSet<_> = if args.locus.is_empty() {
        data.loci.iter().cloned().collect()
    } else {
        args.locus.into_iter().collect()
    };
    let annotations = classify::Table::read(&output, &data, Some(reference))?;
    if wanted.iter().any(|s| !data.index.contains_key(s)) {
        return Err("Unknown selected locus".into());
    }
    if args.vcf.exists() {
        return Err("Export destination already exists".into());
    }
    let parent = args
        .vcf
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    let temporary = Temporary::new(parent, "export")?;
    let path = temporary.path.join("complete.vcf");
    let mut f = io::writer(&path)?;
    let mut metadata=vec!["##fileformat=VCFv4.2".into(),format!("##source=PathAllele-{}",env!("CARGO_PKG_VERSION")),format!("##contig=<ID={chrom},length={}>",reference.len()),
        "##pathallele_locus_scope=\"Independent overlapping complete haplotypes; do not replay all records together\"".into(),
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT allele count\">".into(),"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Called haploid samples\">".into(),"##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency\">".into()];
    metadata.extend(
        variant::CLASS_HEADERS
            .iter()
            .filter(|s| !s.starts_with("##pathallele_hierarchy="))
            .map(|&s| s.to_string()),
    );
    metadata.extend(variant::headers("complete_locus"));
    writeln!(
        f,
        "{}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        metadata.join("\n"),
        data.samples.join("\t")
    )
    .map_err(|e| e.to_string())?;
    let mut order: Vec<_> = (0..data.loci.len()).collect();
    order.sort_by_key(|&i| {
        (
            if data.kinds[i] == "SNP" {
                data.bounds[i].0 + 1
            } else {
                data.bounds[i].0.max(1)
            },
            &data.loci[i],
        )
    });
    for li in order {
        let id = &data.loci[li];
        let seqs = &data.alleles[li];
        if !wanted.contains(id) || seqs.len() == 1 {
            continue;
        }
        let (s, e) = data.bounds[li];
        if reference[s..e] != seqs[0] {
            return Err("Complete REF differs from reference".into());
        }
        let (pos, r, alts) =
            model::vcf_alleles(reference, s, e, &seqs[1..], data.kinds[li] == "SNP")?;
        let mut ac = vec![0; alts.len()];
        let mut an = 0;
        let mut cells = vec![];
        for si in 0..data.samples.len() {
            let g = data.gt(li, si);
            if g >= 0 {
                an += 1;
                if g > 0 {
                    ac[g as usize - 1] += 1;
                }
            }
            cells.push(format!(
                "{}:{}",
                if g < 0 { ".".into() } else { g.to_string() },
                data.status(li, si, false)
            ));
        }
        if an == 0 {
            return Err("Complete ALT catalogue has no called samples".into());
        }
        let counts = ac
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let af = ac
            .iter()
            .map(|&v| format!("{:.8}", v as f64 / an as f64))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            f,
            "{chrom}\t{pos}\t{id}\t{r}\t{}\t.\t.\tAC={counts};AN={an};AF={af};{}\tGT:FT\t{}",
            alts.join(","),
            classify::info(
                &(1..seqs.len())
                    .map(|ai| annotations.get(li, ai))
                    .collect::<Vec<_>>()
            ),
            cells.join("\t")
        )
        .map_err(|e| e.to_string())?;
    }
    f.finish()?;
    File::open(&path)
        .map_err(|e| e.to_string())?
        .sync_all()
        .map_err(|e| e.to_string())?;
    fs::hard_link(path, &args.vcf).map_err(|e| e.to_string())?;
    println!("Exported complete locus haplotypes: {}", args.vcf.display());
    Ok(())
}
fn build(args: Build) -> Result<()> {
    let prefix = canonical(&args.prefix)?;
    let headers = canonical(&args.wfa_headers)?;
    if args.output.exists() {
        return Err("Refusing to overwrite an existing native adapter".into());
    }
    for file in [
        prefix.join("include/vcflib/Variant.h"),
        headers.join("wavefront/wfa.hpp"),
    ] {
        if !file.is_file() {
            return Err(format!("Missing development header: {}", file.display()));
        }
    }
    let compiler = executable(Path::new(&args.cxx))?;
    let parent = args
        .output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    let temp = Temporary::new(parent, "aligner-build")?;
    fs::create_dir(temp.path.join("htslib")).map_err(|e| e.to_string())?;
    std::os::unix::fs::symlink(
        prefix.join("include/htslib"),
        temp.path.join("htslib/htslib"),
    )
    .map_err(|e| e.to_string())?;
    let source = include_bytes!("../native/align_alleles.cpp");
    let source_path = temp.path.join("align_alleles.cpp");
    fs::write(&source_path, source).map_err(|e| e.to_string())?;
    let binary = temp.path.join("pathallele-align");
    let command = vec![
        "-std=c++17".into(),
        "-O2".into(),
        format!("-I{}", temp.path.display()),
        format!("-I{}", prefix.join("include/vcflib").display()),
        format!("-I{}", prefix.join("include").display()),
        format!("-I{}", headers.display()),
        source_path.to_string_lossy().into_owned(),
        format!("-L{}", prefix.join("lib").display()),
        format!("-Wl,-rpath,{}", prefix.join("lib").display()),
        "-lvcflib".into(),
        "-lwfa2".into(),
        "-o".into(),
        binary.to_string_lossy().into_owned(),
    ];
    if !Command::new(&compiler)
        .args(&command)
        .status()
        .map_err(|e| e.to_string())?
        .success()
    {
        return Err("Native adapter compilation failed".into());
    }
    fs::hard_link(&binary, &args.output).map_err(|e| e.to_string())?;
    io::json(
        &args.output.with_file_name(format!(
            "{}.build.json",
            args.output.file_name().unwrap().to_string_lossy()
        )),
        &json!({"compiler":compiler,"command":command,"source_sha256":io::hash_bytes(source),"binary_sha256":io::digest(&binary)?}),
    )?;
    println!("{}", args.output.display());
    Ok(())
}
pub fn run() -> Result<()> {
    match Cli::parse().command {
        Commands::Call(a) => call(a),
        Commands::Validate(a) => {
            println!(
                "{}",
                serde_json::to_string_pretty(&check(a)?).map_err(|e| e.to_string())?
            );
            Ok(())
        }
        Commands::Export(a) => export(a),
        Commands::Doctor { aligner } => {
            println!(
                "{}",
                serde_json::to_string_pretty(&runtime(&aligner)?).map_err(|e| e.to_string())?
            );
            Ok(())
        }
        Commands::BuildAligner(a) => build(a),
    }
}

#[cfg(test)]
mod help_tests {
    use super::*;
    use clap::CommandFactory;

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
        Cli::command().debug_assert();
        check(&Cli::command());
    }

    #[test]
    fn short_options_preserve_call_defaults() {
        let parsed = Cli::try_parse_from([
            "pathallele",
            "call",
            "-g",
            "graph.gfa",
            "-r",
            "ref.fa",
            "-s",
            "samples.txt",
            "-o",
            "out",
            "--aligner",
            "pathallele-align",
        ])
        .unwrap();
        let Commands::Call(call) = parsed.command else {
            panic!("expected call")
        };
        assert_eq!(call.graph, PathBuf::from("graph.gfa"));
        assert_eq!(call.settings.threads, 16);
        assert_eq!(call.settings.jobs, 8);
        assert_eq!(call.settings.minimum_anchor, 31);
        assert_eq!(call.settings.mapping_k, 19);
        assert_eq!(call.settings.sequence_anchor, 31);
        assert_eq!(call.settings.minimum_mapq, 20);
    }
}
