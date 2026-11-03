//! The source stage owns physical paths and emits bounded, typed evidence.
use crate::{
    anchors::{self, Anchor, MatchIndex, SeedCounts},
    context::{self, Block, Oriented, Plan},
    graph::{Graph, Walk},
    io,
    model::{dna, Intervals, Ranges, Result},
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

pub type Stats = BTreeMap<String, u64>;
pub fn count(stats: &mut Stats, key: impl Into<String>, value: usize) {
    *stats.entry(key.into()).or_default() += value as u64;
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Parent {
    pub id: String,
    pub start: usize,
    pub alleles: Vec<String>,
    pub carriers: Vec<Vec<usize>>,
    pub anchors: Vec<Vec<Anchor>>,
    pub kinds: Vec<String>,
}
#[derive(Clone, Serialize, Deserialize)]
pub struct Sources {
    pub samples: Vec<String>,
    pub parents: Vec<Parent>,
    pub coverage: Vec<Intervals>,
    pub uncertain: Vec<Intervals>,
    pub chains: Vec<Intervals>,
    pub structural: Vec<Intervals>,
    pub stats: Stats,
}
#[derive(Clone, Serialize, Deserialize)]
pub struct SourceConfig {
    pub threads: usize,
    pub minimum_anchor: usize,
    pub mapping_k: usize,
    pub maximum_gap: usize,
    pub minimum_mapq: u32,
    pub mapping_batch_bases: usize,
    pub mapping_cache_bytes: usize,
}
type WindowKey = (usize, usize, String, Vec<(usize, usize, usize)>, String);
type Tsv = csv::Writer<io::Output>;

struct SourceInterval<'a> {
    query: (usize, usize),
    reference: (usize, usize),
    allele: (&'a str, usize),
    status: &'a str,
}
fn emit(
    writer: &mut Tsv,
    path: &Walk,
    sample: &str,
    view: &Oriented,
    block: &Block,
    row: SourceInterval<'_>,
) -> Result<()> {
    let (qs, qe) = row.query;
    let (rs, re) = row.reference;
    let (id, allele) = row.allele;
    let status = row.status;
    let name = &path.name;
    let length = path.length;
    let reverse = block.reverse;
    let component = &block.identity;
    let kind = &block.kind;
    let mut start = view.offsets.partition_point(|&v| v <= qs).saturating_sub(1);
    let mut end = view.offsets.partition_point(|&v| v < qe);
    let (mut qs, mut qe) = (qs, qe);
    if reverse {
        (qs, qe) = (length - qe, length - qs);
        (start, end) = (view.walk.len() - end, view.walk.len() - start);
    }
    writer
        .write_record([
            sample.to_string(),
            name.to_string(),
            qs.to_string(),
            qe.to_string(),
            rs.to_string(),
            re.to_string(),
            id.to_string(),
            allele.to_string(),
            start.to_string(),
            end.to_string(),
            status.to_string(),
            if reverse { "-" } else { "+" }.to_string(),
            component.to_string(),
            kind.to_string(),
        ])
        .map_err(|e| e.to_string())
}

struct Collector<'a> {
    graph: &'a Graph,
    config: &'a SourceConfig,
    parents: Vec<Parent>,
    windows: HashMap<WindowKey, usize>,
    coverage: Vec<Intervals>,
    chains: Vec<Intervals>,
    structural: Vec<Intervals>,
    stats: Stats,
    paths: Tsv,
    structures: Tsv,
}
impl Collector<'_> {
    fn accept(&mut self, path: usize, plan: Plan) -> Result<()> {
        let graph = self.graph;
        let pi = &graph.paths[path];
        let si = pi
            .sample
            .ok_or("Reference path cannot be a sample source")?;
        let sample = &graph.samples[si];
        let name = &pi.name;
        let length = plan.oriented[0].query.len();
        let (evidence, unresolved) =
            context::structure_evidence(&plan.blocks, length, graph.reference.len())?;
        self.structural[si].extend(unresolved);
        for (qs, qe, rs, re, status, components) in evidence {
            self.structures
                .write_record([
                    sample.clone(),
                    name.clone(),
                    qs.to_string(),
                    qe.to_string(),
                    rs.to_string(),
                    re.to_string(),
                    status.clone(),
                    components,
                ])
                .map_err(|e| e.to_string())?;
            count(&mut self.stats, format!("structure_{status}"), 1);
        }
        count(&mut self.stats, "paths", 1);
        count(&mut self.stats, "native_bp", length);
        for (lo, hi) in plan.missing {
            let offsets = &plan.oriented[0].offsets;
            self.paths
                .write_record([
                    sample.clone(),
                    name.clone(),
                    lo.to_string(),
                    hi.to_string(),
                    ".".into(),
                    ".".into(),
                    ".".into(),
                    "-1".into(),
                    offsets
                        .partition_point(|&v| v <= lo)
                        .saturating_sub(1)
                        .to_string(),
                    offsets.partition_point(|&v| v < hi).to_string(),
                    "UNBOUNDED".into(),
                    "+".into(),
                    ".".into(),
                    "unresolved".into(),
                ])
                .map_err(|e| e.to_string())?;
            count(&mut self.stats, "unplaced_bp", hi - lo);
        }
        for block in plan.blocks {
            let view = &plan.oriented[block.reverse as usize];
            let runs = &block.runs;
            count(&mut self.stats, format!("{}_components", block.kind), 1);
            for a in runs {
                if graph.reference[a.r..a.r + a.n] != view.query[a.q..a.q + a.n] {
                    return Err(format!("Source match differs from reference: {name}"));
                }
                self.coverage[si].push((a.r, a.r + a.n));
                count(&mut self.stats, "context_matched_bp", a.n);
                emit(
                    &mut self.paths,
                    pi,
                    sample,
                    view,
                    &block,
                    SourceInterval {
                        query: (a.q, a.q + a.n),
                        reference: (a.r, a.r + a.n),
                        allele: (".", 0),
                        status: "MATCH",
                    },
                )?;
            }
            let last = runs.last().ok_or("Empty source component")?;
            self.chains[si].push((runs[0].r, last.r + last.n));
            let index = block
                .record
                .as_ref()
                .map(|r| r.index.clone())
                .unwrap_or_else(|| {
                    MatchIndex::new(runs.iter().map(|r| (r.q, r.q + r.n, r.r)).collect())
                });
            for pair in runs.windows(2) {
                let left = &pair[0];
                let right = &pair[1];
                let (rs, re) = (left.r + left.n, right.r);
                let (qs, qe) = (left.q + left.n, right.q);
                if rs > re || qs > qe {
                    return Err("Source component is not collinear".into());
                }
                self.coverage[si].push((rs, re));
                if graph.reference[rs..re] == view.query[qs..qe] {
                    if qe > qs {
                        emit(
                            &mut self.paths,
                            pi,
                            sample,
                            view,
                            &block,
                            SourceInterval {
                                query: (qs, qe),
                                reference: (rs, re),
                                allele: (".", 0),
                                status: "MATCH",
                            },
                        )?;
                        count(&mut self.stats, "gap_matched_bp", qe - qs);
                    }
                    continue;
                }
                if rs == 0 || qs == 0 || re >= graph.reference.len() || qe >= view.query.len() {
                    return Err("Source allele lacks two flanks".into());
                }
                let reference = &graph.reference[rs - 1..re + 1];
                let alternate = &view.query[qs - 1..qe + 1];
                if reference.as_bytes()[0] != alternate.as_bytes()[0]
                    || reference.as_bytes().last() != alternate.as_bytes().last()
                {
                    return Err("Source allele lacks exact flanks".into());
                }
                let rp = &graph.ref_offsets;
                let qp = &view.offsets;
                let ra = rp.partition_point(|&v| v < rs);
                let rb = rp.partition_point(|&v| v <= re) - 1;
                let qa = qp.partition_point(|&v| v < qs);
                let qb = qp.partition_point(|&v| v <= qe) - 1;
                let ri = if ra <= rb {
                    &graph.refwalk()[ra..rb]
                } else {
                    &[]
                };
                let qi = if qa <= qb { &view.walk[qa..qb] } else { &[] };
                let rnodes: HashSet<_> = ri.iter().copied().collect();
                let mut back = 0;
                let mut forward = 0;
                for n in qi {
                    let len = graph.nodes[n.unsigned_abs() as usize].sequence.len();
                    if rnodes.contains(n) {
                        forward += len;
                    } else if rnodes.contains(&-n) {
                        back += len;
                    }
                }
                let complete = rp[ra] == rs && rp[rb] == re && qp[qa] == qs && qp[qb] == qe;
                let inverted = complete
                    && !ri.is_empty()
                    && ri.len() == qi.len()
                    && ri.iter().rev().zip(qi).all(|(&a, &b)| a == -b);
                let kind = block.native_kind.clone().unwrap_or_else(|| {
                    if inverted {
                        "INV"
                    } else if back >= self.config.minimum_anchor && back > forward {
                        "COMPLEX"
                    } else {
                        "SEQUENCE"
                    }
                    .into()
                });
                let anchors = if kind != "SEQUENCE" {
                    vec![]
                } else {
                    anchors::bounded(
                        reference,
                        alternate,
                        ri,
                        qi,
                        &graph.nodes,
                        rp[ra] - rs + 1,
                        qp[qa] - qs + 1,
                        &index,
                        rs - 1,
                        qs - 1,
                    )?
                };
                let key = (
                    rs - 1,
                    re + 1,
                    alternate.to_string(),
                    anchors.iter().map(|a| (a.r, a.q, a.n)).collect(),
                    kind.clone(),
                );
                let parent_index = if let Some(&i) = self.windows.get(&key) {
                    i
                } else {
                    let i = self.parents.len();
                    let id = format!("PC{i:06}");
                    self.windows.insert(key, i);
                    self.parents.push(Parent {
                        id,
                        start: rs - 1,
                        alleles: vec![reference.into(), alternate.into()],
                        carriers: vec![vec![], vec![]],
                        anchors: vec![vec![], anchors],
                        kinds: vec!["REFERENCE".into(), kind.clone()],
                    });
                    i
                };
                let parent = &mut self.parents[parent_index];
                parent.carriers[1].push(si);
                emit(
                    &mut self.paths,
                    pi,
                    sample,
                    view,
                    &block,
                    SourceInterval {
                        query: (qs, qe),
                        reference: (rs, re),
                        allele: (&parent.id, 1),
                        status: if dna(alternate) {
                            "ALLELE"
                        } else {
                            "PARTIAL_SEQUENCE"
                        },
                    },
                )?;
                count(&mut self.stats, "allele_bp", qe - qs);
                count(&mut self.stats, "allele_observations", 1);
                count(&mut self.stats, format!("native_{kind}"), 1);
                count(
                    &mut self.stats,
                    format!("{}_allele_observations", block.kind),
                    1,
                );
            }
        }
        Ok(())
    }
}

pub fn build(graph: &Graph, work: &Path, output: &Path, config: &SourceConfig) -> Result<Sources> {
    let mut fasta = io::writer(&work.join("paths.fa"))?;
    for paths in &graph.sample_paths {
        for &path in paths {
            writeln!(fasta, ">{}\n{}", graph.paths[path].name, graph.spell(path))
                .map_err(|e| e.to_string())?;
        }
    }
    fasta.finish()?;
    let status = Command::new("minimap2")
        .args([
            "-cx",
            "asm20",
            "-k",
            &config.mapping_k.to_string(),
            "-s",
            "40",
            "--end-bonus",
            "100",
            "--no-hash-name",
            "--eqx",
            "--secondary=yes",
            "-N",
            "20",
            "-K",
            &config.mapping_batch_bases.to_string(),
            "--cap-kalloc",
            &config.mapping_cache_bytes.to_string(),
            "-t",
            &config.threads.to_string(),
        ])
        .arg(work.join("reference.fa"))
        .arg(work.join("paths.fa"))
        .stdout(Stdio::from(
            File::create(work.join("context.paf")).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            File::create(work.join("context.log")).map_err(|e| e.to_string())?,
        ))
        .status()
        .map_err(|e| e.to_string())?;
    if !status.success() {
        return Err("minimap2 failed; inspect work/context.log".into());
    }
    let context = context::read_alignments(&work.join("context.paf"), config.minimum_mapq)?;
    let uncertain = context::copy_intervals(graph, &context)?;
    let n = graph.samples.len();
    let mut paths = io::tsv(&output.join("Graph.paths.tsv.gz"))?;
    paths
        .write_record([
            "sample",
            "path",
            "query_start0",
            "query_end0",
            "ref_start0",
            "ref_end0",
            "block_id",
            "allele_index",
            "start_step0",
            "end_step0",
            "status",
            "strand",
            "component",
            "context_kind",
        ])
        .map_err(|e| e.to_string())?;
    let mut structures = io::tsv(&output.join("Graph.structures.tsv.gz"))?;
    structures
        .write_record([
            "sample",
            "path",
            "query_start0",
            "query_end0",
            "ref_start0",
            "ref_end0",
            "status",
            "components",
        ])
        .map_err(|e| e.to_string())?;
    let mut stats = Stats::new();
    count(
        &mut stats,
        "competing_physical_copy_bp",
        uncertain.iter().flatten().map(|(s, e)| e - s).sum(),
    );
    let mut collector = Collector {
        graph,
        config,
        parents: vec![],
        windows: HashMap::new(),
        coverage: vec![vec![]; n],
        chains: vec![vec![]; n],
        structural: vec![vec![]; n],
        stats,
        paths,
        structures,
    };
    let mut ks = vec![config.minimum_anchor, config.mapping_k];
    ks.sort_unstable();
    ks.dedup();
    let counts: Vec<_> = ks
        .into_iter()
        .map(|k| SeedCounts::new(&graph.reference, k))
        .collect();
    let names: Vec<_> = graph.sample_paths.iter().flatten().copied().collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads.min(names.len()).max(1))
        .build()
        .map_err(|e| e.to_string())?;
    for chunk in names.chunks(config.threads.max(1) * 4) {
        let planned: Vec<_> = pool.install(|| {
            chunk
                .par_iter()
                .map(|&path| {
                    context::plan(
                        graph,
                        path,
                        context
                            .get(&graph.paths[path].name)
                            .map(Vec::as_slice)
                            .unwrap_or(&[]),
                        config.minimum_anchor,
                        config.mapping_k,
                        config.maximum_gap,
                        &counts,
                    )
                })
                .collect()
        });
        for (&path, plan) in chunk.iter().zip(planned) {
            collector.accept(path, plan?)?;
        }
    }
    for parent in &mut collector.parents {
        for carriers in &mut parent.carriers {
            carriers.sort_unstable();
            carriers.dedup();
        }
    }
    for rows in &mut collector.structural {
        *rows = Ranges::new(std::mem::take(rows), graph.reference.len())?.rows;
    }
    let accounted: u64 = [
        "context_matched_bp",
        "gap_matched_bp",
        "allele_bp",
        "unplaced_bp",
    ]
    .iter()
    .map(|&k| collector.stats.get(k).copied().unwrap_or(0))
    .sum();
    if collector.stats.get("native_bp").copied().unwrap_or(0) != accounted {
        return Err("Native bases do not partition into source states".into());
    }
    count(&mut collector.stats, "blocks", collector.parents.len());
    io::finish_tsv(collector.paths)?;
    io::finish_tsv(collector.structures)?;
    let sources = Sources {
        samples: graph.samples.clone(),
        parents: collector.parents,
        coverage: collector.coverage,
        uncertain,
        chains: collector.chains,
        structural: collector.structural,
        stats: collector.stats,
    };
    io::json(&output.join("summary.context.json"), &sources.stats)?;
    io::checkpoint(&work.join("sources.bin"), &sources)?;
    Ok(sources)
}
