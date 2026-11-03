//! Explicit sequence VCF contracts and standard-library-backed VCF reading.
use crate::{
    alleles::Prepared,
    catalogue, classify,
    graph::Graph,
    io,
    model::{self, dna, EditInput, Result},
    project::{Projection, STATUS},
    source::{count, Sources, Stats},
};
use noodles_vcf as vcf;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::Path;
use std::process::{Command, Stdio};

pub const OUTPUT_SCHEMA: u64 = 5;

pub const DESCRIPTIONS: [&str; 11] = [
    "Resolved under this export contract",
    "No unique covering source chain",
    "More than one compatible physical occurrence",
    "Alternative concrete alignments or interacting normalized edits",
    "Non-ACGT source base",
    "Base lies within a non-SNV edit",
    "Edit extends beyond the locus or overlaps another edit",
    "Non-ACGT sequence at a required boundary",
    "Complete locus contains unresolved bases",
    "Complete sequence cannot be represented by this event record",
    "Physical path orientation, order or adjacency is unresolved and no independent complete local source chain is available",
];
pub const CLASS_HEADERS:[&str;7]=[
    "##pathallele_variant_schema=2",
    "##pathallele_hierarchy=\"Unassessed: flat reference projection does not assign nested parents; PID links source contexts in Graph.paths.tsv.gz and Graph.alleles.tsv.gz\"",
    "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Unified allele class: INV (same-interval reverse complement), INV_COMPLEX (local inverted sequence component with other changes), COMPLEX (native structural context), or SEQTYPE; not a unique mutation history\">",
    "##INFO=<ID=SEQTYPE,Number=A,Type=String,Description=\"Sequence shape after trimming common padding: SNV, MNV, INS, DEL or REPLACEMENT\">",
    "##INFO=<ID=NATIVE,Number=A,Type=String,Description=\"Source path orientation class: SEQUENCE, INV, COMPLEX or combined support\">",
    "##INFO=<ID=INVBASIS,Number=A,Type=String,Description=\"Inverted sequence evidence: NONE, EXACT_RC, PATH_AND_RC or LOCAL_RC; coordinates and sequence witnesses in Graph.variant_annotations.tsv.gz\">",
    "##INFO=<ID=SIZECLASS,Number=A,Type=String,Description=\"SV if the larger REF or ALT edit length after trimming common padding is at least 50 bp, otherwise SMALL; independent of TYPE and net length change\">"];
pub fn headers(contract: &str) -> Vec<String> {
    let mut result=vec![format!("##pathallele_output_schema={OUTPUT_SCHEMA}"),format!("##pathallele_genotype_contract={contract}"),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haploid genotype under the declared pathallele_genotype_contract\">".into(),
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample status; definitions in pathallele_status metadata\">".into()];
    result.extend(
        STATUS
            .iter()
            .zip(DESCRIPTIONS)
            .map(|(k, v)| format!("##pathallele_status=<ID={k},Description=\"{v}\">")),
    );
    result
}
pub fn schema(chrom: &str, length: usize, minimum: usize) -> Value {
    json!({"version":OUTPUT_SCHEMA,"coordinate_system":"0-based half-open (VCF uses 1-based POS)","reference":{"contig":chrom,"length":length},"ploidy":1,
    "classification":{"version":2,"annotations":classify::FILE,"minimum_inversion_span":minimum,"inversion_span_contract":"Informative reversed core after removing self reverse-complement-matching flanks; requires a unique seed within that core","sequence_seed_length":minimum.min(31),"inversion_segments":"ref_start0:ref_end0:allele_start0:allele_end0:reference_sequence; pipe-separated witnesses, dot for none; coordinates are exact sequence matches, not inferred historical breakpoints"},
    "event_vcf":{"file":"Graph.vcf.gz","contract":"event_projection","zero":"Reference event backbone; may include nonreference components at other records","links":"Graph.projection_links.tsv.gz"},
    "structural_groups":{"version":1,"genotypes":crate::structure::VCF,"members":crate::structure::MEMBERS,"seed":"Emitted SV ALT alleles","linkage":"Deterministic complete linkage; no phenotype or frequency input","similarity_percent":95,"endpoint_tolerance":"5 percent of shorter span, clamped to 1..50 bp","zero":"Resolved absence at the smallest member locus covering the comparison envelope","scope":"Statistical structural allele groups, not independent ancestral mutations"},
    "complete_loci":{"alleles":"Graph.loci.tsv.gz","genotypes":catalogue::FILE,"index":"Graph.loci.vcf.gz.tbi","authoritative_format":"GT","source_status_format":"FT","missing_value":".","projected_format":"PG","projected_status_format":"PF","allele_index":"GT and PG index the full REF/ALT order in Graph.loci.vcf.gz and Graph.loci.tsv.gz; Graph.vcf.gz ALT indexes may differ","scope":"Independent possibly overlapping locus haplotypes, including reference-only loci; not globally composable records"},
    "statuses":STATUS.iter().zip(DESCRIPTIONS).collect::<BTreeMap<_,_>>(),
    "internal_node_handles":"Graph.nodes.tsv.gz maps signed oriented_node integers to original segment IDs; sign is orientation",
    "path_occurrences":"Graph.paths.tsv.gz: path, query interval, step interval and strand",
    "structural_context":"Graph.structures.tsv.gz: resolved inversions and unresolved adjacency/orientation/order; local calls remain assessable within one independently validated unique source chain, including both locus boundaries"})
}
pub fn classify(reference: &str, alt: &str) -> Result<(&'static str, &'static str)> {
    if reference.is_empty() || alt.is_empty() || !dna(reference) || !dna(alt) {
        return Err("Classification requires explicit nonempty ACGT VCF alleles".into());
    }
    classify::shape(reference, alt)
}

pub struct Record {
    pub chrom: String,
    pub pos: usize,
    pub id: String,
    pub reference: String,
    pub alts: Vec<String>,
    pub info: BTreeMap<String, String>,
    pub gt: Vec<i32>,
    pub ft: Vec<&'static str>,
    pub projected: Option<Vec<i32>>,
    pub projection_status: Option<Vec<&'static str>>,
}
pub struct Reader {
    reader: vcf::io::Reader<Box<dyn BufRead>>,
    pub header: vcf::Header,
    buffer: vcf::Record,
}
impl Reader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut reader = vcf::io::Reader::new(io::reader(path)?);
        let header = reader
            .read_header()
            .map_err(|e| format!("Invalid VCF header: {e}"))?;
        Ok(Self {
            reader,
            header,
            buffer: vcf::Record::default(),
        })
    }
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        if self
            .reader
            .read_record(&mut self.buffer)
            .map_err(|e| e.to_string())?
            == 0
        {
            return Ok(None);
        }
        let r = &self.buffer;
        let pos = usize::from(
            r.variant_start()
                .ok_or("Missing VCF POS")?
                .map_err(|e| e.to_string())?,
        );
        let id = r.ids().as_ref().to_string();
        if id.is_empty() || id.contains(';') {
            return Err("Expected one VCF locus ID".into());
        }
        let mut info = BTreeMap::new();
        let raw = r.info();
        for field in raw.as_ref().split(';') {
            if field.is_empty() {
                continue;
            }
            let (k, v) = field
                .split_once('=')
                .ok_or("Expected valued VCF INFO field")?;
            if info.insert(k.into(), v.into()).is_some() {
                return Err("Duplicate VCF INFO field".into());
            }
        }
        // Header-aware parser checks INFO values, while the restricted public
        // haploid GT:FT contract is decoded directly without per-cell boxing.
        for item in raw.iter(&self.header) {
            item.map_err(|e| e.to_string())?;
        }
        let samples = r.samples();
        let mut fields = samples.as_ref().split('\t');
        let complete = match fields.next() {
            Some("GT:FT") => false,
            Some("GT:FT:PG:PF") => true,
            _ => return Err("Expected haploid GT:FT or GT:FT:PG:PF format".into()),
        };
        let (mut gt, mut ft) = (vec![], vec![]);
        let (mut projected, mut projection_status) = (vec![], vec![]);
        let integer = |g: &str| -> Result<i32> {
            if g == "." {
                return Ok(-1);
            }
            if g.is_empty() || !g.bytes().all(|b| b.is_ascii_digit()) {
                return Err("Non-haploid or invalid genotype".into());
            }
            g.parse::<i32>()
                .map_err(|_| "Genotype integer overflow".into())
        };
        for f in fields {
            let mut cell = f.split(':');
            let value = integer(cell.next().ok_or("Missing GT")?)?;
            let s = cell.next().ok_or("Missing FT")?;
            if s.is_empty() {
                return Err("Invalid GT/FT".into());
            }
            gt.push(value);
            ft.push(
                STATUS
                    .iter()
                    .copied()
                    .find(|&value| value == s)
                    .ok_or("Unknown FT status")?,
            );
            if complete {
                projected.push(integer(cell.next().ok_or("Missing PG")?)?);
                let pf = cell.next().ok_or("Missing PF")?;
                if pf.is_empty() {
                    return Err("Invalid PF".into());
                }
                projection_status.push(
                    STATUS
                        .iter()
                        .copied()
                        .find(|&value| value == pf)
                        .ok_or("Unknown PF status")?,
                );
            }
            if cell.next().is_some() {
                return Err("Unexpected sample subfield".into());
            }
        }
        if gt.len() != self.header.sample_names().len() {
            return Err("VCF sample columns disagree with header".into());
        }
        let alts = r
            .alternate_bases()
            .as_ref()
            .split(',')
            .filter(|s| !s.is_empty() && *s != ".")
            .map(str::to_string)
            .collect();
        Ok(Some(Record {
            chrom: r.reference_sequence_name().into(),
            pos,
            id,
            reference: r.reference_bases().into(),
            alts,
            info,
            gt,
            ft,
            projected: complete.then_some(projected),
            projection_status: complete.then_some(projection_status),
        }))
    }
}

pub fn project(
    prepared: &Prepared,
    sources: &Sources,
    graph: &Graph,
    output: &Path,
    work: &Path,
    minimum: usize,
    threads: usize,
) -> Result<Vec<Vec<EditInput>>> {
    let mut ambiguity = io::tsv(&output.join("Graph.alignment_ambiguity.tsv.gz"))?;
    ambiguity
        .write_record([
            "block_id",
            "allele_index",
            "ambiguous_positions",
            "observed_samples",
        ])
        .map_err(|e| e.to_string())?;
    for (id, ai, n, s) in &prepared.ambiguity_rows {
        ambiguity
            .write_record([id.clone(), ai.to_string(), n.to_string(), s.to_string()])
            .map_err(|e| e.to_string())?;
    }
    io::finish_tsv(ambiguity)?;
    let mut anchors = io::tsv(&output.join("Graph.anchors.tsv.gz"))?;
    anchors
        .write_record([
            "block_id",
            "allele_index",
            "ref_offset0",
            "alt_offset0",
            "length",
            "oriented_node",
            "segment_id",
            "segment_orientation",
        ])
        .map_err(|e| e.to_string())?;
    let mut alleles = io::tsv(&output.join("Graph.alleles.tsv.gz"))?;
    alleles
        .write_record([
            "block_id",
            "allele_index",
            "ref_start0",
            "ref_end0",
            "sequence",
        ])
        .map_err(|e| e.to_string())?;
    for p in &prepared.parents {
        for (ai, (seq, aa)) in p.alleles.iter().zip(&p.anchors).enumerate() {
            alleles
                .write_record([
                    p.id.clone(),
                    ai.to_string(),
                    p.start.to_string(),
                    (p.start + p.alleles[0].len()).to_string(),
                    seq.clone(),
                ])
                .map_err(|e| e.to_string())?;
            for a in aa {
                let (node, name, orientation) = if let Some(n) = a.node {
                    (
                        n.to_string(),
                        graph.nodes[n.unsigned_abs() as usize].name.clone(),
                        if n > 0 { "+" } else { "-" },
                    )
                } else {
                    (String::new(), String::new(), "")
                };
                anchors
                    .write_record([
                        p.id.clone(),
                        ai.to_string(),
                        a.r.to_string(),
                        a.q.to_string(),
                        a.n.to_string(),
                        node,
                        name,
                        orientation.into(),
                    ])
                    .map_err(|e| e.to_string())?;
            }
        }
    }
    io::finish_tsv(anchors)?;
    io::finish_tsv(alleles)?;
    let mut engine = Projection::new(
        graph.reference.clone(),
        sources.samples.clone(),
        prepared
            .models
            .iter()
            .map(|m| m.iter().map(model::Edit::input).collect())
            .collect(),
        prepared
            .known
            .unique
            .iter()
            .map(|r| r.rows.clone())
            .collect(),
        prepared
            .known
            .multiple
            .iter()
            .map(|r| r.rows.clone())
            .collect(),
        prepared.ambiguous.iter().map(|r| r.rows.clone()).collect(),
        sources.structural.clone(),
        sources.chains.clone(),
    )?;
    let mut catalogue =
        catalogue::Writer::new(output, &sources.samples, &graph.chrom, &graph.reference)?;
    let classifier = classify::Classifier::new(&graph.reference, minimum)?;
    let mut annotations = io::tsv(&output.join(classify::FILE))?;
    annotations
        .write_record(classify::COLUMNS)
        .map_err(|e| e.to_string())?;
    let mut f = io::writer(&work.join("Graph.raw.vcf"))?;
    let mut header=vec!["##fileformat=VCFv4.2".into(),format!("##source=PathAllele-{}",env!("CARGO_PKG_VERSION")),format!("##contig=<ID={},length={}>",graph.chrom,graph.reference.len()),
        "##INFO=<ID=PID,Number=.,Type=String,Description=\"Supporting native path source contexts, not nested parent loci\">".into(),
        "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"ALT minus REF length\">".into(),
        "##pathallele_event_reference=\"GT=0 is an event backbone; use Graph.projection_links.tsv.gz and complete-locus gt for full haplotypes\"".into()];
    header.extend(headers("event_projection"));
    header.extend(CLASS_HEADERS.iter().map(|&s| s.into()));
    writeln!(
        f,
        "{}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        header.join("\n"),
        sources.samples.join("\t")
    )
    .map_err(|e| e.to_string())?;
    io::json(
        &output.join("schema.json"),
        &schema(&graph.chrom, graph.reference.len(), minimum),
    )?;
    let mut links = io::writer(&output.join("Graph.projection_links.tsv.gz"))?;
    writeln!(
        links,
        "locus_id\tsample\tcomplete_allele_index\tcomponent_records"
    )
    .map_err(|e| e.to_string())?;
    let mut loci: Vec<_> = prepared
        .origins
        .keys()
        .enumerate()
        .map(|(i, k)| (k, format!("v{i}")))
        .collect();
    loci.sort_by_key(|(k, _)| (k.2 != "SNP", k.1 - k.0, *k));
    let mut seen = HashSet::new();
    let mut counts = Stats::new();
    for ((s, e, kind), id) in loci {
        let row = engine.project(&id, *s, *e, kind == "SNP")?;
        let classes = row
            .alleles
            .iter()
            .skip(1)
            .map(|a| {
                let native = prepared
                    .event_kinds
                    .get(&(*s, *e, a.clone()))
                    .cloned()
                    .unwrap_or_else(|| BTreeSet::from(["SEQUENCE".into()]))
                    .into_iter()
                    .collect::<Vec<_>>()
                    .join("+");
                classifier.annotate(*s, *e, a, &native)
            })
            .collect::<Result<Vec<_>>>()?;
        catalogue.add(&id, *s, *e, kind, &row, &classes)?;
        for (ai, a) in classes.iter().enumerate() {
            annotations
                .write_record(a.fields(&id, ai + 1))
                .map_err(|e| e.to_string())?;
        }
        links
            .write_all(row.links.as_bytes())
            .map_err(|e| e.to_string())?;
        if row.alts.is_empty() {
            continue;
        }
        let (pos, r, alts) =
            model::vcf_alleles(&graph.reference, *s, *e, &row.alts, kind == "SNP")?;
        for alt in &alts {
            if !seen.insert((pos, r.clone(), alt.clone())) {
                return Err("Duplicate projected VCF ALT".into());
            }
        }
        let emitted_classes: Vec<_> = row
            .alts
            .iter()
            .map(|a| {
                let index = row
                    .alleles
                    .iter()
                    .position(|s| s == a)
                    .expect("projected ALT is a catalogue allele");
                &classes[index - 1]
            })
            .collect();
        let pid = prepared.origins[&(*s, *e, kind.clone())]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
            .join(",");
        let svlen = alts
            .iter()
            .map(|a| (a.len() as i64 - r.len() as i64).to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            f,
            "{}\t{pos}\t{id}\t{r}\t{}\t.\t.\tPID={pid};SVLEN={svlen};{}\tGT:FT\t{}",
            graph.chrom,
            alts.join(","),
            classify::info(&emitted_classes),
            row.cells
        )
        .map_err(|e| e.to_string())?;
        count(&mut counts, "records", 1);
        count(&mut counts, "ALT_alleles", alts.len());
        for (i, &c) in row.status_counts.iter().enumerate() {
            if c > 0 {
                count(&mut counts, format!("GT_{}", STATUS[i]), c);
            }
        }
    }
    f.finish()?;
    links.finish()?;
    catalogue.finish(threads)?;
    io::finish_tsv(annotations)?;
    io::json(
        &output.join("summary.projection.json"),
        &json!({"method":"native_path_projection","samples":sources.samples.len(),"counts":counts}),
    )?;
    Ok(engine.positives())
}

pub fn finish_vcf(output: &Path, reference: &Path, work: &Path) -> Result<()> {
    let sort = catalogue::Temporary::new(work, "sort")?;
    let specifications: Vec<Vec<std::ffi::OsString>> = vec![
        vec![
            "norm".into(),
            "-N".into(),
            "-f".into(),
            reference.into(),
            "-c".into(),
            "e".into(),
            work.join("projection/Graph.raw.vcf").into(),
            "-Ou".into(),
        ],
        vec![
            "sort".into(),
            "-T".into(),
            sort.path.join("chunks").into(),
            "-Ou".into(),
        ],
        vec![
            "+fill-tags".into(),
            "-Ou".into(),
            "--".into(),
            "-t".into(),
            "AC,AN,AF".into(),
        ],
        vec![
            "view".into(),
            "-c".into(),
            "1".into(),
            "-Oz".into(),
            "-o".into(),
            output.join("Graph.vcf.gz").into(),
        ],
    ];
    let mut children = vec![];
    let mut previous = None;
    for (i, args) in specifications.iter().enumerate() {
        let mut command = Command::new("bcftools");
        command.args(args).stderr(Stdio::from(
            File::create(work.join(format!("vcf.{i}.log"))).map_err(|e| e.to_string())?,
        ));
        if let Some(pipe) = previous.take() {
            command.stdin(Stdio::from(pipe));
        }
        if i < specifications.len() - 1 {
            command.stdout(Stdio::piped());
        }
        match command.spawn() {
            Ok(mut child) => {
                previous = child.stdout.take();
                children.push(child);
            }
            Err(e) => {
                for c in &mut children {
                    let _ = c.kill();
                    let _ = c.wait();
                }
                return Err(e.to_string());
            }
        }
    }
    let mut failures = vec![];
    for (i, mut child) in children.into_iter().enumerate() {
        if !child.wait().map_err(|e| e.to_string())?.success() {
            failures.push(i);
        }
    }
    if !failures.is_empty() {
        return Err(format!(
            "VCF pipeline failed at stages {failures:?}; inspect work/vcf.*.log"
        ));
    }
    for file in ["Graph.vcf.gz", catalogue::FILE] {
        let status = Command::new("bcftools")
            .args(["index", "-ft"])
            .arg(output.join(file))
            .status()
            .map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("VCF indexing failed: {file}"));
        }
    }
    Ok(())
}
