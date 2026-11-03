//! Validation is independent of genotype inference; no uncertain call is imputed.
use crate::{
    alleles::{self, Prepared},
    catalogue::Catalogue,
    classify,
    graph::Graph,
    io,
    model::{self, dna, intersects, replay, Edit, EditInput, Evidence, Ranges, Result},
    source::{count, Sources, Stats},
    variant,
};
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::Path;
type Key = (usize, usize, String);
fn require(ok: bool, message: &str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(format!("Validation failed: {message}"))
    }
}
fn table(path: &Path, columns: &[&str]) -> Result<csv::Reader<Box<dyn std::io::BufRead>>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(io::reader(path)?);
    require(
        reader
            .headers()
            .map_err(|e| e.to_string())?
            .iter()
            .eq(columns.iter().copied()),
        "table header",
    )?;
    Ok(reader)
}
fn integer(s: &str) -> Result<usize> {
    s.parse().map_err(|_| format!("Invalid table integer: {s}"))
}

fn source_checks(
    output: &Path,
    graph: &Graph,
    sources: &Sources,
    prepared: &Prepared,
) -> Result<Stats> {
    require(sources.samples == graph.samples, "source cohort order")?;
    let mut stats = alleles::validate_models(
        &sources.parents,
        &prepared.models,
        &graph.reference,
        &prepared.known,
        &prepared.ambiguous,
    )?;
    require(
        stats == prepared.model_validation,
        "model replay summary differs",
    )?;
    let parents: HashMap<_, _> = sources
        .parents
        .iter()
        .enumerate()
        .map(|(i, p)| (p.id.as_str(), i))
        .collect();
    let paths: HashMap<_, _> = graph
        .paths
        .iter()
        .enumerate()
        .map(|(i, p)| (p.name.as_str(), i))
        .collect();
    let samples: HashMap<_, _> = graph
        .samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut observed = vec![BTreeSet::new(); sources.parents.len()];
    let mut partitions = vec![vec![]; graph.paths.len()];
    let mut current = None;
    let (mut query, mut reverse) = (String::new(), String::new());
    let mut rows = table(
        &output.join("Graph.paths.tsv.gz"),
        &[
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
        ],
    )?;
    for row in rows.records() {
        let row = row.map_err(|e| e.to_string())?;
        let &pi = paths.get(&row[1]).ok_or("Unknown source path")?;
        let &si = samples.get(&row[0]).ok_or("Unknown source sample")?;
        require(graph.paths[pi].sample == Some(si), "path/sample ownership")?;
        if current != Some(pi) {
            query = graph.spell(pi);
            reverse = io::reverse_complement(&query);
            current = Some(pi);
        }
        let (lo, hi) = (integer(&row[2])?, integer(&row[3])?);
        require(lo <= hi && hi <= query.len(), "physical query bounds")?;
        let (a, b) = (integer(&row[8])?, integer(&row[9])?);
        require(
            a <= b && b <= graph.paths[pi].nodes.len(),
            "physical step bounds",
        )?;
        if lo < hi {
            partitions[pi].push((lo, hi));
        }
        count(&mut stats, format!("source_{}", &row[10]), 1);
        if &row[10] == "UNBOUNDED" {
            continue;
        }
        let (rs, re) = (integer(&row[4])?, integer(&row[5])?);
        require(
            rs <= re && re <= graph.reference.len(),
            "physical reference bounds",
        )?;
        let (sequence, qs, qe) = match &row[11] {
            "+" => (query.as_str(), lo, hi),
            "-" => (reverse.as_str(), query.len() - hi, query.len() - lo),
            _ => return Err("Invalid source strand".into()),
        };
        if &row[10] == "MATCH" {
            require(
                sequence[qs..qe] == graph.reference[rs..re],
                "source exact match",
            )?;
        } else {
            require(
                matches!(&row[10], "ALLELE" | "PARTIAL_SEQUENCE"),
                "unknown source status",
            )?;
            let &i = parents.get(&row[6]).ok_or("Unknown source parent")?;
            let p = &sources.parents[i];
            require(
                rs > 0 && re < graph.reference.len() && qs > 0 && qe < sequence.len(),
                "source flanking bases",
            )?;
            require(
                p.start == rs - 1 && p.alleles[0] == graph.reference[rs - 1..re + 1],
                "source REF sequence",
            )?;
            require(
                integer(&row[7])? == 1 && p.alleles[1] == sequence[qs - 1..qe + 1],
                "source ALT sequence",
            )?;
            observed[i].insert(si);
        }
    }
    for paths in &graph.sample_paths {
        for &pi in paths {
            let intervals = &mut partitions[pi];
            intervals.sort_unstable();
            let mut cursor = 0;
            for &(lo, hi) in intervals.iter() {
                require(lo == cursor, "native bases overlap or disappear")?;
                cursor = hi;
            }
            require(cursor == graph.paths[pi].length, "unaccounted path end")?;
            count(&mut stats, "fully_partitioned_paths", 1);
        }
    }
    require(
        sources.parents.len() == prepared.edits.len(),
        "normalized parent dimensions",
    )?;
    for (i, p) in sources.parents.iter().enumerate() {
        require(
            p.carriers[1].iter().copied().collect::<BTreeSet<_>>() == observed[i],
            "source carrier identity",
        )?;
        require(
            p.alleles.len() == prepared.edits[i].len(),
            "normalized allele dimensions",
        )?;
        for (a, es) in p.alleles.iter().zip(&prepared.edits[i]) {
            require(
                replay(
                    &p.alleles[0],
                    0,
                    p.alleles[0].len(),
                    &es.iter().collect::<Vec<_>>(),
                )? == *a,
                "normalized source allele replay",
            )?;
            count(&mut stats, "replayed_source_alleles", 1);
        }
    }
    Ok(stats)
}

pub struct FullContext<'a> {
    pub reference: &'a str,
    pub prepared: &'a Prepared,
    pub structural: Vec<Ranges>,
    pub chains: &'a [crate::model::Intervals],
}
pub struct EventCheck {
    pub stats: Stats,
    pub portable: Stats,
    pub emitted: Vec<BTreeSet<Key>>,
}
pub fn events(
    output: &Path,
    data: &Catalogue,
    chrom: &str,
    length: usize,
    full: Option<&FullContext>,
) -> Result<EventCheck> {
    let annotations = classify::Table::read(output, data, full.map(|c| c.reference))?;
    if let Some(c) = full {
        for (li, seqs) in data.alleles.iter().enumerate() {
            let (s, e) = data.bounds[li];
            for (ai, a) in seqs.iter().enumerate().skip(1) {
                let native = c
                    .prepared
                    .event_kinds
                    .get(&(s, e, a.clone()))
                    .cloned()
                    .unwrap_or_else(|| BTreeSet::from(["SEQUENCE".into()]))
                    .into_iter()
                    .collect::<Vec<_>>()
                    .join("+");
                require(
                    annotations.get(li, ai).native == native,
                    "native annotation provenance",
                )?;
            }
        }
    }
    let sample_index: HashMap<_, _> = data
        .samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut portable = Stats::from([
        ("loci".into(), data.loci.len() as u64),
        ("samples".into(), data.samples.len() as u64),
    ]);
    let mut stats = Stats::new();
    let mut reader = variant::Reader::open(&output.join("Graph.vcf.gz"))?;
    require(
        reader.header.sample_names().iter().eq(data.samples.iter()),
        "VCF sample order",
    )?;
    require(
        reader
            .header
            .contigs()
            .keys()
            .map(String::as_str)
            .eq([chrom]),
        "VCF contig",
    )?;
    let mut metadata = vec![];
    noodles_vcf::io::Writer::new(&mut metadata)
        .write_header(&reader.header)
        .map_err(|e| e.to_string())?;
    let metadata = String::from_utf8(metadata).map_err(|e| e.to_string())?;
    for required in [
        "##pathallele_genotype_contract=event_projection\n",
        format!(
            "##pathallele_output_schema={}\n",
            crate::variant::OUTPUT_SCHEMA
        )
        .as_str(),
        "##pathallele_variant_schema=2\n",
    ] {
        require(metadata.contains(required), "event contract/schema header")?;
    }
    for key in ["TYPE", "SEQTYPE", "SIZECLASS", "NATIVE", "INVBASIS"] {
        require(
            reader.header.infos().get(key).is_some_and(|v| {
                v.number() == noodles_vcf::header::record::value::map::info::Number::A
            }),
            "classification must be per ALT",
        )?;
    }
    let mut emitted = vec![BTreeSet::new(); data.samples.len()];
    // Store each VCF ALT once; sample cells were checked against the mapped
    // catalogue below. Links are streamed after all VCF witnesses are known.
    let mut components: Vec<Vec<(i32, Edit)>> = vec![vec![]; data.loci.len()];
    let mut seen = HashSet::new();
    let mut alleles_seen = HashSet::new();
    let mut last = 0;
    while let Some(record) = reader.next_record()? {
        let &li = data.index.get(&record.id).ok_or("Unknown VCF locus")?;
        require(seen.insert(record.id.clone()), "duplicate VCF ID")?;
        require(
            record.chrom == chrom && record.pos >= last,
            "VCF contig/order",
        )?;
        last = record.pos;
        let (s, e) = data.bounds[li];
        let seqs = &data.alleles[li];
        let mut alleles = vec![record.reference.as_str()];
        alleles.extend(record.alts.iter().map(String::as_str));
        let unpadded: Vec<_> = if data.kinds[li] == "SNP" {
            require(record.pos == s + 1, "SNP coordinate")?;
            alleles.clone()
        } else {
            require(
                record.pos == s.max(1) && alleles.iter().all(|a| !a.is_empty()),
                "event padding/coordinate",
            )?;
            alleles
                .iter()
                .map(|a| if s > 0 { &a[1..] } else { &a[..a.len() - 1] })
                .collect()
        };
        require(unpadded[0] == seqs[0], "VCF locus REF")?;
        if let Some(c) = full {
            require(
                c.reference
                    .get(record.pos - 1..record.pos - 1 + record.reference.len())
                    == Some(record.reference.as_str()),
                "VCF reference sequence",
            )?;
        }
        let classes: Vec<_> = unpadded
            .iter()
            .skip(1)
            .map(|a| {
                let ai = seqs
                    .iter()
                    .position(|s| s == a)
                    .ok_or("VCF ALT absent from complete catalogue")?;
                if ai == 0 {
                    return Err("VCF ALT equals REF".into());
                }
                Ok(annotations.get(li, ai))
            })
            .collect::<Result<_>>()?;
        for field in classify::info(&classes).split(';') {
            let (k, v) = field.split_once('=').unwrap();
            require(
                record.info.get(k).map(String::as_str) == Some(v),
                "allele classification mismatch",
            )?;
        }
        let mut edits = vec![];
        for alt in &record.alts {
            require(
                alleles_seen.insert((record.pos, record.reference.clone(), alt.clone())),
                "duplicate ALT",
            )?;
            edits.push(model::trim(record.pos - 1, &record.reference, alt));
        }
        for (ai, edit) in edits.iter().enumerate() {
            let index = seqs
                .iter()
                .position(|s| s == unpadded[ai + 1])
                .ok_or("VCF ALT absent from complete catalogue")?;
            components[li].push((index as i32, edit.clone()));
        }
        let mut ac = vec![0usize; record.alts.len()];
        let mut an = 0;
        for (si, emitted_row) in emitted.iter_mut().enumerate() {
            let g = record.gt[si];
            let pg = data.projected(li, si);
            require((g < 0) == (pg < 0), "VCF projection missingness")?;
            require(record.ft[si] == data.status(li, si, true), "VCF status")?;
            if g < 0 {
                count(&mut stats, format!("missing_{}", record.ft[si]), 1);
                continue;
            }
            let ai = g as usize;
            require(
                ai < unpadded.len() && unpadded[ai] == seqs[pg as usize],
                "VCF projected allele",
            )?;
            an += 1;
            if let Some(c) = full {
                let (lo, hi) = if s < e {
                    (s, e)
                } else {
                    (s.saturating_sub(1), (e + 1).min(length))
                };
                require(
                    !c.structural[si].overlaps(lo, hi)
                        || c.chains[si]
                            .iter()
                            .any(|&(a, b)| a <= s.saturating_sub(1) && (e + 1).min(length) <= b),
                    "structural context lacks an independent complete source chain",
                )?;
                require(
                    c.prepared.known.unique[si].covers(lo, hi),
                    "confident GT in unresolved coverage",
                )?;
            }
            if ai > 0 {
                let edit = &edits[ai - 1];
                let key = (edit.start, edit.end, edit.alt.clone());
                require(emitted_row.insert(key), "duplicate ALT genotype")?;
                ac[ai - 1] += 1;
                if let Some(c) = full {
                    let mut x = edit.clone();
                    let found = c.prepared.models[si].contains(&x);
                    x.complex = true;
                    require(
                        found || c.prepared.models[si].contains(&x),
                        "ALT has no sequence-derived model witness",
                    )?;
                }
                count(&mut stats, "source_supported_ALT_genotypes", 1);
            } else {
                count(&mut stats, "reference_genotypes", 1);
            }
        }
        let expected_ac = ac
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        require(
            record.info.get("AC") == Some(&expected_ac)
                && record.info.get("AN") == Some(&an.to_string()),
            "VCF AC/AN",
        )?;
        require(ac.iter().all(|&v| v > 0), "unobserved ALT in VCF")?;
        count(&mut stats, "vcf_records", 1);
        count(&mut stats, "vcf_ALT_alleles", record.alts.len());
        count(&mut portable, "event_records", 1);
    }
    for (li, id) in data.loci.iter().enumerate() {
        require(
            seen.contains(id) || (0..data.samples.len()).all(|si| data.projected(li, si) <= 0),
            "missing positive event record",
        )?;
    }
    let mut seen_links = HashSet::new();
    let mut rows = table(
        &output.join("Graph.projection_links.tsv.gz"),
        &[
            "locus_id",
            "sample",
            "complete_allele_index",
            "component_records",
        ],
    )?;
    for row in rows.records() {
        let r = row.map_err(|e| e.to_string())?;
        let &li = data
            .index
            .get(&r[0])
            .ok_or("Unknown projection link locus")?;
        let &si = sample_index
            .get(&r[1])
            .ok_or("Unknown projection link sample")?;
        require(seen_links.insert((li, si)), "duplicate projection link")?;
        let ai = integer(&r[2])?;
        require(
            ai > 0 && data.gt(li, si) == ai as i32 && data.projected(li, si) == 0,
            "incorrect reference-backbone link",
        )?;
        let ids: Vec<_> = r[3].split(',').collect();
        require(
            !ids.is_empty()
                && !ids.iter().any(|s| s.is_empty())
                && ids.iter().collect::<HashSet<_>>().len() == ids.len(),
            "component identities",
        )?;
        let (s, e) = data.bounds[li];
        let mut local = vec![];
        for id in ids {
            let &ci = data.index.get(id).ok_or("Unknown component locus")?;
            let pg = data.projected(ci, si);
            require(pg > 0, "component has no positive genotype")?;
            let x = &components[ci]
                .iter()
                .find(|(ai, _)| *ai == pg)
                .ok_or("Absent linked VCF component")?
                .1;
            require(x.start >= s && x.end <= e, "out-of-locus component")?;
            local.push(Edit {
                start: x.start - s,
                end: x.end - s,
                alt: x.alt.clone(),
                complex: false,
            });
        }
        let reference = &data.alleles[li][0];
        require(
            replay(
                reference,
                0,
                reference.len(),
                &local.iter().collect::<Vec<_>>(),
            )? == data.alleles[li][data.gt(li, si) as usize],
            "component replay changes complete allele",
        )?;
    }
    let expected_links = (0..data.loci.len())
        .map(|li| {
            (0..data.samples.len())
                .filter(|&si| data.gt(li, si) > 0 && data.projected(li, si) == 0)
                .count()
        })
        .sum::<usize>();
    require(
        expected_links == seen_links.len(),
        "incomplete component-link table",
    )?;
    portable.insert("component_links".into(), seen_links.len() as u64);
    Ok(EventCheck {
        stats,
        portable,
        emitted,
    })
}

pub fn portable(output: &Path, data: &Catalogue, chrom: &str, length: usize) -> Result<Value> {
    let result = events(output, data, chrom, length, None)?;
    let structure = crate::structure::Plan::new(output, data)?;
    structure.validate(output, data, chrom)?;
    Ok(
        json!({"status":"pass","schema":crate::variant::OUTPUT_SCHEMA,"classification_schema":2,"counts":result.portable,"structural_groups":structure.summary(),"scope":"Portable output hashes, complete-locus VCF and sequence catalogue, unified allele annotations and inversion witness replay, haploid event genotypes, complete-locus component reconstruction and structural grouping/genotypes","limitation":"No original graph/source validation or genome-wide inversion uniqueness recheck, and no independent biological accuracy claim"}),
    )
}

pub fn full(
    output: &Path,
    graph: &Graph,
    sources: &Sources,
    prepared: &Prepared,
    positive: &[Vec<EditInput>],
    data: &Catalogue,
) -> Result<Value> {
    require(
        data.samples == graph.samples,
        "joint catalogue cohort order",
    )?;
    let mut stats = source_checks(output, graph, sources, prepared)?;
    let structural: Vec<_> = sources
        .structural
        .iter()
        .map(|r| Ranges::new(r.clone(), graph.reference.len()))
        .collect::<Result<_>>()?;
    let context = FullContext {
        reference: &graph.reference,
        prepared,
        structural,
        chains: &sources.chains,
    };
    let result = events(
        output,
        data,
        &graph.chrom,
        graph.reference.len(),
        Some(&context),
    )?;
    stats.extend(result.stats);
    let structure = crate::structure::Plan::new(output, data)?;
    structure.validate(output, data, &graph.chrom)?;
    let emitted = result.emitted;
    require(
        positive.len() == graph.samples.len(),
        "positive model dimensions",
    )?;
    for (si, model) in prepared.models.iter().enumerate() {
        let expected: BTreeSet<_> = positive[si]
            .iter()
            .map(|(s, e, a, _)| (*s, *e, a.clone()))
            .collect();
        require(
            emitted[si] == expected,
            "final VCF differs from positive projection model",
        )?;
        let model_keys: BTreeSet<_> = model
            .iter()
            .map(|e| (e.start, e.end, e.alt.clone()))
            .collect();
        count(
            &mut stats,
            "unemitted_model_edits",
            model_keys.difference(&emitted[si]).count(),
        );
        let mut previous: Option<(usize, usize)> = None;
        for (s, e, _) in &emitted[si] {
            require(
                previous.is_none_or(|(ps, pe)| *s >= pe && !(*s == *e && *s == ps && ps == pe)),
                "overlapping emitted ALT genotypes",
            )?;
            previous = Some((*s, *e));
        }
    }
    let mut whole = Stats::new();
    for (si, paths) in graph.sample_paths.iter().enumerate() {
        if paths.len() != 1 {
            count(&mut whole, "multiple_paths_not_assessed", 1);
            continue;
        }
        let sequence = graph.spell(paths[0]);
        if !dna(&sequence) {
            count(&mut whole, "ambiguous_sequence_not_assessed", 1);
            continue;
        }
        let es: Vec<_> = emitted[si]
            .iter()
            .map(|(s, e, a)| Edit {
                start: *s,
                end: *e,
                alt: a.clone(),
                complex: false,
            })
            .collect();
        let actual = replay(
            &graph.reference,
            0,
            graph.reference.len(),
            &es.iter().collect::<Vec<_>>(),
        )?;
        count(
            &mut whole,
            if actual == sequence || actual == io::reverse_complement(&sequence) {
                "exact"
            } else {
                "incomplete"
            },
            1,
        );
    }
    let evidence = Evidence::new(
        graph.reference.clone(),
        prepared
            .models
            .iter()
            .map(|m| m.iter().map(Edit::input).collect())
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
    for (li, seqs) in data.alleles.iter().enumerate() {
        let (s, e) = data.bounds[li];
        require(
            graph.reference[s..e] == seqs[0],
            "joint REF differs from reference",
        )?;
        let mut seen = vec![false; seqs.len()];
        count(&mut stats, "joint_alleles", seqs.len() - 1);
        count(&mut stats, "joint_loci", 1);
        for (si, emitted_row) in emitted.iter().enumerate() {
            let g = data.gt(li, si);
            if g < 0 {
                continue;
            }
            seen[g as usize] = true;
            count(&mut stats, "joint_resolved_genotypes", 1);
            let row = &evidence.samples[si];
            let (lo, hi) = (s.saturating_sub(1), (e + 1).min(graph.reference.len()));
            let (inside_lo, inside_hi) = if s < e { (s, e) } else { (lo, hi) };
            require(
                (!row.structural.overlaps(inside_lo, inside_hi) || row.chain_covers(lo, hi))
                    && !row.uncertain.overlaps(inside_lo, inside_hi),
                "resolved complete locus intersects uncertainty",
            )?;
            if data.kinds[li] == "SNP" {
                require(
                    row.unique.covers(s, s + 1)
                        && row.base(&graph.reference, &evidence.pool, s)
                            == seqs[g as usize].as_bytes()[0],
                    "joint SNP/base coverage",
                )?;
                continue;
            }
            require(
                row.chain_covers(lo, hi) && row.unique.covers(inside_lo, inside_hi),
                "joint allele lacks unique covering source chain",
            )?;
            let ids = row.overlapping(&evidence.pool, lo, hi);
            let local: Vec<_> = ids
                .iter()
                .map(|&id| &evidence.pool[id])
                .filter(|x| intersects(s, e, x) && (s != e || x.structural()))
                .collect();
            require(
                local.iter().all(|x| s <= x.start && x.end <= e),
                "joint allele has spanning edits",
            )?;
            let observed = replay(&graph.reference, s, e, &local)?;
            require(observed == seqs[g as usize], "joint sequence replay")?;
            if data.projected(li, si) >= 0 {
                let actual: Vec<_> = local
                    .into_iter()
                    .filter(|x| emitted_row.contains(&(x.start, x.end, x.alt.clone())))
                    .collect();
                require(
                    replay(&graph.reference, s, e, &actual)? == observed,
                    "VCF component replay",
                )?;
                count(&mut stats, "joint_complete_VCF_replays", 1);
                if g > 0 && data.projected(li, si) == 0 {
                    count(&mut stats, "reference_backbone_with_emitted_components", 1);
                }
            }
        }
        require(
            seen.iter().skip(1).all(|&v| v),
            "unobserved complete ALT allele",
        )?;
    }
    Ok(
        json!({"status":"pass","samples":graph.samples.len(),"counts":stats,"whole_path_sequence_replay":whole,"structural_groups":structure.summary(),
        "scope":"Every native path base partitioned; every source allele matches its physical graph occurrence and replays after normalization; every emitted ALT has a source-derived model witness; confident genotypes require resolved coverage; REF, sample order, IDs, sort, AC/AN and projection serialization checked.",
        "limitation":"Unresolved source models remain missing or compound. A pass does not establish sensitivity, specificity, whole-genome haplotype completeness, or independent biological accuracy."}),
    )
}
