//! Deterministic structural allele groups, with exact catalogue provenance.
use crate::{
    catalogue::Catalogue,
    classify, io,
    model::{self, Edit, Result},
    variant,
};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs::File,
    io::Write,
    path::Path,
    process::Command,
};

pub const VCF: &str = "Graph.sv.vcf.gz";
pub const MEMBERS: &str = "Graph.sv_members.tsv.gz";
const COLUMNS: [&str; 5] = [
    "structure_id",
    "locus_id",
    "complete_allele_index",
    "event_allele_index",
    "role",
];
const SIMILARITY: usize = 95;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Allele {
    pub edit: Edit,
    pub reference: String,
    pub kind: String,
    pub reverse: Option<(usize, usize)>,
}
impl Allele {
    fn from_sequence(start: usize, reference: &str, alt: &str, a: &classify::Annotation) -> Self {
        let edit = model::trim(start, reference, alt);
        let offset = edit.start - start;
        let mut axes = vec![];
        for s in &a.segments {
            if s.query < offset || s.query_end < offset || s.query_end - offset > edit.alt.len() {
                continue;
            }
            let right = edit.alt.len() - (s.query_end - offset);
            if let Some(lo) = s.start.checked_sub(right) {
                axes.push((lo, s.end + s.query - offset));
            }
        }
        let tolerance = tolerance(edit.end - edit.start, edit.alt.len());
        let reverse = axes.first().copied().filter(|&(s, e)| {
            axes.iter()
                .all(|&(x, y)| s.abs_diff(x) <= tolerance && e.abs_diff(y) <= tolerance)
        });
        let reference = reference[edit.start - start..edit.end - start].into();
        Self {
            edit,
            reference,
            kind: a.kind.clone(),
            reverse,
        }
    }
    fn base_kind(&self) -> &str {
        if self.kind.starts_with("INV") {
            "INV"
        } else {
            &self.kind
        }
    }
    fn span(&self) -> usize {
        self.reference.len().max(self.edit.alt.len())
    }
}

fn tolerance(a: usize, b: usize) -> usize {
    (a.max(b) * (100 - SIMILARITY) / 100).clamp(1, 50)
}
fn similar_size(a: usize, b: usize) -> bool {
    (a == 0 && b == 0)
        || (a > 0 && b > 0 && a.min(b) as u128 * 100 >= a.max(b) as u128 * SIMILARITY as u128)
}

/// Exact bounded Levenshtein decision using furthest-reaching wavefronts.
/// Trim equal ends first; retained work depends on differences, not the full window area.
pub fn within_distance(mut a: &[u8], mut b: &[u8], budget: usize) -> bool {
    if a.len().abs_diff(b.len()) > budget {
        return false;
    }
    let prefix = a.iter().zip(b).take_while(|(x, y)| x == y).count();
    a = &a[prefix..];
    b = &b[prefix..];
    let suffix = a
        .iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| x == y)
        .count();
    a = &a[..a.len() - suffix];
    b = &b[..b.len() - suffix];
    if a.len().max(b.len()) <= budget {
        return true;
    }
    if a.is_empty() || b.is_empty() || budget == 0 {
        return false;
    }
    let offset = budget + 1;
    let mut previous = vec![-1isize; 2 * budget + 3];
    let mut current = previous.clone();
    previous[offset] = 0;
    for d in 1..=budget {
        for k in -(d as isize)..=d as isize {
            let i = (offset as isize + k) as usize;
            let mut best = -1;
            for x in [previous[i] + 1, previous[i - 1] + 1, previous[i + 1]] {
                let y = x - k;
                if x >= 0 && y >= 0 && x <= a.len() as isize && y <= b.len() as isize {
                    best = best.max(x);
                }
            }
            if best < 0 {
                current[i] = -1;
                continue;
            }
            let mut x = best as usize;
            let mut y = (best - k) as usize;
            while x < a.len() && y < b.len() && a[x] == b[y] {
                x += 1;
                y += 1;
            }
            if x == a.len() && y == b.len() {
                return true;
            }
            current[i] = x as isize;
        }
        std::mem::swap(&mut previous, &mut current);
    }
    false
}

pub fn compatible(a: &Allele, b: &Allele) -> bool {
    if a.base_kind() != b.base_kind() {
        return false;
    }
    if a.edit == b.edit {
        return true;
    }
    let delta = tolerance(a.span().min(b.span()), 0);
    if a.edit.start.abs_diff(b.edit.start) > delta
        || a.edit.end.abs_diff(b.edit.end) > delta
        || !similar_size(a.reference.len(), b.reference.len())
        || !similar_size(a.edit.alt.len(), b.edit.alt.len())
    {
        return false;
    }
    if a.base_kind() == "INV" {
        let (Some((s, e)), Some((x, y))) = (a.reverse, b.reverse) else {
            return false;
        };
        if s.abs_diff(x) > delta || e.abs_diff(y) > delta {
            return false;
        }
    } else if a.kind == "COMPLEX" {
        return false;
    }
    if a.kind == "DEL" {
        return within_distance(
            a.reference.as_bytes(),
            b.reference.as_bytes(),
            a.reference.len().max(b.reference.len()) * (100 - SIMILARITY) / 100,
        );
    }
    let start = a.edit.start.min(b.edit.start);
    let end = a.edit.end.max(b.edit.end);
    // A common reference frame must be known, not interpolated across a gap.
    let mut reference = vec![0; end - start];
    for x in [a, b] {
        for (i, base) in x.reference.bytes().enumerate() {
            let p = x.edit.start - start + i;
            if reference[p] != 0 && reference[p] != base {
                return false;
            }
            reference[p] = base;
        }
    }
    if reference.contains(&0) {
        return false;
    }
    let haplotype = |x: &Allele| {
        let mut v = reference[..x.edit.start - start].to_vec();
        v.extend_from_slice(x.edit.alt.as_bytes());
        v.extend_from_slice(&reference[x.edit.end - start..]);
        v
    };
    let (x, y) = (haplotype(a), haplotype(b));
    let budget = x.len().max(y.len()) * (100 - SIMILARITY) / 100;
    within_distance(&x, &y, budget)
        && !within_distance(&x, &reference, budget)
        && !within_distance(&y, &reference, budget)
}

#[derive(Clone, Debug)]
struct Member {
    locus: usize,
    allele: usize,
    event: Option<usize>,
}
#[derive(Clone, Debug)]
struct Seed {
    allele: Allele,
    members: Vec<Member>,
    anchor: u8,
}
#[derive(Clone, Debug)]
pub struct Group {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub kind: String,
    pub gt: Vec<i32>,
    pub status: Vec<&'static str>,
    seeds: Vec<usize>,
    members: Vec<Member>,
    anchor: u8,
}

fn seeds(output: &Path, data: &Catalogue, annotations: &classify::Table) -> Result<Vec<Seed>> {
    let mut reader = variant::Reader::open(&output.join("Graph.vcf.gz"))?;
    let mut items: BTreeMap<Allele, Seed> = BTreeMap::new();
    while let Some(r) = reader.next_record()? {
        let li = *data.index.get(&r.id).ok_or("Unknown event locus")?;
        let sizes = r
            .info
            .get("SIZECLASS")
            .ok_or("Missing SIZECLASS")?
            .split(',')
            .collect::<Vec<_>>();
        if sizes.len() != r.alts.len() {
            return Err("Invalid SV classification cardinality".into());
        }
        for (i, alt) in r.alts.iter().enumerate().filter(|(i, _)| sizes[*i] == "SV") {
            let edit = model::trim(r.pos - 1, &r.reference, alt);
            let ai = data.alleles[li]
                .iter()
                .position(|s| model::trim(data.bounds[li].0, &data.alleles[li][0], s) == edit)
                .ok_or("SV ALT absent from complete catalogue")?;
            if ai == 0 {
                return Err("SV ALT equals REF".into());
            }
            let allele = Allele::from_sequence(
                data.bounds[li].0,
                &data.alleles[li][0],
                &data.alleles[li][ai],
                annotations.get(li, ai),
            );
            let position = edit.start.saturating_sub(1);
            let anchor = *r
                .reference
                .as_bytes()
                .get(
                    position
                        .checked_sub(r.pos - 1)
                        .ok_or("SV anchor precedes record")?,
                )
                .ok_or("SV reference anchor unavailable")?;
            let m = Member {
                locus: li,
                allele: ai,
                event: Some(i + 1),
            };
            items
                .entry(allele.clone())
                .and_modify(|s| s.members.push(m.clone()))
                .or_insert(Seed {
                    allele,
                    members: vec![m],
                    anchor,
                });
        }
    }
    Ok(items.into_values().collect())
}

/// Complete linkage, deterministic sequence ordering and no carrier-frequency input.
pub fn cluster(alleles: &[Allele]) -> Vec<Vec<usize>> {
    let mut order: Vec<_> = (0..alleles.len()).collect();
    order.sort_by_key(|&i| &alleles[i]);
    let mut result: Vec<Vec<usize>> = vec![];
    for i in order {
        if let Some(group) = result.iter_mut().rev().find(|g| {
            let r = &alleles[g[0]];
            let a = &alleles[i];
            r.edit.start.abs_diff(a.edit.start) <= 50
                && g.iter().all(|&j| compatible(a, &alleles[j]))
        }) {
            group.push(i);
        } else {
            result.push(vec![i]);
        }
    }
    result
}

fn observed(
    data: &Catalogue,
    annotations: &classify::Table,
    li: usize,
    ai: usize,
) -> Option<Allele> {
    (ai > 0).then(|| {
        Allele::from_sequence(
            data.bounds[li].0,
            &data.alleles[li][0],
            &data.alleles[li][ai],
            annotations.get(li, ai),
        )
    })
}

pub fn resolve(positive: bool, states: &[i32], covered: bool) -> (i32, &'static str) {
    let positive = positive || states.contains(&1);
    if positive && covered && states.contains(&0) {
        return (-1, "AMBIGUOUS_ALIGNMENT");
    }
    if positive {
        return (1, "PASS");
    }
    if covered && !states.is_empty() && states.iter().all(|&s| s == 0) {
        return (0, "PASS");
    }
    (-1, "NO_DATA")
}

pub struct Plan {
    pub groups: Vec<Group>,
    seed_count: usize,
}
impl Plan {
    pub fn new(output: &Path, data: &Catalogue) -> Result<Self> {
        let annotations = classify::Table::read(output, data, None)?;
        let seeds = seeds(output, data, &annotations)?;
        let alleles: Vec<_> = seeds.iter().map(|s| s.allele.clone()).collect();
        let clusters = cluster(&alleles);
        let mut nearby: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, g) in clusters.iter().enumerate() {
            nearby.entry(alleles[g[0]].edit.start).or_default().push(i);
        }
        let mut assigned: HashMap<Edit, usize> = HashMap::new();
        for (gi, g) in clusters.iter().enumerate() {
            for &i in g {
                if let Some(old) = assigned.insert(alleles[i].edit.clone(), gi) {
                    if old != gi {
                        return Err("Conflicting structural class for identical edit".into());
                    }
                }
            }
        }
        let mut groups = vec![];
        for (gi, g) in clusters.iter().enumerate() {
            let start = g.iter().map(|&i| alleles[i].edit.start).min().unwrap();
            let end = g.iter().map(|&i| alleles[i].edit.end).max().unwrap();
            let mut kind = if g.iter().any(|&i| alleles[i].kind == "INV_COMPLEX") {
                "INV_COMPLEX"
            } else {
                &alleles[g[0]].kind
            }
            .to_string();
            let mut hash = Sha256::new();
            for &i in g {
                let a = &alleles[i];
                hash.update(format!(
                    "{}\t{}\t{}\t{}\n",
                    a.edit.start, a.edit.end, a.kind, a.edit.alt
                ));
            }
            let id = format!("SV{start}_{:.16}", format!("{:x}", hash.finalize()));
            let mut members: Vec<_> = g.iter().flat_map(|&i| seeds[i].members.clone()).collect();
            let source_loci: BTreeSet<_> = members.iter().map(|m| m.locus).collect();
            let cover: Vec<_> = source_loci
                .iter()
                .copied()
                .filter(|&li| data.bounds[li].0 <= start && data.bounds[li].1 >= end)
                .collect();
            let shortest = cover
                .iter()
                .map(|&li| data.bounds[li].1 - data.bounds[li].0)
                .min();
            let assays: Vec<_> = if let Some(n) = shortest {
                cover
                    .into_iter()
                    .filter(|&li| data.bounds[li].1 - data.bounds[li].0 == n)
                    .collect()
            } else {
                source_loci.iter().copied().collect()
            };
            // Each complete allele is assigned once, independent of sample frequency.
            let mut states = HashMap::new();
            for &li in &assays {
                for ai in 0..data.alleles[li].len() {
                    let state = if let Some(a) = observed(data, &annotations, li, ai) {
                        if let Some(&assigned) = assigned.get(&a.edit) {
                            if assigned == gi {
                                1
                            } else if a.edit.start < start || a.edit.end > end {
                                -1
                            } else {
                                0
                            }
                        } else {
                            let matches: Vec<_> = nearby
                                .range(
                                    a.edit.start.saturating_sub(50)
                                        ..=a.edit.start.saturating_add(50),
                                )
                                .flat_map(|(_, v)| v.iter().copied())
                                .filter(|&j| {
                                    clusters[j].iter().all(|&k| compatible(&a, &alleles[k]))
                                })
                                .collect();
                            if matches.contains(&gi) {
                                if matches.len() == 1 {
                                    1
                                } else {
                                    -1
                                }
                            } else if a.edit.start < start || a.edit.end > end {
                                -1
                            } else {
                                0
                            }
                        }
                    } else {
                        0
                    };
                    states.insert((li, ai), state);
                    if state == 1 && !members.iter().any(|m| m.locus == li && m.allele == ai) {
                        members.push(Member {
                            locus: li,
                            allele: ai,
                            event: None,
                        });
                    }
                }
            }
            // Extra complete-allele matches must also be compatible with each other.
            // Ambiguous extensions remain unknown; they cannot create a transitive bridge.
            let additions: Vec<_> = members
                .iter()
                .filter(|m| m.event.is_none())
                .map(|m| {
                    (
                        (m.locus, m.allele),
                        observed(data, &annotations, m.locus, m.allele).unwrap(),
                    )
                })
                .collect();
            let rejected: BTreeSet<_> = additions
                .iter()
                .filter(|(_, a)| additions.iter().any(|(_, b)| !compatible(a, b)))
                .map(|(key, _)| *key)
                .collect();
            for key in &rejected {
                states.insert(*key, -1);
            }
            members.retain(|m| m.event.is_some() || !rejected.contains(&(m.locus, m.allele)));
            if kind == "INV"
                && members
                    .iter()
                    .any(|m| annotations.get(m.locus, m.allele).kind == "INV_COMPLEX")
            {
                kind = "INV_COMPLEX".into();
            }
            let mut gt = vec![];
            let mut status = vec![];
            for si in 0..data.samples.len() {
                let positive = members
                    .iter()
                    .any(|m| data.gt(m.locus, si) == m.allele as i32);
                let states: Vec<_> = assays
                    .iter()
                    .map(|&li| {
                        let ai = data.gt(li, si);
                        if ai < 0 {
                            -1
                        } else {
                            states[&(li, ai as usize)]
                        }
                    })
                    .collect();
                let (value, mut reason) = resolve(positive, &states, shortest.is_some());
                if value < 0 && reason == "NO_DATA" {
                    let reasons: BTreeSet<_> = assays
                        .iter()
                        .filter(|&&li| data.gt(li, si) < 0)
                        .map(|&li| data.status(li, si, false))
                        .collect();
                    reason = if reasons.len() == 1 && assays.iter().all(|&li| data.gt(li, si) < 0) {
                        crate::project::STATUS
                            .iter()
                            .copied()
                            .find(|s| reasons.contains(s))
                            .unwrap_or("NO_DATA")
                    } else {
                        "UNRESOLVED_COMPONENT"
                    };
                }
                gt.push(value);
                status.push(reason);
            }
            members.sort_by_key(|m| (&data.loci[m.locus], m.allele, m.event));
            members.dedup_by_key(|m| (m.locus, m.allele, m.event));
            let anchor = seeds[g
                .iter()
                .copied()
                .min_by_key(|&i| alleles[i].edit.start)
                .unwrap()]
            .anchor;
            groups.push(Group {
                id,
                start,
                end,
                kind,
                gt,
                status,
                seeds: g.clone(),
                members,
                anchor,
            });
        }
        groups.sort_by(|a, b| {
            (a.start.saturating_sub(1), &a.id).cmp(&(b.start.saturating_sub(1), &b.id))
        });
        Ok(Self {
            groups,
            seed_count: seeds.len(),
        })
    }
    fn member_rows(&self, data: &Catalogue) -> Vec<Vec<String>> {
        self.groups
            .iter()
            .flat_map(|g| {
                g.members.iter().map(|m| {
                    vec![
                        g.id.clone(),
                        data.loci[m.locus].clone(),
                        m.allele.to_string(),
                        m.event.map_or_else(|| ".".into(), |a| a.to_string()),
                        if m.event.is_some() {
                            "SEED"
                        } else {
                            "COMPLETE"
                        }
                        .into(),
                    ]
                })
            })
            .collect()
    }
    fn line(g: &Group, chrom: &str) -> String {
        let symbolic = match g.kind.as_str() {
            "INS" => "INS",
            "DEL" => "DEL",
            "INV" => "INV",
            _ => "CPX",
        };
        let ac = g.gt.iter().filter(|&&v| v == 1).count();
        let an = g.gt.iter().filter(|&&v| v >= 0).count();
        let cells =
            g.gt.iter()
                .zip(&g.status)
                .map(|(&v, s)| format!("{}:{s}", if v < 0 { ".".into() } else { v.to_string() }))
                .collect::<Vec<_>>()
                .join("\t");
        format!("{chrom}\t{}\t{}\t{}\t<{symbolic}>\t.\tPASS\tEND={};START0={};END0={};SVTYPE={symbolic};TYPE={};NALLELES={};AC={ac};AN={an}\tGT:FT\t{cells}",g.start.max(1),g.id,g.anchor as char,g.end.max(g.start.max(1)),g.start,g.end,g.kind,g.members.len())
    }
    pub fn write(
        &self,
        output: &Path,
        data: &Catalogue,
        chrom: &str,
        length: usize,
        threads: usize,
    ) -> Result<()> {
        let mut header=vec!["##fileformat=VCFv4.2".into(),format!("##source=PathAllele-{}",env!("CARGO_PKG_VERSION")),format!("##contig=<ID={chrom},length={length}>"),"##pathallele_structural_groups=\"Phenotype-independent complete-link groups; 95 percent sequence and size similarity; both endpoints within 5 percent (maximum 50 bp); compatible inversion axes. Not independent ancestral mutation counts.\"".into(),"##pathallele_structural_gt=\"1=observed member structure; 0=resolved absence in a covering complete locus; .=unresolved/conflicting. Symbolic ALT does not assert one exact carrier sequence.\"".into()];
        header.extend(variant::headers("structural_allele_group"));
        for (id, description) in [
            ("INS", "Insertion structural allele group"),
            ("DEL", "Deletion structural allele group"),
            ("INV", "Exact inversion structural allele group"),
            (
                "CPX",
                "Compound or sequence-divergent structural allele group; see TYPE",
            ),
        ] {
            header.push(format!("##ALT=<ID={id},Description=\"{description}\">"));
        }
        for (id, typ, desc) in [
            (
                "END",
                "Integer",
                "End of structural comparison envelope; not inferred historical breakpoint",
            ),
            (
                "START0",
                "Integer",
                "Unpadded comparison-envelope start, zero-based",
            ),
            (
                "END0",
                "Integer",
                "Unpadded comparison-envelope end, zero-based exclusive",
            ),
            ("SVTYPE", "String", "Symbolic structural category"),
            (
                "TYPE",
                "String",
                "Preserved evidence class including INV_COMPLEX and REPLACEMENT",
            ),
            (
                "NALLELES",
                "Integer",
                "Complete catalogue memberships, including repeated views",
            ),
            ("AC", "Integer", "Haploid structure carrier count"),
            ("AN", "Integer", "Haploid callable sample count"),
        ] {
            header.push(format!(
                "##INFO=<ID={id},Number=1,Type={typ},Description=\"{desc}\">"
            ));
        }
        let mut writer = noodles_bgzf::MultithreadedWriter::with_worker_count(
            std::num::NonZeroUsize::new(threads).ok_or("Invalid compression threads")?,
            File::create(output.join(VCF)).map_err(|e| e.to_string())?,
        );
        writeln!(
            writer,
            "{}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            header.join("\n"),
            data.samples.join("\t")
        )
        .map_err(|e| e.to_string())?;
        for g in &self.groups {
            writeln!(writer, "{}", Self::line(g, chrom)).map_err(|e| e.to_string())?;
        }
        writer
            .finish()
            .map_err(|e| e.to_string())?
            .sync_all()
            .map_err(|e| e.to_string())?;
        let mut table = io::tsv(&output.join(MEMBERS))?;
        table.write_record(COLUMNS).map_err(|e| e.to_string())?;
        for row in self.member_rows(data) {
            table.write_record(row).map_err(|e| e.to_string())?;
        }
        io::finish_tsv(table)?;
        if !Command::new("bcftools")
            .args(["index", "-t"])
            .arg(output.join(VCF))
            .status()
            .map_err(|e| e.to_string())?
            .success()
        {
            return Err("Structural VCF indexing failed".into());
        }
        Ok(())
    }
    pub fn validate(&self, output: &Path, data: &Catalogue, chrom: &str) -> Result<()> {
        use std::io::BufRead;
        let mut reader = variant::Reader::open(&output.join(VCF))?;
        if !reader.header.sample_names().iter().eq(data.samples.iter()) {
            return Err("Structural sample order differs".into());
        }
        for g in &self.groups {
            let row = reader
                .next_record()?
                .ok_or("Missing structural VCF record")?;
            if row.id != g.id || row.chrom != chrom || row.gt != g.gt || row.ft != g.status {
                return Err(
                    "Structural genotype or record identity disagrees with catalogue".into(),
                );
            }
        }
        if reader.next_record()?.is_some() {
            return Err("Unexpected structural record".into());
        }
        let lines = io::reader(&output.join(VCF))?
            .lines()
            .filter_map(|l| match l {
                Ok(s) if !s.starts_with('#') => Some(Ok(s)),
                Err(e) => Some(Err(e.to_string())),
                _ => None,
            })
            .collect::<Result<Vec<_>>>()?;
        if lines
            != self
                .groups
                .iter()
                .map(|g| Self::line(g, chrom))
                .collect::<Vec<_>>()
        {
            return Err(
                "Structural coordinates, alleles, counts or status contract differs".into(),
            );
        }
        let mut table = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(io::reader(&output.join(MEMBERS))?);
        if !table
            .headers()
            .map_err(|e| e.to_string())?
            .iter()
            .eq(COLUMNS)
        {
            return Err("Invalid structural membership header".into());
        }
        let actual = table
            .records()
            .map(|r| {
                r.map(|r| r.iter().map(str::to_string).collect())
                    .map_err(|e| e.to_string())
            })
            .collect::<Result<Vec<Vec<String>>>>()?;
        if actual != self.member_rows(data) {
            return Err("Structural allele membership differs".into());
        }
        Ok(())
    }
    pub fn summary(&self) -> serde_json::Value {
        serde_json::json!({"schema":1,"groups":self.groups.len(),"unique_exact_seed_alleles":self.seed_count,"merged_groups":self.groups.iter().filter(|g|g.seeds.len()>1).count(),"similarity_percent":SIMILARITY,"scope":"Structural allele analysis units, not independent ancestral events"})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn common_padding_preserves_inversion_identity_and_axes() {
        let mut seed = 77319u64;
        let reference: String = (0..300)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed & 3) as usize] as char
            })
            .collect();
        let alt = io::reverse_complement(&reference[100..184]);
        let padded = format!("{}{alt}{}", &reference[90..100], &reference[184..194]);
        let c = classify::Classifier::new(&reference, 31).unwrap();
        let a = Allele::from_sequence(
            100,
            &reference[100..184],
            &alt,
            &c.annotate(100, 184, &alt, "SEQUENCE").unwrap(),
        );
        let b = Allele::from_sequence(
            90,
            &reference[90..194],
            &padded,
            &c.annotate(90, 194, &padded, "SEQUENCE").unwrap(),
        );
        assert!(a.reverse.is_some());
        assert_eq!(a, b);
    }
}
