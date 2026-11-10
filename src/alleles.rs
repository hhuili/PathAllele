//! Complete allele units, canonical normalization and conservative partial projection.
use crate::{
    alignment::Adapter,
    anchors::{self, Anchor},
    io,
    model::{self, atomic_edits, dna, intersects, replay, Coverage, Edit, Ranges, Result},
    source::{count, Parent, Sources, Stats},
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, Write};
use std::path::Path;
use std::process::{Command, Stdio};

pub type Decomposition = Vec<Vec<Vec<Edit>>>;
pub type Models = Vec<BTreeSet<Edit>>;
type Keys = BTreeSet<(usize, usize)>;
pub type Locus = (usize, usize, String);
pub type OriginMap = BTreeMap<Locus, BTreeSet<String>>;
pub type KindMap = BTreeMap<(usize, usize, String), BTreeSet<String>>;

fn empty(parents: &[Parent]) -> Decomposition {
    parents
        .iter()
        .map(|p| vec![vec![]; p.alleles.len()])
        .collect()
}
fn local_replay(reference: &str, edits: &[Edit]) -> Result<String> {
    replay(
        reference,
        0,
        reference.len(),
        &edits.iter().collect::<Vec<_>>(),
    )
}
fn absolute(e: &Edit, offset: usize) -> Edit {
    Edit {
        start: e.start + offset,
        end: e.end + offset,
        alt: e.alt.clone(),
        complex: e.complex,
    }
}

pub fn reverse_parents(parents: &[Parent]) -> Vec<Parent> {
    parents
        .iter()
        .map(|p| {
            let mut r = p.clone();
            r.alleles = p
                .alleles
                .iter()
                .map(|s| io::reverse_complement(s))
                .collect();
            r.anchors = p
                .anchors
                .iter()
                .zip(&p.alleles)
                .map(|(aa, s)| {
                    aa.iter()
                        .rev()
                        .map(|a| Anchor {
                            r: p.alleles[0].len() - a.r - a.n,
                            q: s.len() - a.q - a.n,
                            n: a.n,
                            node: a.node,
                        })
                        .collect()
                })
                .collect();
            r
        })
        .collect()
}
fn reflect(parents: &[Parent], raw: Decomposition) -> Decomposition {
    parents
        .iter()
        .zip(raw)
        .map(|(p, aa)| {
            aa.into_iter()
                .map(|es| {
                    es.into_iter()
                        .map(|e| Edit {
                            start: p.alleles[0].len() - e.end,
                            end: p.alleles[0].len() - e.start,
                            alt: io::reverse_complement(&e.alt),
                            complex: e.complex,
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}
pub fn decompose(parents: &[Parent]) -> Result<Decomposition> {
    let mut result = empty(parents);
    for (p, ee) in parents.iter().zip(&mut result) {
        for (ai, a) in p.alleles.iter().enumerate().skip(1) {
            if matches!(p.kinds[ai].as_str(), "INV" | "COMPLEX") {
                ee[ai] = vec![Edit {
                    start: 1,
                    end: p.alleles[0].len() - 1,
                    alt: a[1..a.len() - 1].into(),
                    complex: true,
                }];
            } else {
                for (s, r, a) in anchors::gaps(&p.alleles[0], a, &p.anchors[ai])? {
                    ee[ai].extend(atomic_edits(s, r, a, !r.is_empty() && !a.is_empty()));
                }
            }
            if local_replay(&p.alleles[0], &ee[ai])? != *a {
                return Err(format!(
                    "Sequence units failed to replay: {} allele {ai}",
                    p.id
                ));
            }
        }
    }
    Ok(result)
}
fn bounded_indels(reference: &str, edits: &[Edit]) -> Vec<Edit> {
    edits
        .iter()
        .map(|e| {
            let mut x = e.clone();
            if x.start == x.end && !x.alt.is_empty() {
                while x.start > 1
                    && reference.as_bytes()[x.start - 1] == *x.alt.as_bytes().last().unwrap()
                {
                    let last = x.alt.pop().unwrap();
                    x.alt.insert(0, last);
                    x.start -= 1;
                    x.end -= 1;
                }
            } else if x.alt.is_empty() && x.start < x.end {
                while x.start > 1
                    && reference.as_bytes()[x.start - 1] == reference.as_bytes()[x.end - 1]
                {
                    x.start -= 1;
                    x.end -= 1;
                }
            }
            x
        })
        .collect()
}

pub fn normalize(
    parents: &[Parent],
    edits: &Decomposition,
    reference: &str,
    chrom: &str,
    reference_path: &Path,
    work: &Path,
) -> Result<(Decomposition, Keys)> {
    fs::create_dir(work).map_err(|e| e.to_string())?;
    let input = work.join("input.vcf");
    let output = work.join("normalized.vcf");
    let mut f = io::writer(&input)?;
    writeln!(f,"##fileformat=VCFv4.2\n##contig=<ID={chrom},length={}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",reference.len()).map_err(|e|e.to_string())?;
    for (pi, (p, aa)) in parents.iter().zip(edits).enumerate() {
        for (ai, es) in aa.iter().enumerate() {
            for (ei, e) in es.iter().enumerate() {
                let (pos, r, a) = model::vcf_alleles(
                    reference,
                    p.start + e.start,
                    p.start + e.end,
                    &[e.alt.clone()],
                    !e.structural(),
                )?;
                writeln!(f, "{chrom}\t{pos}\te{pi}_{ai}_{ei}\t{r}\t{}\t.\t.\t.", a[0])
                    .map_err(|e| e.to_string())?;
            }
        }
    }
    f.finish()?;
    let status = Command::new("bcftools")
        .args(["norm", "-c", "e", "-f"])
        .arg(reference_path)
        .arg(&input)
        .args(["-Ov", "-o"])
        .arg(&output)
        .stderr(Stdio::from(
            File::create(work.join("normalization.log")).map_err(|e| e.to_string())?,
        ))
        .status()
        .map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!("bcftools normalization failed: {}", work.display()));
    }
    // Parse canonical edits in absolute coordinates: left normalization may
    // cross a parent boundary. Subtract the parent offset only after replay.
    let mut result = empty(parents);
    let mut seen = HashSet::new();
    for line in io::reader(&output)?.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('#') {
            continue;
        }
        let x: Vec<_> = line.split('\t').collect();
        if x.len() < 8 {
            return Err("Truncated normalized VCF".into());
        }
        let id = x[2]
            .strip_prefix('e')
            .ok_or("Unexpected normalized record ID")?;
        let ids: Vec<usize> = id
            .split('_')
            .map(|s| {
                s.parse()
                    .map_err(|_| "Invalid normalized record ID".to_string())
            })
            .collect::<Result<_>>()?;
        if ids.len() != 3 {
            return Err("Invalid normalized record ID".into());
        }
        let (pi, ai, ei) = (ids[0], ids[1], ids[2]);
        let original = edits
            .get(pi)
            .and_then(|aa| aa.get(ai))
            .and_then(|ee| ee.get(ei))
            .ok_or("Unknown normalized record ID")?;
        if !seen.insert((pi, ai, ei)) {
            return Err("Duplicate normalized record ID".into());
        }
        let start = x[1]
            .parse::<usize>()
            .map_err(|_| "Invalid VCF POS")?
            .checked_sub(1)
            .ok_or("VCF POS must be positive")?;
        result[pi][ai].extend(atomic_edits(start, x[3], x[4], original.complex));
    }
    if seen.len() != edits.iter().flatten().map(Vec::len).sum::<usize>() {
        return Err("Normalizer dropped input edits".into());
    }
    let mut failed = Keys::new();
    for (pi, (p, aa)) in parents.iter().zip(&mut result).enumerate() {
        for (ai, es) in aa.iter_mut().enumerate() {
            let lo = es
                .iter()
                .map(|e| e.start)
                .min()
                .unwrap_or(p.start)
                .min(p.start);
            let hi = es
                .iter()
                .map(|e| e.end)
                .max()
                .unwrap_or(0)
                .max(p.start + p.alleles[0].len());
            let expected = format!(
                "{}{}{}",
                &reference[lo..p.start],
                p.alleles[ai],
                &reference[p.start + p.alleles[0].len()..hi]
            );
            let mut valid = replay(reference, lo, hi, &es.iter().collect::<Vec<_>>())
                .is_ok_and(|s| s == expected);
            let bounded = es
                .iter()
                .all(|e| p.start < e.start && e.end < p.start + p.alleles[0].len());
            if valid && !bounded {
                let candidate = bounded_indels(&p.alleles[0], &edits[pi][ai]);
                valid = local_replay(&p.alleles[0], &candidate).is_ok_and(|s| s == p.alleles[ai]);
                if valid {
                    *es = candidate;
                    continue;
                }
            }
            if !valid {
                failed.insert((pi, ai));
                *es = edits[pi][ai].clone();
            } else {
                for e in es {
                    e.start -= p.start;
                    e.end -= p.start;
                }
            }
        }
    }
    Ok((result, failed))
}

pub fn disagreements(
    reference: &str,
    offset: usize,
    left: &[Edit],
    right: &[Edit],
) -> BTreeSet<usize> {
    let l: BTreeSet<_> = left.iter().collect();
    let r: BTreeSet<_> = right.iter().collect();
    if l == r {
        return BTreeSet::new();
    }
    let all: Vec<_> = l.union(&r).copied().collect();
    let lo = (offset + all.iter().map(|e| e.start).min().unwrap()).saturating_sub(1);
    let hi = (offset + all.iter().map(|e| e.end).max().unwrap() + 1).min(reference.len());
    let bases = |es: &[Edit]| {
        let mut out = reference.as_bytes()[lo..hi].to_vec();
        for e in es {
            let (s, t) = (offset + e.start - lo, offset + e.end - lo);
            if e.structural() {
                out[s..t].fill(b'?');
            } else {
                out[s..t].copy_from_slice(e.alt.as_bytes());
            }
        }
        out
    };
    let mut positions: BTreeSet<_> = bases(left)
        .iter()
        .zip(bases(right))
        .enumerate()
        .filter_map(|(i, (&a, b))| (a != b).then_some(lo + i))
        .collect();
    for e in l.symmetric_difference(&r) {
        if e.structural() {
            positions.extend(
                (offset + e.start).saturating_sub(1)..(offset + e.end + 1).min(reference.len()),
            );
        }
    }
    positions
}
fn preserve(
    parents: &[Parent],
    forward: &Decomposition,
    reverse: &Decomposition,
    nf: &Decomposition,
    nr: &Decomposition,
    failed: &Keys,
    reference: &str,
) -> Result<(Decomposition, Decomposition, usize)> {
    let (mut f, mut r) = (forward.clone(), reverse.clone());
    let mut changed = 0;
    for (pi, p) in parents.iter().enumerate() {
        for ai in 1..p.alleles.len() {
            if disagreements(reference, p.start, &nf[pi][ai], &nr[pi][ai]).is_empty()
                && !failed.contains(&(pi, ai))
            {
                continue;
            }
            let raw = disagreements(reference, p.start, &forward[pi][ai], &reverse[pi][ai]);
            let mut components = vec![];
            for (s, rf, al) in anchors::gaps(&p.alleles[0], &p.alleles[ai], &p.anchors[ai])? {
                if failed.contains(&(pi, ai))
                    || raw
                        .range((p.start + s).saturating_sub(1)..=p.start + s + rf.len())
                        .next()
                        .is_some()
                {
                    components.push((s, rf, al));
                }
            }
            if components.is_empty() {
                return Err(format!(
                    "Unstable allele without bounded gap: {}:{ai}",
                    p.id
                ));
            }
            for dst in [&mut f, &mut r] {
                let keep = &mut dst[pi][ai];
                for &(s, rf, al) in &components {
                    let e = s + rf.len();
                    keep.retain(|x| !(s <= x.start && x.end <= e));
                    let c = !(rf.is_empty() || al.is_empty() || rf.len() == 1 && al.len() == 1);
                    keep.extend(atomic_edits(s, rf, al, c));
                }
                if local_replay(&p.alleles[0], keep)? != p.alleles[ai] {
                    return Err("Preserved gap failed source replay".into());
                }
            }
            changed += components.len();
        }
    }
    Ok((f, r, changed))
}
fn close_touching(reference: &str, edits: &[Edit]) -> Result<Vec<Edit>> {
    let ordered: BTreeSet<_> = edits.iter().collect();
    let mut groups: Vec<Vec<&Edit>> = vec![];
    for e in ordered {
        if let Some(g) = groups.last_mut() {
            if e.start <= g.iter().map(|x| x.end).max().unwrap() {
                g.push(e);
                continue;
            }
        }
        groups.push(vec![e]);
    }
    let mut result = vec![];
    for g in groups {
        if g.len() == 1 {
            result.push(g[0].clone());
            continue;
        }
        let s = g.iter().map(|e| e.start).min().unwrap();
        let t = g.iter().map(|e| e.end).max().unwrap();
        let alt = replay(reference, s, t, &g)?;
        result.extend(atomic_edits(s, &reference[s..t], &alt, true));
    }
    if local_replay(reference, &result)? != local_replay(reference, edits)? {
        return Err("Closing normalized units changed sequence".into());
    }
    Ok(result)
}

fn align(
    parents: &[Parent],
    chrom: &str,
    work: &Path,
    jobs: usize,
    executable: &Adapter,
) -> Result<Decomposition> {
    fs::create_dir(work).map_err(|e| e.to_string())?;
    let mut gaps = BTreeMap::new();
    let mut tasks = vec![];
    for (pi, p) in parents.iter().enumerate() {
        for ai in 1..p.alleles.len() {
            if p.kinds[ai] != "SEQUENCE" {
                continue;
            }
            let parts = anchors::gaps(&p.alleles[0], &p.alleles[ai], &p.anchors[ai])?;
            for (segment, &(_, r, a)) in parts.iter().enumerate() {
                if !r.is_empty() && !a.is_empty() && dna(r) && dna(a) {
                    tasks.push((r.len() + a.len(), pi, ai, segment));
                }
            }
            gaps.insert((pi, ai), parts);
        }
    }
    let jobs = jobs.min(tasks.len()).max(1);
    let mut buckets = vec![vec![]; jobs];
    let mut loads = vec![0; jobs];
    tasks.sort_unstable_by(|a, b| b.cmp(a));
    for (weight, pi, ai, seg) in tasks {
        let j = (0..jobs).min_by_key(|&i| loads[i]).unwrap();
        buckets[j].push((pi, ai, seg));
        loads[j] += weight;
    }
    for (j, bucket) in buckets.iter().enumerate() {
        let mut f = io::writer(&work.join(format!("{j}.alleles.vcf")))?;
        writeln!(f,"##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Synthetic allele carrier; ignored in inference\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tallele").map_err(|e|e.to_string())?;
        for &(pi, ai, seg) in bucket {
            let (s, r, a) = gaps[&(pi, ai)][seg];
            writeln!(
                f,
                "{chrom}\t{}\ta{pi}_{ai}_{seg}\t{r}\t{a}\t.\t.\t.\tGT\t1",
                parents[pi].start + 1 + s
            )
            .map_err(|e| e.to_string())?;
        }
        f.finish()?;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .map_err(|e| e.to_string())?;
    let outcomes: Vec<Result<()>> = pool.install(|| {
        (0..jobs)
            .into_par_iter()
            .map(|j| {
                let status = executable
                    .command()
                    .arg(work.join(format!("{j}.alleles.vcf")))
                    .env("OMP_NUM_THREADS", "1")
                    .env("OPENBLAS_NUM_THREADS", "1")
                    .stdout(Stdio::from(
                        File::create(work.join(format!("{j}.edits.tsv")))
                            .map_err(|e| e.to_string())?,
                    ))
                    .stderr(Stdio::from(
                        File::create(work.join(format!("{j}.log"))).map_err(|e| e.to_string())?,
                    ))
                    .status()
                    .map_err(|e| e.to_string())?;
                if status.success() {
                    Ok(())
                } else {
                    Err(format!("BiWFA adapter failed: batch {j}"))
                }
            })
            .collect()
    });
    for x in outcomes {
        x?;
    }
    let mut result = empty(parents);
    for (&(pi, ai), parts) in &gaps {
        for &(s, r, a) in parts {
            if !dna(r) || !dna(a) {
                result[pi][ai].push(Edit {
                    start: s,
                    end: s + r.len(),
                    alt: a.into(),
                    complex: true,
                });
            } else if r.is_empty() || a.is_empty() {
                result[pi][ai].extend(atomic_edits(s, r, a, false));
            }
        }
    }
    for j in 0..jobs {
        for line in io::reader(&work.join(format!("{j}.edits.tsv")))?.lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with('#') {
                continue;
            }
            let x: Vec<_> = line.split('\t').collect();
            if x.len() != 5 {
                return Err("Malformed BiWFA output".into());
            }
            let ids: Vec<usize> = x[0]
                .strip_prefix('a')
                .ok_or("Bad BiWFA record ID")?
                .split('_')
                .map(|s| s.parse().map_err(|_| "Bad BiWFA record ID".into()))
                .collect::<Result<_>>()?;
            if ids.len() != 3 {
                return Err("Bad BiWFA record ID".into());
            }
            let (pi, ai, segment) = (ids[0], ids[1], ids[2]);
            let p = parents.get(pi).ok_or("Unknown BiWFA parent")?;
            let s = x[1]
                .parse::<usize>()
                .map_err(|_| "Bad BiWFA coordinate")?
                .checked_sub(p.start + 1)
                .ok_or("BiWFA output before source")?;
            let (r, a) = (
                if x[2] == "-" { "" } else { x[2] },
                if x[3] == "-" { "" } else { x[3] },
            );
            let &(gs, gr, _) = gaps
                .get(&(pi, ai))
                .and_then(|v| v.get(segment))
                .ok_or("Unknown BiWFA segment")?;
            if gs > s || s + r.len() > gs + gr.len() || p.alleles[0].get(s..s + r.len()) != Some(r)
            {
                return Err("BiWFA output violates source bounds/REF".into());
            }
            if x[4] != "0" && x[4] != "1" {
                return Err("Bad BiWFA complex flag".into());
            }
            result[pi][ai].extend(atomic_edits(s, r, a, x[4] == "1"));
        }
    }
    for (p, aa) in parents.iter().zip(&mut result) {
        for (ai, es) in aa.iter_mut().enumerate() {
            if matches!(p.kinds[ai].as_str(), "INV" | "COMPLEX") {
                *es = vec![Edit {
                    start: 1,
                    end: p.alleles[0].len() - 1,
                    alt: p.alleles[ai][1..p.alleles[ai].len() - 1].into(),
                    complex: true,
                }];
            }
            if local_replay(&p.alleles[0], es)? != p.alleles[ai] {
                return Err(format!(
                    "BiWFA decomposition failed source replay: {}:{ai}",
                    p.id
                ));
            }
        }
    }
    Ok(result)
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Partial {
    pub samples: Vec<usize>,
    pub edits: Vec<Edit>,
    pub positions: BTreeSet<usize>,
}
#[allow(clippy::too_many_arguments)]
fn partial(
    parents: &[Parent],
    complete: &Decomposition,
    known: &Coverage,
    reference: &str,
    chrom: &str,
    reference_path: &Path,
    work: &Path,
    jobs: usize,
    executable: &Adapter,
) -> Result<BTreeMap<(usize, usize), Partial>> {
    let mut requested = BTreeMap::new();
    let mut subset = vec![];
    for (pi, (p, aa)) in parents.iter().zip(complete).enumerate() {
        for (ai, es) in aa.iter().enumerate() {
            let selected: Vec<_> = p.carriers[ai]
                .iter()
                .copied()
                .filter(|&si| {
                    es.iter().any(|e| {
                        e.end > e.start
                            && !known.unique[si].covers(p.start + e.start, p.start + e.end)
                    })
                })
                .collect();
            if selected.is_empty() {
                continue;
            }
            let mut v = p.clone();
            v.alleles = vec![p.alleles[0].clone(), p.alleles[ai].clone()];
            v.anchors = vec![p.anchors[0].clone(), p.anchors[ai].clone()];
            v.kinds = vec!["REFERENCE".into(), p.kinds[ai].clone()];
            subset.push(v);
            requested.insert(
                (pi, ai),
                Partial {
                    samples: selected,
                    edits: vec![],
                    positions: BTreeSet::new(),
                },
            );
        }
    }
    if requested.is_empty() {
        return Ok(requested);
    }
    fs::create_dir(work).map_err(|e| e.to_string())?;
    let f = align(&subset, chrom, &work.join("forward"), jobs, executable)?;
    let r = reflect(
        &subset,
        align(
            &reverse_parents(&subset),
            chrom,
            &work.join("reverse"),
            jobs,
            executable,
        )?,
    );
    let (f, _) = normalize(
        &subset,
        &f,
        reference,
        chrom,
        reference_path,
        &work.join("normalized_forward"),
    )?;
    let (r, _) = normalize(
        &subset,
        &r,
        reference,
        chrom,
        reference_path,
        &work.join("normalized_reverse"),
    )?;
    for (((p, fw), rv), entry) in subset.iter().zip(f).zip(r).zip(requested.values_mut()) {
        entry.positions = disagreements(reference, p.start, &fw[1], &rv[1]);
        entry.edits = fw[1].clone();
    }
    Ok(requested)
}

pub struct Bins {
    rows: HashMap<usize, BTreeSet<Edit>>,
}
impl Bins {
    pub fn new(edits: &BTreeSet<Edit>) -> Self {
        let mut b = Self {
            rows: HashMap::new(),
        };
        for e in edits {
            b.insert(e.clone());
        }
        b
    }
    pub fn insert(&mut self, e: Edit) {
        for k in e.start / 1000..=e.end / 1000 {
            self.rows.entry(k).or_default().insert(e.clone());
        }
    }
    pub fn remove(&mut self, e: &Edit) {
        for k in e.start / 1000..=e.end / 1000 {
            if let Some(v) = self.rows.get_mut(&k) {
                v.remove(e);
            }
        }
    }
    pub fn local(&self, s: usize, e: usize) -> BTreeSet<Edit> {
        let mut out = BTreeSet::new();
        for k in s.saturating_sub(1) / 1000..=(e + 1) / 1000 {
            if let Some(v) = self.rows.get(&k) {
                out.extend(v.iter().filter(|x| intersects(s, e, x)).cloned());
            }
        }
        out
    }
}
fn cohere(
    models: &mut Models,
    reference: &str,
    known: &Coverage,
    complete: &BTreeSet<Edit>,
) -> Result<Vec<(usize, usize, usize, usize)>> {
    let loci: HashSet<_> = models
        .iter()
        .flatten()
        .filter(|e| e.complex)
        .map(|e| (e.start, e.end))
        .collect();
    let mut sequences: BTreeMap<(usize, usize), BTreeSet<String>> = BTreeMap::new();
    for e in models.iter().flatten() {
        if e.structural() && loci.contains(&(e.start, e.end)) && dna(&e.alt) {
            sequences
                .entry((e.start, e.end))
                .or_default()
                .insert(e.alt.clone());
        }
    }
    for e in complete {
        if e.complex && dna(&e.alt) {
            sequences
                .entry((e.start, e.end))
                .or_default()
                .insert(e.alt.clone());
        }
    }
    let mut loci: Vec<_> = sequences.keys().copied().collect();
    loci.sort_by_key(|&(s, e)| (std::cmp::Reverse(e - s), s, e));
    let mut promoted = vec![];
    for (si, model) in models.iter_mut().enumerate() {
        let mut bins = Bins::new(model);
        for &(s, t) in &loci {
            if !known.unique[si].covers(s, t) {
                continue;
            }
            let local = bins.local(s, t);
            if local.is_empty()
                || local.iter().any(|e| e.start < s || e.end > t)
                || local.len() == 1 && local.iter().next().unwrap().complex
            {
                continue;
            }
            let Ok(observed) = replay(reference, s, t, &local.iter().collect::<Vec<_>>()) else {
                continue;
            };
            if !sequences[&(s, t)].contains(&observed) {
                continue;
            }
            let replacement = Edit {
                start: s,
                end: t,
                alt: observed,
                complex: true,
            };
            for e in &local {
                model.remove(e);
                bins.remove(e);
            }
            model.insert(replacement.clone());
            bins.insert(replacement);
            promoted.push((si, s, t, local.len()));
        }
    }
    Ok(promoted)
}
pub fn physical_models(parents: &[Parent], samples: usize) -> Models {
    let mut models = vec![BTreeSet::new(); samples];
    for p in parents {
        let (s, e) = (p.start + 1, p.start + p.alleles[0].len() - 1);
        for (ai, a) in p.alleles.iter().enumerate().skip(1) {
            let x = Edit {
                start: s,
                end: e,
                alt: a[1..a.len() - 1].into(),
                complex: true,
            };
            for &si in &p.carriers[ai] {
                models[si].insert(x.clone());
            }
        }
    }
    models
}
pub fn validate_models(
    parents: &[Parent],
    models: &Models,
    reference: &str,
    known: &Coverage,
    ambiguous: &[Ranges],
) -> Result<Stats> {
    let mut stats = Stats::new();
    let physical = physical_models(parents, models.len());
    for (si, (raw, model)) in physical.iter().zip(models).enumerate() {
        let (left, right) = (Bins::new(raw), Bins::new(model));
        let mut done = HashSet::new();
        for unit in raw {
            let (mut s, mut t) = (unit.start, unit.end);
            let (original, current) = loop {
                let a = left.local(s, t);
                let b = right.local(s, t);
                let lo = a
                    .iter()
                    .chain(&b)
                    .map(|e| e.start)
                    .min()
                    .unwrap_or(s)
                    .min(s);
                let hi = a.iter().chain(&b).map(|e| e.end).max().unwrap_or(t).max(t);
                if lo == s && hi == t {
                    break (a, b);
                }
                (s, t) = (lo, hi);
            };
            if !done.insert((s, t)) {
                continue;
            }
            let (lo, hi) = if s < t {
                (s, t)
            } else {
                (s.saturating_sub(1), (t + 1).min(reference.len()))
            };
            if !known.unique[si].covers(lo, hi) || ambiguous[si].overlaps(lo, hi) {
                count(&mut stats, "source_windows_with_unresolved_territory", 1);
                continue;
            }
            if replay(reference, s, t, &original.iter().collect::<Vec<_>>())?
                != replay(reference, s, t, &current.iter().collect::<Vec<_>>())?
            {
                return Err(format!(
                    "Joint model changed physical sequence: sample {si}, {s}:{t}"
                ));
            }
            count(&mut stats, "joint_source_windows_replayed", 1);
        }
    }
    Ok(stats)
}

#[derive(Serialize, Deserialize)]
pub struct Prepared {
    pub parents: Vec<Parent>,
    pub models: Models,
    pub origins: OriginMap,
    pub event_kinds: KindMap,
    pub known: Coverage,
    pub ambiguous: Vec<Ranges>,
    pub ambiguity_rows: Vec<(String, usize, usize, usize)>,
}
#[allow(clippy::too_many_arguments)]
pub fn prepare(
    sources: &Sources,
    reference: &str,
    chrom: &str,
    reference_path: &Path,
    work: &Path,
    jobs: usize,
    executable: &Adapter,
    sequence_anchor: usize,
) -> Result<Prepared> {
    fs::create_dir(work).map_err(|e| e.to_string())?;
    let known = Coverage::new(&sources.coverage, &sources.uncertain, reference.len())?;
    let mut primitive = sources.parents.clone();
    for p in &mut primitive {
        for ai in 0..p.alleles.len() {
            p.anchors[ai] = anchors::context_runs(
                &p.alleles[0],
                &p.alleles[ai],
                p.anchors[ai].clone(),
                1,
                None,
            )?
            .into_iter()
            .flat_map(|r| r.payload)
            .collect();
        }
    }
    let mut parents = sources.parents.clone();
    for (p, q) in parents.iter_mut().zip(&primitive) {
        for ai in 0..p.alleles.len() {
            if p.kinds[ai] == "SEQUENCE" {
                p.anchors[ai] = anchors::refine(
                    &p.alleles[0],
                    &p.alleles[ai],
                    &q.anchors[ai],
                    sequence_anchor,
                )?;
            }
        }
    }
    let forward = decompose(&parents)?;
    let raw_reverse = reflect(&parents, decompose(&reverse_parents(&parents))?);
    let (mut edits, ff) = normalize(
        &parents,
        &forward,
        reference,
        chrom,
        reference_path,
        &work.join("normalized_forward"),
    )?;
    let (mut reverse, rf) = normalize(
        &parents,
        &raw_reverse,
        reference,
        chrom,
        reference_path,
        &work.join("normalized_reverse"),
    )?;
    let mut failed: Keys = ff.union(&rf).copied().collect();
    let (sf, sr, compound) = preserve(
        &parents,
        &forward,
        &raw_reverse,
        &edits,
        &reverse,
        &failed,
        reference,
    )?;
    if compound > 0 {
        let (f, ff) = normalize(
            &parents,
            &sf,
            reference,
            chrom,
            reference_path,
            &work.join("compound_forward"),
        )?;
        edits = f;
        let (r, rf) = normalize(
            &parents,
            &sr,
            reference,
            chrom,
            reference_path,
            &work.join("compound_reverse"),
        )?;
        reverse = r;
        let mut unresolved: Keys = ff.union(&rf).copied().collect();
        for (pi, p) in parents.iter().enumerate() {
            for ai in 1..p.alleles.len() {
                if !disagreements(reference, p.start, &edits[pi][ai], &reverse[pi][ai]).is_empty() {
                    unresolved.insert((pi, ai));
                }
            }
        }
        if !unresolved.is_empty() {
            let mut subset = vec![];
            let mut raw = vec![];
            for &(pi, ai) in &unresolved {
                let p = &parents[pi];
                let mut q = p.clone();
                q.alleles = vec![p.alleles[0].clone(), p.alleles[ai].clone()];
                subset.push(q);
                raw.push(vec![
                    vec![],
                    vec![Edit {
                        start: 1,
                        end: p.alleles[0].len() - 1,
                        alt: p.alleles[ai][1..p.alleles[ai].len() - 1].into(),
                        complex: true,
                    }],
                ]);
            }
            let (full, fail) = normalize(
                &subset,
                &raw,
                reference,
                chrom,
                reference_path,
                &work.join("compound_parents"),
            )?;
            if !fail.is_empty() {
                return Err("Complete parent failed normalization".into());
            }
            for (&(pi, ai), aa) in unresolved.iter().zip(full) {
                edits[pi][ai] = aa[1].clone();
                reverse[pi][ai] = aa[1].clone();
            }
        }
        failed.clear();
    }
    io::progress(format!(
        "[projection] preserved {compound} complete complex gaps"
    ));
    for (p, (fw, rv)) in parents.iter().zip(edits.iter_mut().zip(&mut reverse)) {
        for ai in 1..fw.len() {
            fw[ai] = close_touching(&p.alleles[0], &fw[ai])?;
            rv[ai] = close_touching(&p.alleles[0], &rv[ai])?;
        }
    }
    let partial = partial(
        &primitive,
        &edits,
        &known,
        reference,
        chrom,
        reference_path,
        &work.join("partial_projection"),
        jobs,
        executable,
    )?;
    let mut masks = vec![vec![]; sources.samples.len()];
    for entry in partial.values() {
        for &si in &entry.samples {
            masks[si].extend(entry.positions.iter().map(|&p| (p, p + 1)));
        }
    }
    let mut ambiguity_rows = vec![];
    for (pi, p) in parents.iter().enumerate() {
        for ai in 0..p.alleles.len() {
            let mut positions = disagreements(reference, p.start, &edits[pi][ai], &reverse[pi][ai]);
            if failed.contains(&(pi, ai)) {
                positions.extend(p.start..p.start + p.alleles[0].len());
            }
            for &si in &p.carriers[ai] {
                masks[si].extend(positions.iter().map(|&p| (p, p + 1)));
            }
            if !positions.is_empty() {
                ambiguity_rows.push((p.id.clone(), ai, positions.len(), p.carriers[ai].len()));
            }
        }
    }
    let mut models = vec![BTreeSet::new(); sources.samples.len()];
    let mut origins = OriginMap::new();
    let mut event_kinds = KindMap::new();
    let mut complete = BTreeSet::new();
    for (pi, p) in parents.iter().enumerate() {
        for (ai, es) in edits[pi].iter().enumerate() {
            let abs: BTreeSet<_> = es.iter().map(|e| absolute(e, p.start)).collect();
            complete.extend(abs.iter().cloned());
            let entry = partial.get(&(pi, ai));
            let fallback: BTreeSet<_> = entry
                .into_iter()
                .flat_map(|v| v.edits.iter())
                .map(|e| absolute(e, p.start))
                .collect();
            for &si in &p.carriers[ai] {
                models[si].extend(
                    if entry.is_some_and(|v| v.samples.contains(&si)) {
                        &fallback
                    } else {
                        &abs
                    }
                    .iter()
                    .cloned(),
                );
            }
            for e in abs.union(&fallback) {
                if e.complex {
                    event_kinds
                        .entry((e.start, e.end, e.alt.clone()))
                        .or_default()
                        .insert(p.kinds[ai].clone());
                }
                origins
                    .entry((
                        e.start,
                        e.end,
                        if e.structural() { "EDIT" } else { "SNP" }.into(),
                    ))
                    .or_default()
                    .insert(p.id.clone());
            }
        }
    }
    let promoted = cohere(&mut models, reference, &known, &complete)?;
    io::json(&work.join("representation_changes.json"), &promoted)?;
    for (si, mm) in models.iter().enumerate() {
        let mut previous: Option<&Edit> = None;
        for e in mm {
            if let Some(p) = previous {
                if e.start < p.end || e.start == e.end && e.end == p.start && p.start == p.end {
                    masks[si].push((
                        e.start.min(p.start).saturating_sub(1),
                        (e.end.max(p.end) + 1).min(reference.len()),
                    ));
                }
            }
            if previous.is_none_or(|p| e.end > p.end) {
                previous = Some(e);
            }
        }
    }
    let ambiguous: Vec<_> = masks
        .into_iter()
        .map(|v| Ranges::new(v, reference.len()))
        .collect::<Result<_>>()?;
    let model_validation = validate_models(&parents, &models, reference, &known, &ambiguous)?;
    io::progress(format!(
        "[projection] {} partial sources; {} model promotions; {}",
        partial.len(),
        promoted.len(),
        serde_json::to_string(&model_validation).map_err(|e| e.to_string())?
    ));
    Ok(Prepared {
        parents,
        models,
        origins,
        event_kinds,
        known,
        ambiguous,
        ambiguity_rows,
    })
}
