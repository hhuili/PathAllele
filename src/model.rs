//! Owned evidence, half-open coordinates and exact edit replay.
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub type EditInput = (usize, usize, String, bool);
pub type Intervals = Vec<(usize, usize)>;
pub type Result<T> = std::result::Result<T, String>;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct Edit {
    pub start: usize,
    pub end: usize,
    pub alt: String,
    pub complex: bool,
}

impl Edit {
    pub fn structural(&self) -> bool {
        self.complex || self.end - self.start != self.alt.len()
    }
    pub fn input(&self) -> EditInput {
        (self.start, self.end, self.alt.clone(), self.complex)
    }
}

pub fn intersects(s: usize, e: usize, x: &Edit) -> bool {
    if s == e {
        x.start <= s && s <= x.end
    } else if x.start == x.end {
        s <= x.start && x.start <= e
    } else {
        s < x.end && x.start < e
    }
}

pub fn dna(value: &str) -> bool {
    value
        .bytes()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

pub fn replay(reference: &str, s: usize, e: usize, edits: &[&Edit]) -> Result<String> {
    if s > e || e > reference.len() {
        return Err("Replay bounds outside reference".into());
    }
    let mut ordered = edits.to_vec();
    ordered.sort_unstable();
    ordered.dedup();
    let mut cursor = s;
    let mut previous_insertion = None;
    let mut result = String::with_capacity(e - s);
    for edit in ordered {
        if edit.start < cursor || edit.end > e || edit.start > edit.end {
            return Err("Overlapping/out-of-bounds edits".into());
        }
        if edit.start == edit.end && previous_insertion == Some(edit.start) {
            return Err("Multiple insertions at one boundary have no defined order".into());
        }
        result.push_str(&reference[cursor..edit.start]);
        result.push_str(&edit.alt);
        cursor = edit.end;
        previous_insertion = (edit.start == edit.end).then_some(edit.start);
    }
    result.push_str(&reference[cursor..e]);
    Ok(result)
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Ranges {
    pub rows: Intervals,
}

impl Ranges {
    pub fn new(mut rows: Intervals, length: usize) -> Result<Self> {
        if rows.iter().any(|&(s, e)| s > e || e > length) {
            return Err("Invalid interval bounds".into());
        }
        rows.sort_unstable();
        let mut merged: Intervals = Vec::new();
        for (s, e) in rows {
            if s == e {
                continue;
            }
            if let Some(last) = merged.last_mut() {
                if s <= last.1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        Ok(Self { rows: merged })
    }
    pub fn contains(&self, position: usize) -> bool {
        let i = self.rows.partition_point(|&(s, _)| s <= position);
        i > 0 && position < self.rows[i - 1].1
    }
    pub fn covers(&self, s: usize, e: usize) -> bool {
        if s == e {
            return true;
        }
        let i = self.rows.partition_point(|&(start, _)| start <= s);
        i > 0 && self.rows[i - 1].1 >= e
    }
    pub fn overlaps(&self, s: usize, e: usize) -> bool {
        let i = self.rows.partition_point(|&(_, end)| end <= s);
        s < e && i < self.rows.len() && self.rows[i].0 < e
    }
}

pub fn trim(start: usize, reference: &str, alternate: &str) -> Edit {
    let (r, a) = (reference.as_bytes(), alternate.as_bytes());
    let mut left = 0;
    while left < r.len().min(a.len()) && r[left] == a[left] {
        left += 1;
    }
    let (mut re, mut ae) = (r.len(), a.len());
    while re > left && ae > left && r[re - 1] == a[ae - 1] {
        re -= 1;
        ae -= 1;
    }
    Edit {
        start: start + left,
        end: start + re,
        alt: alternate[left..ae].into(),
        complex: false,
    }
}
pub fn atomic_edits(start: usize, reference: &str, alternate: &str, complex: bool) -> Vec<Edit> {
    let mut e = trim(start, reference, alternate);
    if e.start == e.end && e.alt.is_empty() {
        return vec![];
    }
    if e.end - e.start == 1 && e.alt.len() == 1 {
        return vec![e];
    }
    if complex {
        e.complex = true;
        return vec![e];
    }
    if e.end - e.start == e.alt.len() {
        return e
            .alt
            .bytes()
            .enumerate()
            .filter(|&(i, a)| (a != reference.as_bytes()[e.start - start + i]))
            .map(|(i, a)| Edit {
                start: e.start + i,
                end: e.start + i + 1,
                alt: (a as char).to_string(),
                complex: false,
            })
            .collect();
    }
    vec![e]
}
pub fn vcf_alleles(
    reference: &str,
    s: usize,
    e: usize,
    alts: &[String],
    snp: bool,
) -> Result<(usize, String, Vec<String>)> {
    if s > e || e > reference.len() {
        return Err("VCF coordinates outside reference".into());
    }
    if snp {
        return Ok((s + 1, reference[s..e].into(), alts.to_vec()));
    }
    if s > 0 {
        let anchor = &reference[s - 1..s];
        return Ok((
            s,
            reference[s - 1..e].into(),
            alts.iter().map(|a| format!("{anchor}{a}")).collect(),
        ));
    }
    if e < reference.len() {
        let anchor = &reference[e..e + 1];
        return Ok((
            1,
            reference[..e + 1].into(),
            alts.iter().map(|a| format!("{a}{anchor}")).collect(),
        ));
    }
    Err("Cannot pad an event covering the entire reference".into())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Coverage {
    pub unique: Vec<Ranges>,
    pub multiple: Vec<Ranges>,
}
impl Coverage {
    pub fn new(coverage: &[Intervals], uncertain: &[Intervals], length: usize) -> Result<Self> {
        if coverage.len() != uncertain.len() {
            return Err("Coverage sample dimensions differ".into());
        }
        let (mut unique, mut multiple) = (Vec::new(), Vec::new());
        for (a, b) in coverage.iter().zip(uncertain) {
            let mut events = std::collections::BTreeMap::<usize, [i64; 2]>::new();
            for (channel, rows) in [a, b].iter().enumerate() {
                for &(s, e) in *rows {
                    if s > e || e > length {
                        return Err("Coverage outside reference".into());
                    }
                    if s < e {
                        events.entry(s).or_default()[channel] += 1;
                        events.entry(e).or_default()[channel] -= 1;
                    }
                }
            }
            let (mut active, mut previous) = ([0i64; 2], 0usize);
            let (mut u, mut m) = (vec![], vec![]);
            for (pos, delta) in events {
                if previous < pos {
                    if active[1] > 0 || active[0] >= 2 {
                        m.push((previous, pos));
                    } else if active[0] == 1 {
                        u.push((previous, pos));
                    }
                }
                for k in 0..2 {
                    active[k] += delta[k];
                    if active[k] < 0 {
                        return Err("Unbalanced coverage".into());
                    }
                }
                previous = pos;
            }
            if active != [0, 0] {
                return Err("Unbalanced coverage sweep".into());
            }
            unique.push(Ranges::new(u, length)?);
            multiple.push(Ranges::new(m, length)?);
        }
        Ok(Self { unique, multiple })
    }
    pub fn state(&self, sample: usize, p: usize) -> u8 {
        if self.multiple[sample].contains(p) {
            2
        } else {
            self.unique[sample].contains(p) as u8
        }
    }
}

#[derive(Debug)]
pub struct Sample {
    pub edits: Vec<usize>,
    pub members: Vec<usize>,
    pub maximum_ends: Vec<usize>,
    pub unique: Ranges,
    pub multiple: Ranges,
    pub uncertain: Ranges,
    pub structural: Ranges,
    pub chain_starts: Vec<usize>,
    pub chain_ends: Vec<usize>,
}

impl Sample {
    pub fn overlapping(&self, pool: &[Edit], s: usize, e: usize) -> Vec<usize> {
        let lo = self.maximum_ends.partition_point(|&end| end < s);
        let hi = self.edits.partition_point(|&i| pool[i].start <= e);
        self.edits[lo.min(hi)..hi]
            .iter()
            .copied()
            .filter(|&i| intersects(s, e, &pool[i]))
            .collect()
    }
    pub fn base(&self, reference: &str, pool: &[Edit], position: usize) -> u8 {
        let mut result = reference.as_bytes()[position];
        for id in self.overlapping(pool, position, position + 1) {
            let edit = &pool[id];
            if edit.start <= position && position < edit.end {
                result = if edit.structural() {
                    b'-'
                } else {
                    edit.alt.as_bytes()[position - edit.start]
                };
            }
        }
        result
    }
    pub fn chain_covers(&self, s: usize, e: usize) -> bool {
        let i = self.chain_starts.partition_point(|&start| start <= s);
        i > 0 && self.chain_ends[i - 1] >= e
    }
}

pub struct Evidence {
    pub reference: String,
    pub pool: Vec<Edit>,
    pub lookup: HashMap<Edit, usize>,
    pub samples: Vec<Sample>,
}

impl Evidence {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: String,
        models: Vec<Vec<EditInput>>,
        unique: Vec<Intervals>,
        multiple: Vec<Intervals>,
        uncertain: Vec<Intervals>,
        structural: Vec<Intervals>,
        mut chains: Vec<Intervals>,
    ) -> Result<Self> {
        let length = reference.len();
        if !dna(&reference) || length == 0 {
            return Err("Reference must be nonempty ACGT".into());
        }
        let n = models.len();
        if [
            unique.len(),
            multiple.len(),
            uncertain.len(),
            structural.len(),
            chains.len(),
        ]
        .iter()
        .any(|&v| v != n)
        {
            return Err("Evidence sample dimensions differ".into());
        }
        let mut pool = Vec::new();
        let mut lookup = HashMap::new();
        let mut samples = Vec::with_capacity(n);
        for (si, model) in models.into_iter().enumerate() {
            let mut ids = Vec::with_capacity(model.len());
            for (start, end, alt, complex) in model {
                if start > end || end > length || !alt.is_ascii() {
                    return Err("Edit bounds or sequence encoding are invalid".into());
                }
                let edit = Edit {
                    start,
                    end,
                    alt,
                    complex,
                };
                let next = pool.len();
                let id = *lookup.entry(edit.clone()).or_insert_with(|| {
                    pool.push(edit);
                    next
                });
                ids.push(id);
            }
            ids.sort_unstable_by(|&a, &b| pool[a].cmp(&pool[b]));
            ids.dedup();
            let mut members = ids.clone();
            members.sort_unstable();
            let mut high = 0;
            let maximum_ends = ids
                .iter()
                .map(|&id| {
                    high = high.max(pool[id].end);
                    high
                })
                .collect();
            if chains[si].iter().any(|&(s, e)| s > e || e > length) {
                return Err("Chain bounds are invalid".into());
            }
            chains[si].sort_unstable();
            high = 0;
            let chain_starts = chains[si].iter().map(|&(s, _)| s).collect();
            let chain_ends = chains[si]
                .iter()
                .map(|&(_, e)| {
                    high = high.max(e);
                    high
                })
                .collect();
            samples.push(Sample {
                edits: ids,
                members,
                maximum_ends,
                unique: Ranges::new(unique[si].clone(), length)?,
                multiple: Ranges::new(multiple[si].clone(), length)?,
                uncertain: Ranges::new(uncertain[si].clone(), length)?,
                structural: Ranges::new(structural[si].clone(), length)?,
                chain_starts,
                chain_ends,
            });
        }
        Ok(Self {
            reference,
            pool,
            lookup,
            samples,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn closed_insertion_intersection_and_half_open_coverage() {
        let x = Edit {
            start: 2,
            end: 4,
            alt: "A".into(),
            complex: false,
        };
        assert!(intersects(2, 2, &x));
        assert!(intersects(4, 4, &x));
        assert!(!intersects(4, 5, &x));
        let r = Ranges::new(vec![(0, 2), (2, 4)], 5).unwrap();
        assert!(r.covers(0, 4));
        assert!(!r.contains(4));
    }
    #[test]
    fn replay_rejects_conflicting_insertions() {
        let a = Edit {
            start: 1,
            end: 1,
            alt: "T".into(),
            complex: false,
        };
        let b = Edit {
            start: 1,
            end: 1,
            alt: "G".into(),
            complex: false,
        };
        assert!(replay("AC", 0, 2, &[&a, &b]).is_err());
        assert_eq!(replay("AC", 0, 2, &[&a, &a]).unwrap(), "ATC");
    }
}
