//! Exact homology, occurrence support and ties in collinear chains.
use crate::{
    graph::{Graph, Node, NodeId},
    io::reverse_complement,
    model::{dna, Result},
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct Anchor {
    pub r: usize,
    pub q: usize,
    pub n: usize,
    pub node: Option<NodeId>,
}
#[derive(Clone, Debug)]
pub struct Run {
    pub r: usize,
    pub q: usize,
    pub n: usize,
    pub payload: Vec<Anchor>,
}
#[derive(Clone, Debug)]
pub struct MatchIndex {
    pub runs: Vec<(usize, usize, usize)>,
    pub starts: Vec<usize>,
    pub ends: Vec<usize>,
}
impl MatchIndex {
    pub fn new(runs: Vec<(usize, usize, usize)>) -> Self {
        let starts = runs.iter().map(|x| x.0).collect();
        let mut end = 0;
        let ends = runs
            .iter()
            .map(|x| {
                end = end.max(x.1);
                end
            })
            .collect();
        Self { runs, starts, ends }
    }
    pub fn locations(&self, q: usize, n: usize) -> HashSet<i64> {
        let need = 31.min(n);
        let mut i = self.starts.partition_point(|&s| s <= q + n - need);
        let mut found = HashSet::new();
        while i > 0 && self.ends[i - 1] >= q + need {
            i -= 1;
            let (qs, qe, rs) = self.runs[i];
            if (q + n).min(qe) >= q.max(qs) + need {
                found.insert(rs as i64 + q as i64 - qs as i64);
            }
        }
        found
    }
}

type Score = (usize, BigUint);
fn combine(a: &Score, b: &Score) -> Score {
    match a.0.cmp(&b.0) {
        std::cmp::Ordering::Greater => a.clone(),
        std::cmp::Ordering::Less => b.clone(),
        std::cmp::Ordering::Equal => (a.0, &a.1 + &b.1),
    }
}
fn scores(items: &[(usize, usize, usize)]) -> (Vec<Score>, Score) {
    let mut coords: Vec<_> = items.iter().map(|&(_, q, n)| q + n).collect();
    coords.sort_unstable();
    coords.dedup();
    let mut tree = vec![(0, BigUint::zero()); coords.len() + 1];
    let mut pending = BinaryHeap::new();
    let mut values: Vec<Score> = Vec::new();
    for (i, &(r, q, n)) in items.iter().enumerate() {
        while pending.peek().is_some_and(|Reverse((end, _, _))| *end <= r) {
            let Reverse((_, j, qe)) = pending.pop().unwrap();
            let mut k = coords.partition_point(|&v| v < qe) + 1;
            while k < tree.len() {
                tree[k] = combine(&tree[k], &values[j]);
                k += k & k.wrapping_neg();
            }
        }
        let mut best = (0, BigUint::one());
        let mut k = coords.partition_point(|&v| v <= q);
        while k > 0 {
            best = combine(&best, &tree[k]);
            k -= k & k.wrapping_neg();
        }
        values.push((best.0 + n, best.1));
        pending.push(Reverse((r + n, i, q + n)));
    }
    let total = values
        .iter()
        .fold((0, BigUint::zero()), |a, b| combine(&a, b));
    (values, total)
}
pub fn obligatory(items: &[(usize, usize, usize)]) -> Result<Vec<bool>> {
    if items.iter().any(|x| x.2 == 0) {
        return Err("Anchor lengths must be positive".into());
    }
    if items.is_empty() {
        return Ok(vec![]);
    }
    let mut order: Vec<_> = (0..items.len()).collect();
    order.sort_unstable_by_key(|&i| items[i]);
    let sorted: Vec<_> = order.iter().map(|&i| items[i]).collect();
    let (f, total) = scores(&sorted);
    let limit_r = sorted.iter().map(|x| x.0 + x.2).max().unwrap();
    let limit_q = sorted.iter().map(|x| x.1 + x.2).max().unwrap();
    let mut reflected: Vec<_> = sorted
        .iter()
        .enumerate()
        .map(|(i, &(r, q, n))| ((limit_r - r - n, limit_q - q - n, n), i))
        .collect();
    reflected.sort_unstable();
    let (b, _) = scores(&reflected.iter().map(|x| x.0).collect::<Vec<_>>());
    let mut back = vec![(0, BigUint::zero()); items.len()];
    for (j, item) in reflected.iter().enumerate() {
        back[item.1] = b[j].clone();
    }
    let mut selected = vec![false; items.len()];
    for (i, &original) in order.iter().enumerate() {
        selected[original] =
            f[i].0 + back[i].0 - sorted[i].2 == total.0 && &f[i].1 * &back[i].1 == total.1;
    }
    Ok(selected)
}

/// An exact forward k-mer counter. Ambiguous bases reset the rolling key.
pub struct SeedCounts {
    pub k: usize,
    counts: HashMap<u64, u32>,
}
fn encode(word: &str) -> Option<u64> {
    let mut value = 0;
    for b in word.bytes() {
        value = (value << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(value)
}
impl SeedCounts {
    pub fn new(sequence: &str, k: usize) -> Self {
        let mut counts = HashMap::new();
        if k == 0 || k > 31 {
            return Self { k, counts };
        }
        let mask = (1u64 << (2 * k)) - 1;
        let (mut value, mut valid) = (0u64, 0usize);
        for b in sequence.bytes() {
            let bits = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    valid = 0;
                    value = 0;
                    continue;
                }
            };
            value = ((value << 2) | bits) & mask;
            valid += 1;
            if valid >= k {
                let c = counts.entry(value).or_insert(0u32);
                *c = c.saturating_add(1);
            }
        }
        Self { k, counts }
    }
    pub fn count(&self, word: &str) -> Option<u32> {
        if word.len() != self.k || self.k > 31 || self.k == 0 {
            return None;
        }
        Some(
            encode(word)
                .and_then(|x| self.counts.get(&x).copied())
                .unwrap_or(0),
        )
    }
}
pub struct Uniqueness<'a> {
    pub reference: &'a str,
    pub query: &'a str,
    pub reference_counts: &'a [SeedCounts],
    pub query_counts: Vec<SeedCounts>,
}
impl<'a> Uniqueness<'a> {
    pub fn new(reference: &'a str, query: &'a str, reference_counts: &'a [SeedCounts]) -> Self {
        let query_counts = reference_counts
            .iter()
            .map(|c| SeedCounts::new(query, c.k))
            .collect();
        Self {
            reference,
            query,
            reference_counts,
            query_counts,
        }
    }
    fn counts(&self, word: &str) -> Option<(u32, u32)> {
        self.reference_counts
            .iter()
            .zip(&self.query_counts)
            .find_map(|(a, b)| Some((a.count(word)?, b.count(word)?)))
    }
    pub fn oriented(&self, word: &str) -> bool {
        let rc = reverse_complement(word);
        if !dna(word) || word == rc {
            return false;
        }
        if let (Some((r, q)), Some((rr, qr))) = (self.counts(word), self.counts(&rc)) {
            return r == 1 && q == 1 && rr == 0 && qr == 0;
        }
        [self.reference, self.query]
            .iter()
            .all(|s| occurs_once(s, word) && !s.contains(&rc))
    }
    /// The reference word occurs once, in reverse orientation only in query.
    pub fn inverted(&self, word: &str) -> bool {
        let rc = reverse_complement(word);
        if !dna(word) || word == rc {
            return false;
        }
        if let (Some((r, q)), Some((rr, qr))) = (self.counts(word), self.counts(&rc)) {
            return r == 1 && q == 0 && rr == 0 && qr == 1;
        }
        occurs_once(self.reference, word)
            && !self.reference.contains(&rc)
            && occurs_once(self.query, &rc)
            && !self.query.contains(word)
    }
    pub fn forward(&self, r: usize, q: usize, n: usize) -> bool {
        // A unique contained word proves uniqueness of an already exact match.
        for c in self.reference_counts {
            let k = c.k;
            if k > n || k > 31 || k == 0 {
                continue;
            }
            for offset in [0, n - k] {
                if self.counts(&self.reference[r + offset..r + offset + k]) == Some((1, 1)) {
                    return true;
                }
            }
        }
        let word = &self.reference[r..r + n];
        self.reference.find(word) == Some(r)
            && !self.reference[r + 1..].contains(word)
            && self.query.find(word) == Some(q)
            && !self.query[q + 1..].contains(word)
    }
}
pub fn occurs_once(sequence: &str, word: &str) -> bool {
    sequence
        .find(word)
        .is_some_and(|p| !sequence[p + 1..].contains(word))
}

/// Maximal unique exact matches. Counts and positions are exact, including ties.
pub fn exact_runs(reference: &str, alternate: &str, minimum: usize) -> Result<Vec<Anchor>> {
    if minimum == 0 {
        return Err("Sequence anchor length must be positive".into());
    }
    if reference.len().min(alternate.len()) < minimum || !dna(reference) || !dna(alternate) {
        return Ok(vec![]);
    }
    fn index(s: &str, k: usize) -> HashMap<&str, (usize, usize)> {
        let mut x = HashMap::new();
        for i in 0..=s.len() - k {
            let v = x.entry(&s[i..i + k]).or_insert((0, i));
            v.0 += 1;
        }
        x
    }
    let ri = index(reference, minimum);
    let qi = index(alternate, minimum);
    let mut seeds = Vec::new();
    for (word, &(count, q)) in &qi {
        if count == 1 {
            if let Some(&(1, r)) = ri.get(word) {
                seeds.push((q, r));
            }
        }
    }
    seeds.sort_unstable();
    let mut merged: Vec<(usize, usize, usize)> = vec![];
    for (q, r) in seeds {
        if let Some(a) = merged.last_mut() {
            if r as i64 - q as i64 == a.0 as i64 - a.1 as i64 && q <= a.1 + a.2 {
                a.2 = a.2.max(q + minimum - a.1);
                continue;
            }
        }
        merged.push((r, q, minimum));
    }
    let (rr, ra) = (reverse_complement(reference), reverse_complement(alternate));
    let mut matches = Vec::new();
    for (mut r, mut q, mut n) in merged {
        while r > 0 && q > 0 && reference.as_bytes()[r - 1] == alternate.as_bytes()[q - 1] {
            r -= 1;
            q -= 1;
            n += 1;
        }
        while r + n < reference.len()
            && q + n < alternate.len()
            && reference.as_bytes()[r + n] == alternate.as_bytes()[q + n]
        {
            n += 1;
        }
        let word = &reference[r..r + n];
        if rr.contains(word) || ra.contains(word) {
            continue;
        }
        if reference.find(word) != Some(r)
            || reference[r + 1..].contains(word)
            || alternate.find(word) != Some(q)
            || alternate[q + 1..].contains(word)
        {
            continue;
        }
        matches.push((r, q, n));
    }
    matches.sort_unstable();
    matches.dedup();
    let chosen = obligatory(&matches)?;
    Ok(matches
        .into_iter()
        .zip(chosen)
        .filter_map(|((r, q, n), yes)| {
            yes.then_some(Anchor {
                r,
                q,
                n,
                node: None,
            })
        })
        .collect())
}
pub fn refine(
    reference: &str,
    alternate: &str,
    anchors: &[Anchor],
    minimum: usize,
) -> Result<Vec<Anchor>> {
    let mut ordered = vec![Anchor {
        r: 0,
        q: 0,
        n: 1,
        node: None,
    }];
    ordered.extend_from_slice(anchors);
    ordered.push(Anchor {
        r: reference.len() - 1,
        q: alternate.len() - 1,
        n: 1,
        node: None,
    });
    let mut result = anchors.to_vec();
    let (mut re, mut qe) = (0, 0);
    for a in ordered {
        if a.r < re || a.q < qe || reference.get(a.r..a.r + a.n) != alternate.get(a.q..a.q + a.n) {
            return Err("Invalid physical anchors".into());
        }
        if reference[re..a.r] != alternate[qe..a.q] {
            for mut x in exact_runs(&reference[re..a.r], &alternate[qe..a.q], minimum)? {
                x.r += re;
                x.q += qe;
                result.push(x);
            }
        }
        re = a.r + a.n;
        qe = a.q + a.n;
    }
    result.sort_by_key(|a| (a.r, a.q));
    Ok(result)
}
pub fn gaps<'a>(
    reference: &'a str,
    alternate: &'a str,
    anchors: &[Anchor],
) -> Result<Vec<(usize, &'a str, &'a str)>> {
    if reference.is_empty()
        || alternate.is_empty()
        || reference.as_bytes()[0] != alternate.as_bytes()[0]
        || reference.as_bytes().last() != alternate.as_bytes().last()
    {
        return Err("Source allele lacks identical padding".into());
    }
    let mut ordered = vec![(0, 0, 1)];
    ordered.extend(anchors.iter().map(|a| (a.r, a.q, a.n)));
    ordered.push((reference.len() - 1, alternate.len() - 1, 1));
    let (mut re, mut qe) = (0, 0);
    let mut out = vec![];
    for (r, q, n) in ordered {
        if r < re || q < qe || n == 0 {
            return Err("Overlapping/noncollinear native graph anchors".into());
        }
        let a = reference.get(r..r + n);
        if a.is_none() || a != alternate.get(q..q + n) {
            return Err("Native graph anchor sequence mismatch".into());
        }
        if reference[re..r] != alternate[qe..q] {
            out.push((re, &reference[re..r], &alternate[qe..q]));
        }
        re = r + n;
        qe = q + n;
    }
    if re != reference.len() || qe != alternate.len() {
        return Err("Incomplete native anchors".into());
    }
    Ok(out)
}

pub fn context_runs(
    reference: &str,
    query: &str,
    mut anchors: Vec<Anchor>,
    minimum: usize,
    unique: Option<&Uniqueness>,
) -> Result<Vec<Run>> {
    anchors.sort_by_key(|a| (a.q, a.r));
    let mut runs: Vec<Run> = Vec::new();
    for a in anchors {
        if reference.get(a.r..a.r + a.n) != query.get(a.q..a.q + a.n)
            || reference.get(a.r..a.r + a.n).is_none()
        {
            return Err("Graph homology differs from sequence".into());
        }
        if let Some(last) = runs.last_mut() {
            if a.r == last.r + last.n && a.q == last.q + last.n {
                last.n += a.n;
                last.payload.push(a);
                continue;
            }
        }
        runs.push(Run {
            r: a.r,
            q: a.q,
            n: a.n,
            payload: vec![a],
        });
    }
    runs.retain(|a| {
        a.n >= minimum
            && dna(&reference[a.r..a.r + a.n])
            && unique.is_none_or(|u| u.forward(a.r, a.q, a.n))
    });
    if runs.windows(2).any(|w| w[0].q + w[0].n > w[1].q) {
        return Err("Overlapping query anchors".into());
    }
    let selected = obligatory(&runs.iter().map(|a| (a.r, a.q, a.n)).collect::<Vec<_>>())?;
    Ok(runs
        .into_iter()
        .zip(selected)
        .filter_map(|(a, ok)| ok.then_some(a))
        .collect())
}
pub fn occurrence_seeds(
    graph: &Graph,
    walk: &[NodeId],
    offsets: &[usize],
    index: &MatchIndex,
) -> Vec<Anchor> {
    let mut result = Vec::new();
    for (i, &node) in walk.iter().enumerate() {
        let q = offsets[i];
        let n = graph.nodes[node.unsigned_abs() as usize].sequence.len();
        let loc = index.locations(q, n);
        if loc.len() == 1 {
            let r = *loc.iter().next().unwrap();
            if r >= 0 && graph.ref_occurrences.contains_key(&(node, r as usize)) {
                result.push(Anchor {
                    r: r as usize,
                    q,
                    n,
                    node: Some(node),
                });
            }
        }
    }
    result
}
#[allow(clippy::too_many_arguments)]
pub fn bounded(
    reference: &str,
    query: &str,
    rw: &[NodeId],
    qw: &[NodeId],
    nodes: &[Node],
    rbase: usize,
    qbase: usize,
    index: &MatchIndex,
    ro: usize,
    qo: usize,
) -> Result<Vec<Anchor>> {
    fn positions(walk: &[NodeId], nodes: &[Node]) -> HashMap<NodeId, Vec<usize>> {
        let mut map: HashMap<NodeId, Vec<usize>> = HashMap::new();
        let mut p = 0;
        for &n in walk {
            map.entry(n).or_default().push(p);
            p += nodes[n.unsigned_abs() as usize].sequence.len();
        }
        map
    }
    let ri = positions(rw, nodes);
    let qi = positions(qw, nodes);
    let mut anchors = Vec::new();
    for (node, positions) in qi {
        if positions.len() != 1 {
            continue;
        }
        if let Some(targets) = ri.get(&node) {
            if targets.len() == 1 {
                let a = Anchor {
                    r: targets[0] + rbase,
                    q: positions[0] + qbase,
                    n: nodes[node.unsigned_abs() as usize].sequence.len(),
                    node: Some(node),
                };
                if index.locations(a.q + qo, a.n) == HashSet::from([(a.r + ro) as i64]) {
                    anchors.push(a);
                }
            }
        }
    }
    Ok(context_runs(reference, query, anchors, 1, None)?
        .into_iter()
        .flat_map(|r| r.payload)
        .collect())
}
pub fn extend_terminal(
    reference: &str,
    query: &str,
    mut runs: Vec<Run>,
    index: &MatchIndex,
    minimum: usize,
    unique: &Uniqueness,
) -> Result<Vec<Run>> {
    let mut events: BTreeMap<usize, Vec<(i64, i32)>> = BTreeMap::new();
    for &(q, e, r) in &index.runs {
        events.entry(q).or_default().push((r as i64 - q as i64, 1));
        events.entry(e).or_default().push((r as i64 - q as i64, -1));
    }
    let mut active = HashMap::new();
    let mut segments: Vec<Anchor> = Vec::new();
    let mut previous = None;
    for (q, deltas) in events {
        if let Some(p) = previous {
            if q > p && active.len() == 1 {
                let r = p as i64 + *active.keys().next().unwrap();
                if r < 0 {
                    return Err("Negative alignment coordinate".into());
                }
                let r = r as usize;
                let n = q - p;
                if let Some(last) = segments.last_mut() {
                    if last.r + last.n == r && last.q + last.n == p {
                        last.n += n;
                    } else {
                        segments.push(Anchor {
                            r,
                            q: p,
                            n,
                            node: None,
                        });
                    }
                } else {
                    segments.push(Anchor {
                        r,
                        q: p,
                        n,
                        node: None,
                    });
                }
            }
        }
        for (offset, delta) in deltas {
            let count = active.entry(offset).or_insert(0);
            *count += delta;
            if *count == 0 {
                active.remove(&offset);
            }
        }
        previous = Some(q);
    }
    if runs.is_empty() {
        return context_runs(reference, query, segments, minimum, Some(unique));
    }
    let first = &runs[0];
    let last = runs.last().unwrap();
    let mut before = Vec::new();
    let mut after = Vec::new();
    for a in segments {
        let n = (a.n as i64)
            .min(first.r as i64 - a.r as i64)
            .min(first.q as i64 - a.q as i64);
        if n > 0 {
            before.push(Anchor {
                r: a.r,
                q: a.q,
                n: n as usize,
                node: None,
            });
        }
        let skip = (last.r + last.n)
            .saturating_sub(a.r)
            .max((last.q + last.n).saturating_sub(a.q));
        if skip < a.n {
            after.push(Anchor {
                r: a.r + skip,
                q: a.q + skip,
                n: a.n - skip,
                node: None,
            });
        }
    }
    let mut out = context_runs(reference, query, before, minimum, Some(unique))?;
    out.append(&mut runs);
    out.extend(context_runs(
        reference,
        query,
        after,
        minimum,
        Some(unique),
    )?);
    Ok(out)
}
