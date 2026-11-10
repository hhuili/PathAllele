//! Physical mapping ownership and conservative graph-only source components.
use crate::{
    anchors::{self, Anchor, MatchIndex, Run, SeedCounts, Uniqueness},
    graph::{offsets, Graph, NodeId},
    io,
    model::{dna, Intervals, Ranges, Result},
};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Alignment {
    pub identity: String,
    pub reverse: bool,
    pub qs: usize,
    pub qe: usize,
    pub rs: usize,
    pub re: usize,
    pub native_qs: usize,
    pub native_qe: usize,
    pub index: MatchIndex,
}
#[derive(Clone, Debug)]
pub struct Block {
    pub identity: String,
    pub kind: String,
    pub native_kind: Option<String>,
    pub reverse: bool,
    pub runs: Vec<Run>,
    pub record: Option<Alignment>,
}
pub struct Oriented {
    pub walk: Vec<NodeId>,
    pub offsets: Vec<usize>,
    pub query: String,
}
pub struct Plan {
    pub blocks: Vec<Block>,
    pub missing: Intervals,
    pub oriented: [Oriented; 2],
}

pub fn read_alignments(path: &Path, min_mapq: u32) -> Result<HashMap<String, Vec<Alignment>>> {
    let mut result: HashMap<String, Vec<Alignment>> = HashMap::new();
    for (i, line) in io::reader(path)?.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let f: Vec<_> = line.trim_end().split('\t').collect();
        if f.len() < 12 {
            return Err("Malformed PAF row".into());
        }
        let number = |i: usize| {
            f[i].parse::<usize>()
                .map_err(|_| "Invalid PAF coordinate".to_string())
        };
        if number(11)? < min_mapq as usize || !f[12..].contains(&"tp:A:P") {
            continue;
        }
        let cigar = f[12..]
            .iter()
            .find_map(|v| v.strip_prefix("cg:Z:"))
            .ok_or("PAF requires an =/X CIGAR")?;
        let reverse = match f[4] {
            "+" => false,
            "-" => true,
            _ => return Err("Invalid PAF orientation".into()),
        };
        let qlen = number(1)?;
        let (qs, qe, rs, re) = (number(2)?, number(3)?, number(7)?, number(8)?);
        if qs > qe || qe > qlen || rs > re {
            return Err("Invalid PAF bounds".into());
        }
        let mut q = if reverse { qlen - qe } else { qs };
        let start = q;
        let mut r = rs;
        let mut runs = Vec::new();
        let (mut count, mut digits) = (0usize, false);
        for op in cigar.bytes() {
            if op.is_ascii_digit() {
                count = count
                    .checked_mul(10)
                    .and_then(|n| n.checked_add((op - b'0') as usize))
                    .ok_or("CIGAR length overflow")?;
                digits = true;
                continue;
            }
            if !digits || !b"IDNSHP=X".contains(&op) {
                return Err("Exact =/X CIGAR operations are required".into());
            }
            if op == b'=' {
                runs.push((q, q.checked_add(count).ok_or("CIGAR overflow")?, r));
            }
            if b"IS=X".contains(&op) {
                q = q.checked_add(count).ok_or("CIGAR overflow")?;
            }
            if b"DN=X".contains(&op) {
                r = r.checked_add(count).ok_or("CIGAR overflow")?;
            }
            count = 0;
            digits = false;
        }
        if digits || q != if reverse { qlen - qs } else { qe } || r != re {
            return Err("PAF CIGAR does not consume declared intervals".into());
        }
        result.entry(f[0].into()).or_default().push(Alignment {
            identity: format!("M{:08}", i + 1),
            reverse,
            qs: start,
            qe: q,
            rs,
            re,
            native_qs: qs,
            native_qe: qe,
            index: MatchIndex::new(runs),
        });
    }
    Ok(result)
}
pub fn intervals(mask: &[bool]) -> Intervals {
    let mut rows = Vec::new();
    let mut start = None;
    for (i, &v) in mask.iter().enumerate() {
        if v && start.is_none() {
            start = Some(i)
        } else if !v {
            if let Some(s) = start.take() {
                rows.push((s, i));
            }
        }
    }
    if let Some(s) = start {
        rows.push((s, mask.len()));
    }
    rows
}
pub fn copy_intervals(
    graph: &Graph,
    context: &HashMap<String, Vec<Alignment>>,
) -> Result<Vec<Intervals>> {
    let mut component = 0usize;
    let mut output = Vec::new();
    for paths in &graph.sample_paths {
        let mut seen = vec![0usize; graph.reference.len()];
        let mut conflict = vec![false; seen.len()];
        for &pi in paths {
            let mut records: Vec<_> = context
                .get(&graph.paths[pi].name)
                .into_iter()
                .flatten()
                .collect();
            records.sort_by_key(|r| (r.native_qs, r.native_qe, r.rs, r.re));
            let mut end = None;
            for r in records {
                if end.is_none_or(|e| r.native_qs >= e) {
                    component += 1;
                }
                end = Some(end.unwrap_or(0).max(r.native_qe));
                if r.re > seen.len() {
                    return Err("PAF reference interval exceeds reference".into());
                }
                for p in r.rs..r.re {
                    if seen[p] == 0 {
                        seen[p] = component;
                    } else if seen[p] != component {
                        conflict[p] = true;
                    }
                }
            }
        }
        output.push(intervals(&conflict));
    }
    Ok(output)
}
fn graph_chains(
    graph: &Graph,
    view: &Oriented,
    minimum: usize,
    max_gap: usize,
    unique: &Uniqueness,
) -> Vec<Vec<Run>> {
    let mut counts = HashMap::new();
    for &n in &view.walk {
        *counts.entry(n.abs()).or_insert(0usize) += 1;
    }
    let mut chunks = Vec::new();
    let mut run: Vec<Run> = Vec::new();
    let mut q = 0;
    for &node in &view.walk {
        let sequence = &graph.nodes[node.unsigned_abs() as usize].sequence;
        let n = sequence.len();
        let Some(targets) = graph.ref_places.get(&node.abs()) else {
            q += n;
            continue;
        };
        let valid = targets.len() == 1
            && counts[&node.abs()] == 1
            && graph.refwalk()[targets[0]] == node
            && dna(sequence);
        let r = if valid {
            graph.ref_offsets[targets[0]]
        } else {
            0
        };
        if let Some(last) = run.last() {
            if !valid
                || r < last.r + last.n
                || r.saturating_sub(last.r + last.n).max(q - last.q - last.n) > max_gap
                || !dna(&view.query[last.q + last.n..q])
            {
                chunks.push(std::mem::take(&mut run));
            }
        }
        if valid {
            let a = Anchor {
                r,
                q,
                n,
                node: Some(node),
            };
            if let Some(last) = run.last_mut() {
                if (r, q) == (last.r + last.n, last.q + last.n) {
                    last.n += n;
                    last.payload.push(a);
                    q += n;
                    continue;
                }
            }
            run.push(Run {
                r,
                q,
                n,
                payload: vec![a],
            });
        }
        q += n;
    }
    if !run.is_empty() {
        chunks.push(run);
    }
    chunks
        .into_iter()
        .filter_map(|chunk| {
            let boundaries: Vec<_> = chunk
                .iter()
                .enumerate()
                .filter(|(_, a)| {
                    a.n >= minimum
                        && unique.oriented(&view.query[a.q..a.q + minimum])
                        && unique.oriented(&view.query[a.q + a.n - minimum..a.q + a.n])
                })
                .map(|(i, _)| i)
                .collect();
            if boundaries.is_empty() {
                None
            } else {
                Some(chunk[boundaries[0]..=*boundaries.last().unwrap()].to_vec())
            }
        })
        .collect()
}
pub fn clip(
    runs: &[Run],
    left: usize,
    right: usize,
    anchor: usize,
    unique: &Uniqueness,
) -> Vec<Run> {
    let clipped: Vec<_> = runs
        .iter()
        .filter_map(|a| {
            let lo = left.max(a.q);
            let hi = right.min(a.q + a.n);
            (lo < hi).then(|| {
                let b = Anchor {
                    r: a.r + lo - a.q,
                    q: lo,
                    n: hi - lo,
                    node: None,
                };
                Run {
                    r: b.r,
                    q: b.q,
                    n: b.n,
                    payload: vec![b],
                }
            })
        })
        .collect();
    let valid: Vec<_> = clipped
        .iter()
        .enumerate()
        .filter(|(_, a)| {
            a.n >= anchor
                && unique.oriented(&unique.query[a.q..a.q + anchor])
                && unique.oriented(&unique.query[a.q + a.n - anchor..a.q + a.n])
        })
        .map(|(i, _)| i)
        .collect();
    if valid.is_empty() {
        vec![]
    } else {
        clipped[valid[0]..=*valid.last().unwrap()].to_vec()
    }
}
pub fn native_bounds(block: &Block, length: usize) -> (usize, usize) {
    let a = &block.runs[0];
    let b = block.runs.last().unwrap();
    if block.reverse {
        (length - b.q - b.n, length - a.q)
    } else {
        (a.q, b.q + b.n)
    }
}

pub fn plan(
    graph: &Graph,
    path: usize,
    records: &[Alignment],
    minimum: usize,
    mapping_k: usize,
    max_gap: usize,
    ref_counts: &[SeedCounts],
) -> Result<Plan> {
    let sequence = graph.spell(path);
    let length = sequence.len();
    let walk = graph.paths[path].nodes.clone();
    let reverse_walk: Vec<_> = walk.iter().rev().map(|n| -n).collect();
    let reverse_query = io::reverse_complement(&sequence);
    let oriented = [
        Oriented {
            offsets: offsets(&graph.nodes, &walk),
            walk,
            query: sequence,
        },
        Oriented {
            offsets: offsets(&graph.nodes, &reverse_walk),
            walk: reverse_walk,
            query: reverse_query,
        },
    ];
    let unique = [
        Uniqueness::new(&graph.reference, &oriented[0].query, ref_counts),
        Uniqueness::new(&graph.reference, &oriented[1].query, ref_counts),
    ];
    let mut mapped = Vec::new();
    for strand in 0..2 {
        let view = &oriented[strand];
        for record in records.iter().filter(|r| r.reverse == (strand == 1)) {
            let seeds = anchors::occurrence_seeds(graph, &view.walk, &view.offsets, &record.index)
                .into_iter()
                .filter(|a| {
                    record.qs <= a.q
                        && a.q + a.n <= record.qe
                        && record.rs <= a.r
                        && a.r + a.n <= record.re
                })
                .collect();
            let runs = anchors::context_runs(
                &graph.reference,
                &view.query,
                seeds,
                minimum,
                Some(&unique[strand]),
            )?;
            let runs = anchors::extend_terminal(
                &graph.reference,
                &view.query,
                runs,
                &record.index,
                mapping_k,
                &unique[strand],
            )?;
            if !runs.is_empty() {
                mapped.push(Block {
                    identity: record.identity.clone(),
                    kind: "mapped".into(),
                    native_kind: None,
                    reverse: strand == 1,
                    runs,
                    record: Some(record.clone()),
                });
            }
        }
    }
    let mut owners = vec![0u32; length];
    for b in &mapped {
        let (s, e) = native_bounds(b, length);
        for v in &mut owners[s..e] {
            *v += 1;
        }
    }
    let mut accepted = Vec::new();
    let mut occupied = vec![false; length];
    for block in mapped {
        let strand = block.reverse as usize;
        let rs = &block.runs;
        let lo = rs[0].q;
        let hi = rs.last().unwrap().q + rs.last().unwrap().n;
        let allowed: Vec<_> = (lo..hi)
            .map(|q| owners[if block.reverse { length - 1 - q } else { q }] == 1)
            .collect();
        for (a, b) in intervals(&allowed) {
            let (left, right) = (a + lo, b + lo);
            let runs = if (left, right) == (lo, hi) {
                rs.clone()
            } else {
                clip(rs, left, right, mapping_k, &unique[strand])
            };
            if runs.is_empty() {
                continue;
            }
            let item = Block {
                runs,
                ..block.clone()
            };
            let (x, y) = native_bounds(&item, length);
            occupied[x..y].fill(true);
            accepted.push(item);
        }
    }
    let available: Vec<_> = occupied
        .iter()
        .zip(&owners)
        .map(|(&o, &n)| !o && n < 2)
        .collect();
    let mut supplements = Vec::new();
    for strand in 0..2 {
        for (j, runs) in graph_chains(graph, &oriented[strand], minimum, max_gap, &unique[strand])
            .into_iter()
            .enumerate()
        {
            let lo = runs[0].q;
            let last = runs.last().unwrap();
            let hi = last.q + last.n;
            let mask: Vec<_> = (lo..hi)
                .map(|q| available[if strand == 1 { length - 1 - q } else { q }])
                .collect();
            for (a, b) in intervals(&mask) {
                let selected = clip(&runs, a + lo, b + lo, minimum, &unique[strand]);
                if !selected.is_empty() {
                    supplements.push(Block {
                        identity: format!("G{j:06}{}", if strand == 0 { '+' } else { '-' }),
                        kind: "graph".into(),
                        native_kind: None,
                        reverse: strand == 1,
                        runs: selected,
                        record: None,
                    });
                }
            }
        }
    }
    let mut count = vec![0u32; length];
    for b in &supplements {
        let (s, e) = native_bounds(b, length);
        for v in &mut count[s..e] {
            *v += 1;
        }
    }
    for b in supplements {
        let (s, e) = native_bounds(&b, length);
        if count[s..e].iter().all(|&v| v == 1) {
            accepted.push(b);
        }
    }
    accepted = reconcile_inversions(
        graph, &oriented, &unique, accepted, &owners, minimum, mapping_k,
    )?;
    occupied.fill(false);
    for b in &accepted {
        let (s, e) = native_bounds(b, length);
        occupied[s..e].fill(true);
    }
    let missing = intervals(&occupied.iter().map(|v| !v).collect::<Vec<_>>());
    drop(unique);
    Ok(Plan {
        blocks: accepted,
        missing,
        oriented,
    })
}

fn bounded_inversions(
    graph: &Graph,
    oriented: &[Oriented; 2],
    unique: &[Uniqueness; 2],
    anchor: usize,
) -> Vec<Block> {
    let mut candidates = Vec::new();
    for strand in 0..2 {
        let view = &oriented[strand];
        let mut counts = HashMap::new();
        for n in &view.walk {
            *counts.entry(n.abs()).or_insert(0usize) += 1;
        }
        let location = |node: NodeId| {
            graph
                .ref_places
                .get(&node.abs())
                .and_then(|p| (p.len() == 1 && counts.get(&node.abs()) == Some(&1)).then_some(p[0]))
        };
        let mut i = 0;
        while i < view.walk.len() {
            let high = location(view.walk[i]);
            if high.is_none_or(|h| graph.refwalk()[h] != -view.walk[i]) {
                i += 1;
                continue;
            }
            let high = high.unwrap();
            let mut j = i + 1;
            while j < view.walk.len()
                && j - i <= high
                && location(view.walk[j]) == Some(high - (j - i))
                && graph.refwalk()[high - (j - i)] == -view.walk[j]
            {
                j += 1;
            }
            let low = high + 1 - (j - i);
            let (mut li, mut lr) = (i, low);
            let (mut ri, mut rr) = (j, high + 1);
            while li > 0 && lr > 0 && view.walk[li - 1] == graph.refwalk()[lr - 1] {
                li -= 1;
                lr -= 1;
            }
            while ri < view.walk.len()
                && rr < graph.refwalk().len()
                && view.walk[ri] == graph.refwalk()[rr]
            {
                ri += 1;
                rr += 1;
            }
            let (qs, qe) = (view.offsets[i], view.offsets[j]);
            let (rs, re) = (graph.ref_offsets[low], graph.ref_offsets[high + 1]);
            if qs - view.offsets[li] >= anchor
                && view.offsets[ri] - qe >= anchor
                && view.query[qs..qe] != graph.reference[rs..re]
                && view.query[qs..qe] == io::reverse_complement(&graph.reference[rs..re])
                && unique[strand].oriented(&view.query[qs - anchor..qs])
                && unique[strand].oriented(&view.query[qe..qe + anchor])
            {
                let sign = if strand == 0 { '+' } else { '-' };
                let id = io::hash_bytes(format!("{sign}:{qs}:{qe}:{rs}:{re}").as_bytes());
                let anchors = [
                    Anchor {
                        r: rs - anchor,
                        q: qs - anchor,
                        n: anchor,
                        node: None,
                    },
                    Anchor {
                        r: re,
                        q: qe,
                        n: anchor,
                        node: None,
                    },
                ];
                candidates.push(Block {
                    identity: format!("I{}", &id[..16]),
                    kind: "inversion".into(),
                    native_kind: Some("INV".into()),
                    reverse: strand == 1,
                    record: None,
                    runs: anchors
                        .into_iter()
                        .map(|a| Run {
                            r: a.r,
                            q: a.q,
                            n: a.n,
                            payload: vec![a],
                        })
                        .collect(),
                });
            }
            i = j;
        }
    }
    candidates
}
fn reconcile_inversions(
    graph: &Graph,
    oriented: &[Oriented; 2],
    unique: &[Uniqueness; 2],
    blocks: Vec<Block>,
    owners: &[u32],
    minimum: usize,
    mapping_k: usize,
) -> Result<Vec<Block>> {
    let length = oriented[0].query.len();
    let mut candidates = bounded_inversions(graph, oriented, unique, minimum);
    // PGGB can spell an inverted interval using new forward-oriented segments.
    // Recover only exact same-interval reversals between two independent unique
    // forward anchors. This is not permission to stitch arbitrary PAF gaps.
    let graph_bounds: Vec<_> = candidates
        .iter()
        .map(|b| native_bounds(b, length))
        .collect();
    for strand in 0..2 {
        let view = &oriented[strand];
        let mut runs: Vec<_> = blocks
            .iter()
            .filter(|b| b.reverse == (strand == 1))
            .flat_map(|b| b.runs.iter())
            .filter(|r| r.n >= minimum)
            .collect();
        runs.sort_by_key(|r| (r.q, r.r, r.n));
        runs.dedup_by_key(|r| (r.q, r.r, r.n));
        for pair in runs.windows(2) {
            let (left, right) = (pair[0], pair[1]);
            let (rs, re, qs, qe) = (left.r + left.n, right.r, left.q + left.n, right.q);
            if re <= rs || qe <= qs || re - rs != qe - qs || re - rs < mapping_k {
                continue;
            }
            let (native_s, native_e) = if strand == 0 {
                (qs - minimum, qe + minimum)
            } else {
                (length - qe - minimum, length - qs + minimum)
            };
            if graph_bounds
                .iter()
                .any(|&(s, e)| native_s < e && s < native_e)
                || graph.reference[rs..re] == view.query[qs..qe]
                || view.query[qs..qe] != io::reverse_complement(&graph.reference[rs..re])
                || !unique[strand].oriented(&view.query[qs - minimum..qs])
                || !unique[strand].oriented(&view.query[qe..qe + minimum])
            {
                continue;
            }
            let k = mapping_k.min(31);
            if !(rs..=re - k).any(|p| unique[strand].inverted(&graph.reference[p..p + k])) {
                continue;
            }
            let sign = if strand == 0 { '+' } else { '-' };
            let id = io::hash_bytes(format!("sequence:{sign}:{qs}:{qe}:{rs}:{re}").as_bytes());
            let runs = [(rs - minimum, qs - minimum), (re, qe)]
                .into_iter()
                .map(|(r, q)| Run {
                    r,
                    q,
                    n: minimum,
                    payload: vec![Anchor {
                        r,
                        q,
                        n: minimum,
                        node: None,
                    }],
                })
                .collect();
            candidates.push(Block {
                identity: format!("I{}", &id[..16]),
                kind: "inversion".into(),
                native_kind: Some("INV".into()),
                reverse: strand == 1,
                record: None,
                runs,
            });
        }
    }
    let bounds: Vec<_> = candidates
        .iter()
        .map(|b| native_bounds(b, length))
        .collect();
    let accepted: Vec<_> = candidates
        .into_iter()
        .enumerate()
        .filter_map(|(i, b)| {
            let (s, e) = bounds[i];
            (!owners[s..e].iter().any(|&v| v > 1)
                && !bounds
                    .iter()
                    .enumerate()
                    .any(|(j, &(a, b))| i != j && s < b && a < e))
            .then_some(b)
        })
        .collect();
    if accepted.is_empty() {
        return Ok(blocks);
    }
    let mut occupied: Vec<_> = accepted.iter().map(|b| native_bounds(b, length)).collect();
    occupied.sort_unstable();
    let mut result = Vec::new();
    for block in blocks {
        let (start, end) = native_bounds(&block, length);
        let mut pieces = vec![(start, end)];
        for &(a, b) in &occupied {
            pieces = pieces
                .into_iter()
                .flat_map(|(lo, hi)| {
                    if lo < b && a < hi {
                        vec![(lo, hi.min(a)), (lo.max(b), hi)]
                    } else {
                        vec![(lo, hi)]
                    }
                })
                .filter(|&(s, e)| s < e)
                .collect();
        }
        if pieces == vec![(start, end)] {
            result.push(block);
            continue;
        }
        for (s, e) in pieces {
            let (lo, hi) = if block.reverse {
                (length - e, length - s)
            } else {
                (s, e)
            };
            let runs = clip(
                &block.runs,
                lo,
                hi,
                mapping_k,
                &unique[block.reverse as usize],
            );
            if !runs.is_empty() {
                result.push(Block {
                    runs,
                    ..block.clone()
                });
            }
        }
    }
    result.extend(accepted);
    Ok(result)
}
pub type StructureRow = (usize, usize, usize, usize, String, String);
pub fn structure_evidence(
    blocks: &[Block],
    length: usize,
    ref_length: usize,
) -> Result<(Vec<StructureRow>, Intervals)> {
    let mut ordered: Vec<_> = blocks.iter().collect();
    ordered.sort_by_key(|b| native_bounds(b, length));
    let mut rows = Vec::new();
    let mut masks = Vec::new();
    for b in &ordered {
        if b.kind == "inversion" {
            let left = &b.runs[0];
            let right = b.runs.last().unwrap();
            let (mut qs, mut qe) = (left.q + left.n, right.q);
            if b.reverse {
                (qs, qe) = (length - qe, length - qs);
            }
            rows.push((
                qs,
                qe,
                left.r + left.n,
                right.r,
                "RESOLVED_INVERSION".into(),
                b.identity.clone(),
            ));
        }
    }
    for pair in ordered.windows(2) {
        let (left, right) = (pair[0], pair[1]);
        let (ls, le) = native_bounds(left, length);
        let (rs, re) = native_bounds(right, length);
        let lr = (
            left.runs[0].r,
            left.runs.last().unwrap().r + left.runs.last().unwrap().n,
        );
        let rr = (
            right.runs[0].r,
            right.runs.last().unwrap().r + right.runs.last().unwrap().n,
        );
        let reason = if left.reverse != right.reverse {
            Some("UNRESOLVED_ORIENTATION")
        } else if (!left.reverse && lr.1 > rr.0) || (left.reverse && rr.1 > lr.0) {
            Some("UNRESOLVED_ORDER")
        } else {
            None
        };
        let ids = format!("{},{}", left.identity, right.identity);
        if let Some(reason) = reason {
            masks.extend([lr, rr]);
            rows.push((ls, re, lr.0.min(rr.0), lr.1.max(rr.1), reason.into(), ids));
        } else if left.identity != right.identity {
            let gap = if left.reverse {
                (rr.1, lr.0)
            } else {
                (lr.1, rr.0)
            };
            if rs > le || gap.1 > gap.0 {
                rows.push((le, rs, gap.0, gap.1, "UNRESOLVED_ADJACENCY".into(), ids));
                masks.push(gap);
            }
        }
    }
    Ok((rows, Ranges::new(masks, ref_length)?.rows))
}
