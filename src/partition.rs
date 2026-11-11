//! Sequence units are bounded by matches compatible with structural evidence.
//!
//! A forward graph match may locate a source without being a valid cut through
//! an inverted sequence relation. Establish those relations on complete source
//! alleles before choosing cuts. Both direct and partial projection consume the
//! resulting units; neither may split a protected relation during alignment.
use crate::{
    anchors::{self, Anchor},
    classify::{Classifier, Segment},
    io,
    model::{dna, Result},
    source::Parent,
};

#[derive(Clone, Debug)]
pub struct Unit {
    pub r: usize,
    pub q: usize,
    pub rn: usize,
    pub qn: usize,
    pub preserve: bool,
}

#[derive(Clone)]
pub struct Plan {
    pub parents: Vec<Parent>,
    pub units: Vec<Vec<Vec<Unit>>>,
}

impl Plan {
    pub fn reversed(&self) -> Self {
        let mut units = self.units.clone();
        for (p, aa) in self.parents.iter().zip(&mut units) {
            for (ai, uu) in aa.iter_mut().enumerate() {
                for u in uu.iter_mut() {
                    u.r = p.alleles[0].len() - u.r - u.rn;
                    u.q = p.alleles[ai].len() - u.q - u.qn;
                }
                uu.reverse();
            }
        }
        Self {
            parents: reverse_parents(&self.parents),
            units,
        }
    }

    pub fn select(&self, selected: &[(usize, usize)]) -> Self {
        let mut parents = Vec::with_capacity(selected.len());
        let mut units = Vec::with_capacity(selected.len());
        for &(pi, ai) in selected {
            let p = &self.parents[pi];
            let mut v = p.clone();
            v.alleles = vec![p.alleles[0].clone(), p.alleles[ai].clone()];
            v.anchors = vec![p.anchors[0].clone(), p.anchors[ai].clone()];
            v.kinds = vec!["REFERENCE".into(), p.kinds[ai].clone()];
            v.carriers = vec![p.carriers[0].clone(), p.carriers[ai].clone()];
            parents.push(v);
            units.push(vec![vec![], self.units[pi][ai].clone()]);
        }
        Self { parents, units }
    }
}

fn reverse_parents(parents: &[Parent]) -> Vec<Parent> {
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
// Retain only portions of a forward match that lie on the same side of
// every relation on both axes. A match crossing a boundary may still provide
// useful flanking sequence; rejecting it wholesale would join unrelated edits.
fn compatible_parts(a: Anchor, start: usize, relations: &[Segment]) -> Vec<Anchor> {
    let mut spans = vec![(0, a.n)];
    for s in relations {
        let before = s
            .start
            .saturating_sub(start + a.r)
            .min(s.query.saturating_sub(a.q))
            .min(a.n);
        let after = s
            .end
            .saturating_sub(start + a.r)
            .max(s.query_end.saturating_sub(a.q))
            .min(a.n);
        spans = spans
            .into_iter()
            .flat_map(|(lo, hi)| {
                [(lo, hi.min(before)), (lo.max(after), hi)]
                    .into_iter()
                    .filter(|(l, h)| l < h)
            })
            .collect();
    }
    spans
        .into_iter()
        .map(|(lo, hi)| Anchor {
            r: a.r + lo,
            q: a.q + lo,
            n: hi - lo,
            node: a.node,
        })
        .collect()
}

fn units(p: &Parent, ai: usize, relations: &[Segment]) -> Result<Vec<Unit>> {
    let (reference, alternate) = (&p.alleles[0], &p.alleles[ai]);
    if matches!(p.kinds[ai].as_str(), "INV" | "COMPLEX") {
        return Ok(vec![Unit {
            r: 1,
            q: 1,
            rn: reference.len() - 2,
            qn: alternate.len() - 2,
            preserve: true,
        }]);
    }
    // Validate the same padding, exact-match and monotonicity contract used by
    // source alleles, then retain both coordinates for forward/reverse rescue.
    anchors::gaps(reference, alternate, &p.anchors[ai])?;
    let mut cuts = vec![(0, 0, 1)];
    cuts.extend(p.anchors[ai].iter().map(|a| (a.r, a.q, a.n)));
    cuts.push((reference.len() - 1, alternate.len() - 1, 1));
    let (mut re, mut qe) = (0, 0);
    let mut out = vec![];
    for (r, q, n) in cuts {
        if reference[re..r] != alternate[qe..q] {
            let preserve = relations.iter().any(|s| {
                p.start + re < s.end && s.start < p.start + r && qe < s.query_end && s.query < q
            });
            out.push(Unit {
                r: re,
                q: qe,
                rn: r - re,
                qn: q - qe,
                preserve,
            });
        }
        re = r + n;
        qe = q + n;
    }
    Ok(out)
}

/// Complete and partial plans share structural constraints. The complete plan
/// additionally uses exact-sequence anchors; partial rescue keeps its existing
/// graph-context boundaries and aligns only unprotected units.
pub fn build(
    sources: &[Parent],
    reference: &str,
    sequence_anchor: usize,
    minimum: usize,
) -> Result<(Plan, Plan)> {
    let classifier = Classifier::new(reference, minimum)?;
    let mut complete = Plan {
        parents: sources.to_vec(),
        units: vec![],
    };
    let mut partial = complete.clone();
    for (p, q) in complete.parents.iter_mut().zip(&mut partial.parents) {
        let mut pu = vec![vec![]; p.alleles.len()];
        let mut qu = pu.clone();
        for ai in 0..p.alleles.len() {
            let primitive: Vec<_> = anchors::context_runs(
                &p.alleles[0],
                &p.alleles[ai],
                p.anchors[ai].clone(),
                1,
                None,
            )?
            .into_iter()
            .flat_map(|r| r.payload)
            .collect();
            if ai == 0 {
                q.anchors[ai] = primitive;
                continue;
            }
            let relations = if p.kinds[ai] == "SEQUENCE"
                && p.alleles[0] != p.alleles[ai]
                && dna(&p.alleles[0])
                && dna(&p.alleles[ai])
            {
                if reference.get(p.start..p.start + p.alleles[0].len()) != Some(&p.alleles[0]) {
                    return Err("Source allele REF disagrees with reference".into());
                }
                classifier
                    .annotate(
                        p.start,
                        p.start + p.alleles[0].len(),
                        &p.alleles[ai],
                        &p.kinds[ai],
                    )?
                    .segments
            } else {
                vec![]
            };
            q.anchors[ai] = primitive
                .into_iter()
                .flat_map(|a| compatible_parts(a, p.start, &relations))
                .collect();
            if p.kinds[ai] == "SEQUENCE" {
                p.anchors[ai] = anchors::refine(
                    &p.alleles[0],
                    &p.alleles[ai],
                    &q.anchors[ai],
                    sequence_anchor,
                )?
                .into_iter()
                .flat_map(|a| compatible_parts(a, p.start, &relations))
                .collect();
            }
            pu[ai] = units(p, ai, &relations)?;
            qu[ai] = units(q, ai, &relations)?;
        }
        complete.units.push(pu);
        partial.units.push(qu);
    }
    Ok((complete, partial))
}
