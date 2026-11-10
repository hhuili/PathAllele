//! One locus per boundary crossing; all sample-level work stays in Rust.
use crate::model::{dna, intersects, replay, Edit, EditInput, Evidence, Intervals, Result};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

pub const STATUS: [&str; 11] = [
    "PASS",
    "NO_DATA",
    "MULTIPLE_OCCURRENCES",
    "AMBIGUOUS_ALIGNMENT",
    "AMBIGUOUS_BASE",
    "SPANNING_EDIT",
    "OVERLAPPING_EDIT",
    "UNCERTAIN_BREAKPOINT",
    "PARTIAL_SEQUENCE",
    "UNRESOLVED_COMPONENT",
    "UNRESOLVED_STRUCTURE",
];

pub struct ProjectionRow {
    pub alleles: Vec<String>,
    pub alts: Vec<String>,
    pub cells: String,
    pub links: String,
    pub status_counts: Vec<usize>,
    pub gt: Vec<i32>,
    pub projected: Vec<i32>,
    pub source: Vec<u8>,
    pub export: Vec<u8>,
}

pub struct Projection {
    evidence: Evidence,
    names: Vec<String>,
    positive: Vec<HashSet<usize>>,
    components: HashMap<(usize, usize, String), String>,
}

impl Projection {
    fn observe(
        &self,
        si: usize,
        s: usize,
        e: usize,
        snp: bool,
    ) -> (Option<Cow<'_, str>>, u8, Vec<usize>) {
        let evidence = &self.evidence;
        let row = &evidence.samples[si];
        if snp {
            if row.structural.contains(s)
                && !row.chain_covers(s.saturating_sub(1), (s + 2).min(evidence.reference.len()))
            {
                return (None, 10, vec![]);
            }
            if row.multiple.contains(s) {
                return (None, 2, vec![]);
            }
            if row.uncertain.contains(s) {
                return (None, 3, vec![]);
            }
            if !row.unique.contains(s) {
                return (None, 1, vec![]);
            }
            let base = row.base(&evidence.reference, &evidence.pool, s);
            return match base {
                b'A' => (Some(Cow::Borrowed("A")), 0, vec![]),
                b'C' => (Some(Cow::Borrowed("C")), 0, vec![]),
                b'G' => (Some(Cow::Borrowed("G")), 0, vec![]),
                b'T' => (Some(Cow::Borrowed("T")), 0, vec![]),
                b'-' => (None, 5, vec![]),
                b'?' => (None, 1, vec![]),
                _ => (None, 4, vec![]),
            };
        }
        let lo = s.saturating_sub(1);
        let hi = (e + 1).min(evidence.reference.len());
        let (inside_lo, inside_hi) = if s < e { (s, e) } else { (lo, hi) };
        if row.structural.overlaps(inside_lo, inside_hi) && !row.chain_covers(lo, hi) {
            return (None, 10, vec![]);
        }
        if row.multiple.overlaps(inside_lo, inside_hi) {
            return (None, 2, vec![]);
        }
        if row.uncertain.overlaps(inside_lo, inside_hi) {
            return (None, 3, vec![]);
        }
        if !row.chain_covers(lo, hi) || !row.unique.covers(inside_lo, inside_hi) {
            return (None, 1, vec![]);
        }
        let all = row.overlapping(&evidence.pool, lo, hi);
        let local: Vec<_> = all
            .iter()
            .copied()
            .filter(|&i| {
                let edit = &evidence.pool[i];
                intersects(s, e, edit) && (s != e || edit.structural())
            })
            .collect();
        if local
            .iter()
            .any(|&i| evidence.pool[i].start < s || evidence.pool[i].end > e)
        {
            return (None, 6, local);
        }
        for &id in &all {
            let x = &evidence.pool[id];
            let flank = (s > 0 && x.start < s && s - 1 < x.end)
                || (e < evidence.reference.len() && x.start <= e && e < x.end);
            if flank && !dna(&x.alt) {
                return (None, 7, local);
            }
        }
        if local.is_empty() {
            return (Some(Cow::Borrowed(&evidence.reference[s..e])), 0, local);
        }
        let edits: Vec<_> = local.iter().map(|&i| &evidence.pool[i]).collect();
        match replay(&evidence.reference, s, e, &edits) {
            Ok(value) if dna(&value) => (Some(Cow::Owned(value)), 0, local),
            Ok(_) => (None, 8, local),
            Err(_) => (None, 6, local),
        }
    }

    pub fn project(&mut self, ident: &str, s: usize, e: usize, snp: bool) -> Result<ProjectionRow> {
        if s > e || e > self.evidence.reference.len() || (snp && e != s + 1) {
            return Err("Invalid locus bounds/kind".into());
        }
        if ident.is_empty() || ident.bytes().any(|b| b.is_ascii_whitespace()) {
            return Err("Invalid locus identifier".into());
        }
        let reference = &self.evidence.reference[s..e];
        let n = self.names.len();
        let mut sequences = vec![reference.to_string()];
        let mut index: HashMap<String, i32> = HashMap::from([(reference.to_string(), 0)]);
        let mut gt = Vec::with_capacity(n);
        let mut projected = Vec::with_capacity(n);
        let mut source = Vec::with_capacity(n);
        let mut export = Vec::with_capacity(n);
        let mut factored = Vec::new();
        for si in 0..n {
            let (observed, status, local) = self.observe(si, s, e, snp);
            let Some(observed) = observed else {
                gt.push(-1);
                projected.push(-1);
                source.push(status);
                export.push(status);
                continue;
            };
            let ai = if let Some(&i) = index.get(observed.as_ref()) {
                i
            } else {
                let i = i32::try_from(sequences.len()).map_err(|_| "Too many alleles")?;
                let value = observed.to_string();
                index.insert(value.clone(), i);
                sequences.push(value);
                i
            };
            gt.push(ai);
            source.push(0);
            if snp || ai == 0 {
                projected.push(ai);
                export.push(0);
                continue;
            }
            if !local.is_empty() && local.iter().all(|id| self.positive[si].contains(id)) {
                // observe() already reconstructed this allele from these edits.
                projected.push(0);
                export.push(0);
                let mut components = Vec::new();
                for &id in &local {
                    let x = &self.evidence.pool[id];
                    let key = (x.start, x.end, x.alt.clone());
                    let component = self
                        .components
                        .get(&key)
                        .ok_or("Missing emitted component identifier")?;
                    components.push(component.clone());
                }
                components.sort_unstable();
                components.dedup();
                factored.push((si, ai, components));
            } else if local.len() == 1 {
                let x = &self.evidence.pool[local[0]];
                if (x.start, x.end, x.alt.as_str()) == (s, e, observed.as_ref()) {
                    projected.push(ai);
                    export.push(0);
                } else {
                    projected.push(-1);
                    export.push(9);
                }
            } else {
                projected.push(-1);
                export.push(9);
            }
        }
        let mut order: Vec<usize> = (1..sequences.len()).collect();
        order.sort_unstable_by(|&a, &b| {
            (sequences[a].len(), &sequences[a]).cmp(&(sequences[b].len(), &sequences[b]))
        });
        order.insert(0, 0);
        let mut remap = vec![0i32; order.len()];
        let alleles: Vec<_> = order
            .iter()
            .enumerate()
            .map(|(new, &old)| {
                remap[old] = new as i32;
                sequences[old].clone()
            })
            .collect();
        for values in [&mut gt, &mut projected] {
            for value in values {
                if *value >= 0 {
                    *value = remap[*value as usize];
                }
            }
        }
        let mut used = vec![false; alleles.len()];
        for &value in &projected {
            if value > 0 {
                used[value as usize] = true;
            }
        }
        let alt_indices: Vec<_> = (1..alleles.len()).filter(|&i| used[i]).collect();
        let alts: Vec<_> = alt_indices.iter().map(|&i| alleles[i].clone()).collect();
        let mut vcf_index = vec![0; alleles.len()];
        let mut model_ids = vec![Vec::new(); alleles.len()];
        for (i, &ai) in alt_indices.iter().enumerate() {
            vcf_index[ai] = i + 1;
            self.components
                .insert((s, e, alleles[ai].clone()), ident.to_string());
            for complex in [false, true] {
                if let Some(&id) = self.evidence.lookup.get(&Edit {
                    start: s,
                    end: e,
                    alt: alleles[ai].clone(),
                    complex,
                }) {
                    model_ids[ai].push(id);
                }
            }
        }
        let mut cells = String::new();
        let mut status_counts = vec![0; STATUS.len()];
        if !alts.is_empty() {
            for si in 0..n {
                if si > 0 {
                    cells.push('\t');
                }
                let value = projected[si];
                if value < 0 {
                    cells.push('.');
                } else {
                    write!(&mut cells, "{}", vcf_index[value as usize]).expect("write to String");
                }
                cells.push(':');
                cells.push_str(STATUS[export[si] as usize]);
                status_counts[export[si] as usize] += 1;
                if value > 0 {
                    for &id in &model_ids[value as usize] {
                        if self.evidence.samples[si].members.binary_search(&id).is_ok() {
                            self.positive[si].insert(id);
                        }
                    }
                }
            }
        }
        let mut links = String::new();
        for (si, ai, components) in factored {
            writeln!(
                &mut links,
                "{}\t{}\t{}\t{}",
                ident,
                self.names[si],
                remap[ai as usize],
                components.join(",")
            )
            .expect("write to String");
        }
        Ok(ProjectionRow {
            alleles,
            alts,
            cells,
            links,
            status_counts,
            gt,
            projected,
            source,
            export,
        })
    }
}

impl Projection {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: String,
        names: Vec<String>,
        models: Vec<Vec<EditInput>>,
        unique: Vec<Intervals>,
        multiple: Vec<Intervals>,
        uncertain: Vec<Intervals>,
        structural: Vec<Intervals>,
        chains: Vec<Intervals>,
    ) -> Result<Self> {
        if names.len() != models.len()
            || names
                .iter()
                .any(|n| n.is_empty() || n.bytes().any(|b| b.is_ascii_whitespace()))
        {
            return Err("Invalid sample names or dimensions".into());
        }
        let evidence = Evidence::new(
            reference, models, unique, multiple, uncertain, structural, chains,
        )?;
        let positive = (0..names.len()).map(|_| HashSet::new()).collect();
        Ok(Self {
            evidence,
            names,
            positive,
            components: HashMap::new(),
        })
    }
}
