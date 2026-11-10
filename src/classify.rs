//! One sequence/structure annotation contract for every complete catalogue ALT.
//!
//! Reverse homology is evidence of an inverted sequence component, not a unique
//! mutational history. Only exact reversal of the same interval is a pure INV.
use crate::{
    catalogue::Catalogue,
    io,
    model::{self, dna, Result},
};
use std::{collections::HashMap, path::Path};

pub const FILE: &str = "Graph.variant_annotations.tsv.gz";
pub const COLUMNS: [&str; 8] = [
    "locus_id",
    "allele_index",
    "type",
    "sequence_type",
    "size_class",
    "native",
    "inversion_basis",
    "inversion_segments",
];

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Segment {
    pub start: usize,
    pub end: usize,
    pub query: usize,
    pub query_end: usize,
    pub template: String,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Annotation {
    pub kind: String,
    pub sequence: String,
    pub size: String,
    pub native: String,
    pub basis: String,
    pub segments: Vec<Segment>,
}

pub fn shape(reference: &str, alternate: &str) -> Result<(&'static str, &'static str)> {
    if (!reference.is_empty() && !dna(reference)) || (!alternate.is_empty() && !dna(alternate)) {
        return Err("Classification requires explicit ACGT alleles".into());
    }
    let edit = model::trim(0, reference, alternate);
    let (r, a) = (edit.end - edit.start, edit.alt.len());
    let kind = match (r, a) {
        (0, 0) => return Err("ALT is identical to REF".into()),
        (0, _) => "INS",
        (_, 0) => "DEL",
        (1, 1) => "SNV",
        _ if r == a => "MNV",
        _ => "REPLACEMENT",
    };
    Ok((kind, if r.max(a) >= 50 { "SV" } else { "SMALL" }))
}

// Exact two-bit canonical keys: no hash collisions can establish homology.
fn seeds(sequence: &str, k: usize) -> Vec<(u64, usize, bool)> {
    let (mut f, mut r, mut valid) = (0u64, 0u64, 0usize);
    let mask = (1u64 << (2 * k)) - 1;
    let mut result = Vec::with_capacity(sequence.len().saturating_sub(k - 1));
    for (i, b) in sequence.bytes().enumerate() {
        let bits = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                f = 0;
                r = 0;
                valid = 0;
                continue;
            }
        };
        f = ((f << 2) | bits) & mask;
        r = (r >> 2) | ((3 - bits) << (2 * (k - 1)));
        valid += 1;
        if valid >= k && f != r {
            result.push((f.min(r), i + 1 - k, f < r));
        }
    }
    result
}

pub struct Classifier<'a> {
    reference: &'a str,
    minimum: usize,
    k: usize,
    unique: HashMap<u64, Option<(usize, bool)>>,
}
impl<'a> Classifier<'a> {
    pub fn new(reference: &'a str, minimum: usize) -> Result<Self> {
        if minimum < 19 {
            return Err("Inversion evidence requires a minimum span of at least 19 bp".into());
        }
        let k = minimum.min(31);
        let mut unique = HashMap::new();
        for (key, pos, orientation) in seeds(reference, k) {
            unique
                .entry(key)
                .and_modify(|v| *v = None)
                .or_insert(Some((pos, orientation)));
        }
        Ok(Self {
            reference,
            minimum,
            k,
            unique,
        })
    }

    pub fn annotate(
        &self,
        start: usize,
        end: usize,
        alternate: &str,
        native: &str,
    ) -> Result<Annotation> {
        let reference = self
            .reference
            .get(start..end)
            .ok_or("Invalid annotation bounds")?;
        let (sequence, size) = shape(reference, alternate)?;
        check_native(native)?;
        let edit = model::trim(start, reference, alternate);
        let mut segments = vec![];
        if edit.end - edit.start >= self.minimum && edit.alt.len() >= self.minimum {
            let words = seeds(alternate, self.k);
            let lo = edit.start - start;
            let hi = lo + edit.alt.len();
            let mut candidates = vec![];
            for &(key, q, orientation) in &words {
                if q < lo || q + self.k > hi {
                    continue;
                }
                if let Some(Some((r, forward))) = self.unique.get(&key) {
                    // A distant reverse-copy insertion is not a local inversion.
                    if *forward != orientation && *r >= edit.start && r + self.k <= edit.end {
                        candidates.push((key, r + q, q, *r));
                    }
                }
            }
            // Most alleles have no local reverse seed. Only count candidate
            // keys, but across the complete allele including its padding.
            let mut counts: HashMap<_, u32> = candidates.iter().map(|x| (x.0, 0)).collect();
            if !counts.is_empty() {
                for (key, _, _) in words {
                    if let Some(n) = counts.get_mut(&key) {
                        *n += 1;
                    }
                }
            }
            let mut hits: Vec<_> = candidates
                .into_iter()
                .filter_map(|(key, d, q, r)| (counts[&key] == 1).then_some((d, q, r)))
                .collect();
            hits.sort_unstable();
            let mut i = 0;
            while i < hits.len() {
                let (diagonal, q, r) = hits[i];
                let mut j = i + 1;
                let mut qe = q + self.k;
                let mut rs = r;
                while j < hits.len() && hits[j].0 == diagonal && hits[j].1 <= qe {
                    qe = qe.max(hits[j].1 + self.k);
                    rs = rs.min(hits[j].2);
                    j += 1;
                }
                let (mut qs, mut re) = (q, r + self.k);
                let (rb, ab) = (self.reference.as_bytes(), alternate.as_bytes());
                while qs > lo && re < rb.len() && complement(rb[re]) == ab[qs - 1] {
                    qs -= 1;
                    re += 1;
                }
                while qe < hi && rs > 0 && complement(rb[rs - 1]) == ab[qe] {
                    rs -= 1;
                    qe += 1;
                }
                // Self reverse-complement-matching flanks do not establish
                // orientation. Require an informative core and a unique seed
                // wholly inside it, so palindrome padding cannot promote a
                // short substitution in a larger compound locus to inversion.
                let template = &self.reference[rs..re];
                let core = model::trim(rs, template, &io::reverse_complement(template));
                let informative = core.end - core.start >= self.minimum
                    && hits[i..j]
                        .iter()
                        .any(|&(_, _, r)| r >= core.start && r + self.k <= core.end);
                if informative && rs < edit.end && edit.start < re {
                    segments.push(Segment {
                        start: rs,
                        end: re,
                        query: qs,
                        query_end: qe,
                        template: self.reference[rs..re].into(),
                    });
                }
                i = j;
            }
            segments.sort_by_key(|s| (s.query, s.start, s.end));
            segments.dedup();
        }
        // A graph/path-bounded inversion may be shorter than the sequence-only
        // threshold. It still must replay as the exact same-interval RC.
        let core = &self.reference[edit.start..edit.end];
        if segments.is_empty()
            && native.split('+').any(|k| k == "INV")
            && !core.is_empty()
            && edit.alt == io::reverse_complement(core)
        {
            segments.push(Segment {
                start: edit.start,
                end: edit.end,
                query: edit.start - start,
                query_end: edit.start - start + edit.alt.len(),
                template: core.into(),
            });
        }
        let mut result = Annotation {
            kind: sequence.into(),
            sequence: sequence.into(),
            size: size.into(),
            native: native.into(),
            basis: "NONE".into(),
            segments,
        };
        result.assign(start, reference, alternate);
        result.verify(start, reference, alternate, self.minimum)?;
        Ok(result)
    }
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => 0,
    }
}
fn check_native(value: &str) -> Result<()> {
    let keys: Vec<_> = value.split('+').collect();
    if keys.is_empty()
        || keys
            .iter()
            .any(|k| !matches!(*k, "SEQUENCE" | "INV" | "COMPLEX"))
        || keys.windows(2).any(|w| w[0] >= w[1])
    {
        return Err("Invalid native orientation class".into());
    }
    Ok(())
}
impl Annotation {
    fn assign(&mut self, start: usize, reference: &str, alternate: &str) {
        let e = model::trim(start, reference, alternate);
        let exact = !self.segments.is_empty()
            && self.segments.iter().any(|s| {
                s.start == e.start
                    && s.end == e.end
                    && s.query == e.start - start
                    && s.query_end == e.start - start + e.alt.len()
            });
        if exact {
            self.kind = "INV".into();
            self.basis = if self.native.split('+').any(|k| k == "INV") {
                "PATH_AND_RC"
            } else {
                "EXACT_RC"
            }
            .into();
        } else if !self.segments.is_empty() {
            self.kind = "INV_COMPLEX".into();
            self.basis = "LOCAL_RC".into();
        } else if self
            .native
            .split('+')
            .any(|k| matches!(k, "INV" | "COMPLEX"))
        {
            self.kind = "COMPLEX".into();
            self.basis = "NONE".into();
        }
    }
    pub fn evidence(&self) -> String {
        if self.segments.is_empty() {
            return ".".into();
        }
        self.segments
            .iter()
            .map(|s| {
                format!(
                    "{}:{}:{}:{}:{}",
                    s.start, s.end, s.query, s.query_end, s.template
                )
            })
            .collect::<Vec<_>>()
            .join("|")
    }
    pub fn fields(&self, id: &str, allele: usize) -> Vec<String> {
        vec![
            id.into(),
            allele.to_string(),
            self.kind.clone(),
            self.sequence.clone(),
            self.size.clone(),
            self.native.clone(),
            self.basis.clone(),
            self.evidence(),
        ]
    }
    fn verify(&self, start: usize, reference: &str, alternate: &str, minimum: usize) -> Result<()> {
        check_native(&self.native)?;
        let (sequence, size) = shape(reference, alternate)?;
        if self.sequence != sequence || self.size != size {
            return Err("Allele sequence classification mismatch".into());
        }
        let e = model::trim(start, reference, alternate);
        let mut seen = std::collections::HashSet::new();
        for s in &self.segments {
            let reversed = io::reverse_complement(&s.template);
            let core = model::trim(0, &s.template, &reversed);
            if e.start == e.end
                || e.alt.is_empty()
                || s.start >= s.end
                || s.query >= s.query_end
                || s.end - s.start != s.query_end - s.query
                || s.end - s.start != s.template.len()
                || !dna(&s.template)
                || s.query < e.start - start
                || s.query_end > e.start - start + e.alt.len()
                || s.start >= e.end
                || s.end <= e.start
                || (self.basis != "PATH_AND_RC" && core.end - core.start < minimum)
                || alternate.get(s.query..s.query_end) != Some(reversed.as_str())
                || !seen.insert((s.start, s.end, s.query, s.query_end))
            {
                return Err("Invalid inversion sequence witness".into());
            }
            // Check every base of the template that is locally available even
            // in portable validation without the original reference FASTA.
            let lo = s.start.max(start);
            let hi = s.end.min(start + reference.len());
            if reference[lo - start..hi - start] != s.template[lo - s.start..hi - s.start] {
                return Err("Inversion template disagrees with catalogue REF".into());
            }
        }
        let mut expected = self.clone();
        expected.kind = sequence.into();
        expected.basis = "NONE".into();
        expected.assign(start, reference, alternate);
        if self.kind != expected.kind || self.basis != expected.basis {
            return Err("Inversion classification disagrees with sequence witnesses".into());
        }
        Ok(())
    }
}
pub fn info(rows: &[&Annotation]) -> String {
    [
        (
            "TYPE",
            rows.iter().map(|a| a.kind.as_str()).collect::<Vec<_>>(),
        ),
        (
            "SEQTYPE",
            rows.iter().map(|a| a.sequence.as_str()).collect(),
        ),
        ("SIZECLASS", rows.iter().map(|a| a.size.as_str()).collect()),
        ("NATIVE", rows.iter().map(|a| a.native.as_str()).collect()),
        ("INVBASIS", rows.iter().map(|a| a.basis.as_str()).collect()),
    ]
    .iter()
    .map(|(k, v)| format!("{k}={}", v.join(",")))
    .collect::<Vec<_>>()
    .join(";")
}

pub struct Table(pub Vec<Vec<Option<Annotation>>>);
impl Table {
    pub fn read(output: &Path, data: &Catalogue, reference: Option<&str>) -> Result<Self> {
        let schema: serde_json::Value =
            serde_json::from_reader(io::reader(&output.join("schema.json"))?)
                .map_err(|e| e.to_string())?;
        if schema["version"] != crate::variant::OUTPUT_SCHEMA
            || schema["classification"]["version"] != 2
        {
            return Err("Unsupported variant annotation schema".into());
        }
        let minimum = schema["classification"]["minimum_inversion_span"]
            .as_u64()
            .ok_or("Missing inversion threshold")? as usize;
        if minimum < 19 {
            return Err("Invalid inversion threshold".into());
        }
        let classifier = reference.map(|r| Classifier::new(r, minimum)).transpose()?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(io::reader(&output.join(FILE))?);
        if !reader
            .headers()
            .map_err(|e| e.to_string())?
            .iter()
            .eq(COLUMNS)
        {
            return Err("Invalid variant annotation header".into());
        }
        let mut rows: Vec<Vec<Option<Annotation>>> =
            data.alleles.iter().map(|a| vec![None; a.len()]).collect();
        for record in reader.records() {
            let r = record.map_err(|e| e.to_string())?;
            let li = *data.index.get(&r[0]).ok_or("Unknown annotated locus")?;
            let ai: usize = r[1].parse().map_err(|_| "Invalid annotated allele index")?;
            if ai == 0 || ai >= rows[li].len() || rows[li][ai].is_some() {
                return Err("Duplicate or invalid annotated allele".into());
            }
            let mut segments = vec![];
            if &r[7] != "." {
                for text in r[7].split('|') {
                    let p: Vec<_> = text.split(':').collect();
                    if p.len() != 5 {
                        return Err("Invalid inversion witness encoding".into());
                    }
                    let v: Vec<usize> = p[..4]
                        .iter()
                        .map(|x| {
                            x.parse()
                                .map_err(|_| "Invalid inversion witness coordinate".into())
                        })
                        .collect::<Result<_>>()?;
                    segments.push(Segment {
                        start: v[0],
                        end: v[1],
                        query: v[2],
                        query_end: v[3],
                        template: p[4].into(),
                    });
                }
            }
            let a = Annotation {
                kind: r[2].into(),
                sequence: r[3].into(),
                size: r[4].into(),
                native: r[5].into(),
                basis: r[6].into(),
                segments,
            };
            let length = schema["reference"]["length"]
                .as_u64()
                .ok_or("Invalid schema reference length")? as usize;
            if a.segments.iter().any(|s| s.end > length) {
                return Err("Inversion template outside reference".into());
            }
            let (s, e) = data.bounds[li];
            a.verify(s, &data.alleles[li][0], &data.alleles[li][ai], minimum)?;
            if let Some(c) = &classifier {
                if a != c.annotate(s, e, &data.alleles[li][ai], &a.native)? {
                    return Err("Annotation differs from full reference evidence".into());
                }
            }
            rows[li][ai] = Some(a);
        }
        if rows.iter().any(|r| r.iter().skip(1).any(Option::is_none)) {
            return Err("Missing catalogue allele annotation".into());
        }
        Ok(Self(rows))
    }
    pub fn get(&self, locus: usize, allele: usize) -> &Annotation {
        self.0[locus][allele].as_ref().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn witness_validation_rejects_corrupted_sequence_bounds_and_type() {
        let mut state = 28391u64;
        let reference: String = (0..500)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize] as char
            })
            .collect();
        let alternate = io::reverse_complement(&reference[104..304]);
        let original = Classifier::new(&reference, 31)
            .unwrap()
            .annotate(100, 300, &alternate, "SEQUENCE")
            .unwrap();
        assert_eq!(original.kind, "INV_COMPLEX");
        for change in 0..4 {
            let mut bad = original.clone();
            match change {
                0 => bad.segments[0].end += 1,
                1 => bad.segments[0].template.replace_range(0..1, "N"),
                2 => bad.kind = "INV".into(),
                _ => bad.segments.clear(),
            }
            assert!(bad
                .verify(100, &reference[100..300], &alternate, 31)
                .is_err());
        }
    }
}
