//! Complete haplotypes in indexed VCF; temporary numeric mappings are private.
use crate::{
    classify, io,
    model::{self, dna, Result},
    project::{ProjectionRow, STATUS},
    variant,
};
use memmap2::{Mmap, MmapOptions};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

pub const FILE: &str = "Graph.loci.vcf.gz";

pub struct Temporary {
    pub path: PathBuf,
}
impl Temporary {
    pub fn new(root: &Path, prefix: &str) -> Result<Self> {
        static ID: AtomicU64 = AtomicU64::new(0);
        for _ in 0..100 {
            let path = root.join(format!(
                ".{prefix}-{}-{}",
                std::process::id(),
                ID.fetch_add(1, Ordering::Relaxed)
            ));
            let mut builder = fs::DirBuilder::new();
            #[cfg(unix)]
            {
                use std::os::unix::fs::DirBuilderExt;
                builder.mode(0o700);
            }
            match builder.create(&path) {
                Ok(()) => return Ok(Self { path }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e.to_string()),
            }
        }
        Err("Cannot create private temporary directory".into())
    }
}
impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// One record is buffered at a time. Sorting only keeps file offsets in memory.
pub struct Writer<'a> {
    temporary: Temporary,
    stream: BufWriter<File>,
    order: Vec<(usize, String, u64, u64)>,
    offset: u64,
    samples: usize,
    reference: &'a str,
    chrom: &'a str,
    header: String,
    table: csv::Writer<io::Output>,
    output: PathBuf,
}
impl<'a> Writer<'a> {
    pub fn new(
        output: &Path,
        samples: &[String],
        chrom: &'a str,
        reference: &'a str,
    ) -> Result<Self> {
        let temporary = Temporary::new(output, "catalogue")?;
        let stream = BufWriter::new(
            File::create(temporary.path.join("rows.vcf")).map_err(|e| e.to_string())?,
        );
        let mut table = io::tsv(&output.join("Graph.loci.tsv.gz"))?;
        table
            .write_record([
                "locus_id",
                "allele_index",
                "ref_start0",
                "ref_end0",
                "sequence",
            ])
            .map_err(|e| e.to_string())?;
        let mut header = vec!["##fileformat=VCFv4.2".into(), format!("##source=PathAllele-{}", env!("CARGO_PKG_VERSION")), format!("##contig=<ID={chrom},length={}>", reference.len()),
            "##pathallele_locus_scope=\"Independent overlapping complete haplotypes, including reference-only loci; do not replay all records together\"".into(),
            "##INFO=<ID=LSTART,Number=1,Type=Integer,Description=\"Unpadded complete locus start, zero-based inclusive\">".into(),
            "##INFO=<ID=LEND,Number=1,Type=Integer,Description=\"Unpadded complete locus end, zero-based exclusive\">".into(),
            "##INFO=<ID=LKIND,Number=1,Type=String,Description=\"Catalogue locus kind: SNP or EDIT\">".into(),
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Complete ALT allele count\">".into(),
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Called complete haploid genotypes\">".into(),
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Complete ALT allele frequency\">".into(),
            "##FORMAT=<ID=PG,Number=1,Type=Integer,Description=\"Projected genotype indexing the complete catalogue, not Graph.vcf.gz ALT order; dot if unresolved\">".into(),
            "##FORMAT=<ID=PF,Number=1,Type=String,Description=\"Projected genotype status; definitions in pathallele_status metadata\">".into()];
        header.extend(variant::headers("complete_locus"));
        header.extend(
            variant::CLASS_HEADERS
                .iter()
                .filter(|s| !s.starts_with("##pathallele_hierarchy="))
                .map(|&s| s.to_string()),
        );
        let header = format!(
            "{}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            header.join("\n"),
            samples.join("\t")
        );
        Ok(Self {
            temporary,
            stream,
            order: vec![],
            offset: 0,
            samples: samples.len(),
            reference,
            chrom,
            header,
            table,
            output: output.into(),
        })
    }
    pub fn add(
        &mut self,
        id: &str,
        s: usize,
        e: usize,
        kind: &str,
        row: &ProjectionRow,
        classes: &[classify::Annotation],
    ) -> Result<()> {
        if row.alleles.is_empty()
            || classes.len() + 1 != row.alleles.len()
            || [&row.gt, &row.projected]
                .iter()
                .any(|g| g.len() != self.samples)
            || [&row.source, &row.export]
                .iter()
                .any(|s| s.len() != self.samples)
        {
            return Err("Invalid complete locus row dimensions".into());
        }
        if self.reference.get(s..e) != Some(row.alleles[0].as_str())
            || e > i32::MAX as usize
            || !matches!(kind, "SNP" | "EDIT")
        {
            return Err("Invalid complete locus reference/bounds/kind".into());
        }
        for (ai, seq) in row.alleles.iter().enumerate() {
            self.table
                .write_record([
                    id.to_string(),
                    ai.to_string(),
                    s.to_string(),
                    e.to_string(),
                    seq.clone(),
                ])
                .map_err(|e| e.to_string())?;
        }
        let (pos, reference, alts) =
            model::vcf_alleles(self.reference, s, e, &row.alleles[1..], kind == "SNP")?;
        let mut ac = vec![0usize; alts.len()];
        let mut an = 0usize;
        let mut cells = String::new();
        for si in 0..self.samples {
            let (g, p, ft, pf) = (
                row.gt[si],
                row.projected[si],
                row.source[si] as usize,
                row.export[si] as usize,
            );
            for (v, status) in [(g, ft), (p, pf)] {
                if v < -1
                    || v as i64 >= row.alleles.len() as i64
                    || status >= STATUS.len()
                    || (v >= 0) != (status == 0)
                {
                    return Err("Genotype index/status disagreement".into());
                }
            }
            if p >= 0 && g < 0 {
                return Err("Projection invents observation".into());
            }
            if g >= 0 {
                an += 1;
                if g > 0 {
                    ac[g as usize - 1] += 1;
                }
            }
            use std::fmt::Write as _;
            cells.push('\t');
            write_genotype(&mut cells, g)?;
            write!(cells, ":{}:", STATUS[ft]).map_err(|e| e.to_string())?;
            write_genotype(&mut cells, p)?;
            write!(cells, ":{}", STATUS[pf]).map_err(|e| e.to_string())?;
        }
        let mut info = format!("LSTART={s};LEND={e};LKIND={kind};AN={an}");
        if !alts.is_empty() {
            let counts = ac
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let af = ac
                .iter()
                .map(|&v| {
                    if an == 0 {
                        ".".into()
                    } else {
                        format!("{:.8}", v as f64 / an as f64)
                    }
                })
                .collect::<Vec<_>>()
                .join(",");
            info.push_str(&format!(
                ";AC={counts};AF={af};{}",
                classify::info(&classes.iter().collect::<Vec<_>>())
            ));
        }
        let line = format!(
            "{}\t{pos}\t{id}\t{reference}\t{}\t.\t.\t{info}\tGT:FT:PG:PF{cells}\n",
            self.chrom,
            if alts.is_empty() {
                ".".into()
            } else {
                alts.join(",")
            }
        );
        self.stream
            .write_all(line.as_bytes())
            .map_err(|e| e.to_string())?;
        self.order
            .push((pos, id.into(), self.offset, line.len() as u64));
        self.offset = self
            .offset
            .checked_add(line.len() as u64)
            .ok_or("Locus spool length overflow")?;
        Ok(())
    }
    pub fn finish(mut self, threads: usize) -> Result<()> {
        io::finish_tsv(self.table)?;
        self.stream.flush().map_err(|e| e.to_string())?;
        self.order
            .sort_unstable_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
        if self
            .order
            .iter()
            .map(|r| &r.1)
            .collect::<HashSet<_>>()
            .len()
            != self.order.len()
        {
            return Err("Duplicate locus ID".into());
        }
        let mut input = BufReader::new(
            File::open(self.temporary.path.join("rows.vcf")).map_err(|e| e.to_string())?,
        );
        let pending = self.temporary.path.join(FILE);
        let mut out = noodles_bgzf::MultithreadedWriter::with_worker_count(
            std::num::NonZeroUsize::new(threads).ok_or("Compression threads must be positive")?,
            File::create(&pending).map_err(|e| e.to_string())?,
        );
        out.write_all(self.header.as_bytes())
            .map_err(|e| e.to_string())?;
        for (_, _, offset, length) in self.order {
            input
                .seek(SeekFrom::Start(offset))
                .map_err(|e| e.to_string())?;
            let copied = std::io::copy(&mut (&mut input).take(length), &mut out)
                .map_err(|e| e.to_string())?;
            if copied != length {
                return Err("Truncated locus spool".into());
            }
        }
        out.finish()
            .map_err(|e| e.to_string())?
            .sync_all()
            .map_err(|e| e.to_string())?;
        fs::rename(pending, self.output.join(FILE)).map_err(|e| e.to_string())
    }
}
fn write_genotype(out: &mut String, g: i32) -> Result<()> {
    use std::fmt::Write;
    if g < 0 {
        out.push('.');
        Ok(())
    } else {
        write!(out, "{g}").map_err(|e| e.to_string())
    }
}

struct Matrices {
    map: Option<Mmap>,
    _temporary: Temporary,
}
impl Matrices {
    fn integer(&self, offset: usize) -> i32 {
        i32::from_le_bytes(
            self.map.as_ref().unwrap()[offset..offset + 4]
                .try_into()
                .unwrap(),
        )
    }
    fn byte(&self, offset: usize) -> usize {
        self.map.as_ref().unwrap()[offset] as usize
    }
}

pub struct Catalogue {
    data: Matrices,
    pub samples: Vec<String>,
    pub loci: Vec<String>,
    pub kinds: Vec<String>,
    pub bounds: Vec<(usize, usize)>,
    pub alleles: Vec<Vec<String>>,
    pub index: HashMap<String, usize>,
}
impl Catalogue {
    pub fn open(output: &Path, length: usize) -> Result<Self> {
        let mut reader = variant::Reader::open(&output.join(FILE))?;
        let mut metadata = vec![];
        noodles_vcf::io::Writer::new(&mut metadata)
            .write_header(&reader.header)
            .map_err(|e| e.to_string())?;
        let metadata = String::from_utf8(metadata).map_err(|e| e.to_string())?;
        for required in [
            "##pathallele_genotype_contract=complete_locus\n",
            format!(
                "##pathallele_output_schema={}\n",
                crate::variant::OUTPUT_SCHEMA
            )
            .as_str(),
            "##pathallele_variant_schema=2\n",
        ] {
            if !metadata.contains(required) {
                return Err("Wrong complete-locus VCF contract/schema".into());
            }
        }
        if reader.header.contigs().len() != 1 {
            return Err("Expected one complete-locus contig".into());
        }
        let chrom = reader.header.contigs().keys().next().unwrap().clone();
        let samples: Vec<_> = reader.header.sample_names().iter().cloned().collect();
        if samples.is_empty()
            || samples
                .iter()
                .any(|s| s.is_empty() || s.chars().any(char::is_whitespace))
            || samples.iter().collect::<HashSet<_>>().len() != samples.len()
        {
            return Err("Invalid/duplicate sample identities".into());
        }
        let temporary = Temporary::new(&std::env::temp_dir(), "pathallele-matrix")?;
        let path = temporary.path.join("rows");
        let mut file = BufWriter::new(File::create(&path).map_err(|e| e.to_string())?);
        let (mut loci, mut kinds, mut bounds, mut alleles) = (vec![], vec![], vec![], vec![]);
        let mut index = HashMap::new();
        let mut previous_pos = 0;
        while let Some(r) = reader.next_record()? {
            if r.chrom != chrom {
                return Err("Unexpected complete-locus contig".into());
            }
            let get = |key: &str| {
                r.info
                    .get(key)
                    .ok_or_else(|| format!("Missing locus INFO/{key}"))
            };
            let s = get("LSTART")?
                .parse::<usize>()
                .map_err(|_| "Invalid LSTART")?;
            let e = get("LEND")?.parse::<usize>().map_err(|_| "Invalid LEND")?;
            let kind = get("LKIND")?;
            if s > e
                || e > length
                || !matches!(kind.as_str(), "SNP" | "EDIT")
                || r.pos < previous_pos
            {
                return Err("Invalid or unsorted complete locus bounds/kind".into());
            }
            previous_pos = r.pos;
            if index.insert(r.id.clone(), loci.len()).is_some() {
                return Err("Duplicate locus ID".into());
            }
            let padded: Vec<_> = std::iter::once(r.reference.clone())
                .chain(r.alts.iter().cloned())
                .collect();
            let mut seqs = Vec::new();
            for seq in &padded {
                if seq.is_empty() || !dna(seq) {
                    return Err("Invalid explicit complete allele".into());
                }
                let value = if kind == "SNP" {
                    seq.as_str()
                } else if s > 0 {
                    &seq[1..]
                } else {
                    &seq[..seq.len() - 1]
                };
                seqs.push(value.to_string());
            }
            if seqs[0].len() != e - s
                || seqs.iter().collect::<HashSet<_>>().len() != seqs.len()
                || (kind == "SNP" && (e - s != 1 || seqs.iter().any(|s| s.len() != 1)))
            {
                return Err("Invalid complete sequence dimensions".into());
            }
            let expected_pos = if kind == "SNP" { s + 1 } else { s.max(1) };
            if r.pos != expected_pos
                || (kind != "SNP"
                    && padded.iter().any(|a| {
                        if s > 0 {
                            a.as_bytes()[0] != padded[0].as_bytes()[0]
                        } else {
                            a.as_bytes().last() != padded[0].as_bytes().last()
                        }
                    }))
            {
                return Err("Invalid complete locus padding/POS".into());
            }
            let an = r.gt.iter().filter(|&&g| g >= 0).count();
            if get("AN")?.parse::<usize>().map_err(|_| "Invalid AN")? != an {
                return Err("Complete AN disagrees with GT".into());
            }
            if seqs.len() > 1 {
                let ac = get("AC")?
                    .split(',')
                    .map(|v| v.parse::<usize>().map_err(|_| "Invalid AC"))
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                if ac.len() + 1 != seqs.len() {
                    return Err("Wrong AC cardinality".into());
                }
                let mut observed = vec![0; ac.len()];
                for &g in &r.gt {
                    if g > 0 {
                        *observed
                            .get_mut(g as usize - 1)
                            .ok_or("GT outside catalogue")? += 1;
                    }
                }
                if ac != observed {
                    return Err("Complete AC disagrees with GT".into());
                }
            }
            let pg = r.projected.as_ref().ok_or("Missing FORMAT/PG")?;
            let pf = r.projection_status.as_ref().ok_or("Missing FORMAT/PF")?;
            let mut row = Vec::with_capacity(
                samples
                    .len()
                    .checked_mul(10)
                    .ok_or("Matrix row size overflow")?,
            );
            for values in [&r.gt, pg] {
                for &g in values {
                    row.extend(g.to_le_bytes());
                }
            }
            for (values, statuses) in [(&r.gt, &r.ft), (pg, pf)] {
                for (&g, status) in values.iter().zip(statuses) {
                    let code = STATUS
                        .iter()
                        .position(|s| s == status)
                        .ok_or("Unknown genotype status")?;
                    if g < -1 || g as i64 >= seqs.len() as i64 || (g >= 0) != (code == 0) {
                        return Err("Genotype index/status disagreement".into());
                    }
                    row.push(code as u8);
                }
            }
            if pg.iter().zip(&r.gt).any(|(&p, &g)| p >= 0 && g < 0) {
                return Err("Projection invents observation".into());
            }
            file.write_all(&row).map_err(|e| e.to_string())?;
            loci.push(r.id);
            kinds.push(kind.clone());
            bounds.push((s, e));
            alleles.push(seqs);
        }
        file.flush().map_err(|e| e.to_string())?;
        drop(file);
        let file = File::open(&path).map_err(|e| e.to_string())?;
        let expected = loci
            .len()
            .checked_mul(samples.len())
            .and_then(|n| n.checked_mul(10))
            .ok_or("Matrix size overflow")?;
        if file.metadata().map_err(|e| e.to_string())?.len() != expected as u64 {
            return Err("Truncated private matrix".into());
        }
        // The file is in a private directory, fully flushed, immutable for this
        // object's lifetime, and removed only after the mapping is dropped.
        let map = if expected == 0 {
            None
        } else {
            Some(unsafe { MmapOptions::new().map(&file) }.map_err(|e| e.to_string())?)
        };
        let data = Matrices {
            map,
            _temporary: temporary,
        };
        let mut table = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(io::reader(&output.join("Graph.loci.tsv.gz"))?);
        if !table.headers().map_err(|e| e.to_string())?.iter().eq([
            "locus_id",
            "allele_index",
            "ref_start0",
            "ref_end0",
            "sequence",
        ]) {
            return Err("Invalid locus table header".into());
        }
        let mut seen = vec![0; loci.len()];
        for row in table.records() {
            let row = row.map_err(|e| e.to_string())?;
            let li = *index
                .get(&row[0])
                .ok_or("Unknown locus in sequence table")?;
            let number = |i: usize| {
                row[i]
                    .parse::<usize>()
                    .map_err(|_| "Invalid locus table integer")
            };
            let ai = number(1)?;
            if ai != seen[li]
                || (number(2)?, number(3)?) != bounds[li]
                || alleles[li].get(ai).map(String::as_str) != Some(&row[4])
            {
                return Err("VCF/catalogue sequence disagreement".into());
            }
            seen[li] += 1;
        }
        if seen.iter().zip(&alleles).any(|(&n, a)| n != a.len()) {
            return Err("Incomplete locus sequence table".into());
        }
        Ok(Self {
            data,
            samples,
            loci,
            kinds,
            bounds,
            alleles,
            index,
        })
    }
    pub fn gt(&self, li: usize, si: usize) -> i32 {
        self.data.integer(li * self.samples.len() * 10 + si * 4)
    }
    pub fn projected(&self, li: usize, si: usize) -> i32 {
        self.data
            .integer(li * self.samples.len() * 10 + self.samples.len() * 4 + si * 4)
    }
    pub fn status(&self, li: usize, si: usize, export: bool) -> &str {
        STATUS[self.data.byte(
            li * self.samples.len() * 10 + self.samples.len() * if export { 9 } else { 8 } + si,
        )]
    }
}
