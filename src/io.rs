//! Streaming files and atomic, versioned work checkpoints.
use crate::model::Result;
use bincode::Options;
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use serde::{de::DeserializeOwned, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
static LOG: std::sync::Mutex<Option<File>> = std::sync::Mutex::new(None);
pub fn set_log(path: &Path) -> Result<()> {
    *LOG.lock().map_err(|e| e.to_string())? = Some(
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .map_err(|e| e.to_string())?,
    );
    Ok(())
}
pub fn progress(message: impl AsRef<str>) {
    let text = message.as_ref();
    eprintln!("{text}");
    if let Ok(mut guard) = LOG.lock() {
        if let Some(log) = guard.as_mut() {
            let _ = writeln!(log, "{text}");
        }
    }
}

pub fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if path.extension().is_some_and(|s| s == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}
pub enum Output {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}
impl Write for Output {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(bytes),
            Self::Gzip(w) => w.write(bytes),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
        }
    }
}
impl Output {
    /// Complete compression and surface footer/flush errors before recording hashes.
    pub fn finish(self) -> Result<()> {
        let writer = match self {
            Self::Plain(w) => w,
            Self::Gzip(w) => w.finish().map_err(|e| e.to_string())?,
        };
        let file = writer.into_inner().map_err(|e| e.to_string())?;
        file.sync_all().map_err(|e| e.to_string())
    }
}
pub fn finish_tsv(writer: csv::Writer<Output>) -> Result<()> {
    writer.into_inner().map_err(|e| e.to_string())?.finish()
}
pub fn writer(path: &Path) -> Result<Output> {
    let file = File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if path.extension().is_some_and(|s| s == "gz") {
        Ok(Output::Gzip(GzEncoder::new(
            BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Output::Plain(BufWriter::new(file)))
    }
}
pub fn read_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    let mut records: Vec<(String, String)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for line in reader(path)?.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let text = line.trim();
        if text.is_empty() {
            continue;
        }
        if let Some(header) = text.strip_prefix('>') {
            let name = header
                .split_whitespace()
                .next()
                .ok_or("Empty FASTA header")?;
            if !seen.insert(name.to_string()) {
                return Err("Duplicate FASTA identifier".into());
            }
            records.push((name.to_string(), String::new()));
        } else {
            let row = records
                .last_mut()
                .ok_or("FASTA sequence precedes its header")?;
            let seq = text.to_ascii_uppercase();
            if !iupac(&seq) {
                return Err(format!("Unsupported FASTA symbols: {}", row.0));
            }
            row.1.push_str(&seq);
        }
    }
    if records.is_empty() || records.iter().any(|r| r.1.is_empty()) {
        return Err("FASTA must contain nonempty records".into());
    }
    Ok(records)
}
pub fn iupac(s: &str) -> bool {
    s.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
}
pub fn reverse_complement(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            _ => 'N',
        })
        .collect()
}
pub fn hash_bytes(value: &[u8]) -> String {
    format!("{:x}", Sha256::digest(value))
}
pub fn digest(path: &Path) -> Result<String> {
    let mut file = File::open(path).map_err(|e| e.to_string())?;
    let mut value = Sha256::new();
    let mut buffer = vec![0; 1 << 20];
    loop {
        let n = file.read(&mut buffer).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        value.update(&buffer[..n]);
    }
    Ok(format!("{:x}", value.finalize()))
}
pub fn json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?;
    atomic(path, &bytes)
}
pub fn atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let temp = path.with_extension("tmp");
    {
        let mut out = File::create(&temp).map_err(|e| e.to_string())?;
        out.write_all(bytes).map_err(|e| e.to_string())?;
        out.sync_all().map_err(|e| e.to_string())?;
    }
    std::fs::rename(temp, path).map_err(|e| e.to_string())
}
pub fn checkpoint<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let temporary = path.with_extension("tmp");
    let mut out = BufWriter::new(File::create(&temporary).map_err(|e| e.to_string())?);
    out.write_all(b"PATHALLELE\0WORK\x02")
        .map_err(|e| e.to_string())?;
    bincode::serialize_into(&mut out, value).map_err(|e| e.to_string())?;
    out.flush().map_err(|e| e.to_string())?;
    out.get_ref().sync_all().map_err(|e| e.to_string())?;
    drop(out);
    std::fs::rename(temporary, path).map_err(|e| e.to_string())
}
pub fn load_checkpoint<T: DeserializeOwned>(path: &Path) -> Result<T> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let length = file.metadata().map_err(|e| e.to_string())?.len();
    let mut data = BufReader::new(file);
    let prefix = b"PATHALLELE\0WORK\x02";
    let mut magic = [0; b"PATHALLELE\0WORK\x02".len()];
    data.read_exact(&mut magic).map_err(|e| e.to_string())?;
    if magic != *prefix {
        return Err("Unsupported work checkpoint schema".into());
    }
    let value = bincode::DefaultOptions::new()
        .with_fixint_encoding()
        .with_limit(length)
        .deserialize_from(&mut data)
        .map_err(|e| e.to_string())?;
    if data.read(&mut [0u8]).map_err(|e| e.to_string())? != 0 {
        return Err("Trailing data in work checkpoint".into());
    }
    Ok(value)
}
pub fn tsv(path: &Path) -> Result<csv::Writer<Output>> {
    Ok(csv::WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(csv::Terminator::Any(b'\n'))
        .from_writer(writer(path)?))
}
