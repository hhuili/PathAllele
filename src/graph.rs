//! Strict sequence-resolved GFA S/L/P graphs, retaining external segment names.
use crate::{
    io,
    model::{dna, Result},
};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;

pub type NodeId = i32;
pub struct Node {
    pub name: String,
    pub sequence: String,
}
pub struct Walk {
    pub name: String,
    pub sample: Option<usize>,
    pub nodes: Vec<NodeId>,
    pub length: usize,
}
pub struct Graph {
    pub chrom: String,
    pub reference: String,
    pub samples: Vec<String>,
    pub nodes: Vec<Node>,
    pub paths: Vec<Walk>,
    pub sample_paths: Vec<Vec<usize>>,
    pub reference_path: usize,
    pub ref_offsets: Vec<usize>,
    pub ref_occurrences: HashMap<(NodeId, usize), usize>,
    pub ref_places: HashMap<NodeId, Vec<usize>>,
    inventory: Vec<(String, String, usize, String)>,
}

fn valid_name(s: &str) -> bool {
    !s.is_empty()
        && !s.starts_with(['*', '='])
        && s.bytes().all(|b| (33..=126).contains(&b))
        && !s.contains("+,")
        && !s.contains("-,")
}
fn rows(path: &Path, mut callback: impl FnMut(usize, &[&str]) -> Result<()>) -> Result<()> {
    let mut input = io::reader(path)?;
    let mut line = String::new();
    let mut number = 0;
    loop {
        line.clear();
        if input.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
            break;
        }
        number += 1;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
        let min = match fields[0] {
            "H" => 1,
            "S" => 3,
            "L" => 6,
            "P" => 4,
            _ => {
                return Err(format!(
                    "GFA line {number}: unsupported record {}; sequence-resolved S/L/P only",
                    fields[0]
                ))
            }
        };
        if fields.len() < min {
            return Err(format!(
                "GFA line {number}: {} requires {min} fields",
                fields[0]
            ));
        }
        callback(number, &fields).map_err(|e| format!("GFA line {number}: {e}"))?;
    }
    Ok(())
}
fn handle(
    name: &str,
    names: &mut HashMap<String, NodeId>,
    nodes: &mut Vec<Node>,
) -> Result<NodeId> {
    if let Some(&id) = names.get(name) {
        return Ok(id);
    }
    let id = i32::try_from(nodes.len()).map_err(|_| "Too many graph segments")?;
    nodes.push(Node {
        name: name.into(),
        sequence: String::new(),
    });
    names.insert(name.into(), id);
    Ok(id)
}

impl Graph {
    pub fn read(
        gfa: &Path,
        reference: &Path,
        sample_file: &Path,
        reference_name: Option<&str>,
        path_map: Option<&Path>,
    ) -> Result<Self> {
        let fasta = io::read_fasta(reference)?;
        if fasta.len() != 1 {
            return Err("Exactly one reference FASTA record is supported".into());
        }
        let (chrom, reference) = fasta.into_iter().next().unwrap();
        if !dna(&reference)
            || chrom
                .bytes()
                .any(|b| b",<>\"".contains(&b) || !(33..=126).contains(&b))
        {
            return Err("Invalid reference sequence/VCF contig identifier".into());
        }
        let mut samples = Vec::new();
        let mut sample_index = HashMap::new();
        for line in io::reader(sample_file)?.lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let name = line.trim();
            if name == chrom
                || name.chars().any(char::is_whitespace)
                || sample_index.contains_key(name)
            {
                return Err("Samples must be unique non-reference haploid identifiers".into());
            }
            sample_index.insert(name.to_string(), samples.len());
            samples.push(name.to_string());
        }
        if samples.is_empty() {
            return Err("Samples must be nonempty".into());
        }
        let explicit = if let Some(path) = path_map {
            let mut input = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_reader(io::reader(path)?);
            let headers = input.headers().map_err(|e| e.to_string())?.clone();
            let pi = headers
                .iter()
                .position(|x| x == "path")
                .ok_or("Path map needs path/sample columns")?;
            let si = headers
                .iter()
                .position(|x| x == "sample")
                .ok_or("Path map needs path/sample columns")?;
            let mut map = HashMap::new();
            for row in input.records() {
                let row = row.map_err(|e| e.to_string())?;
                let sample = sample_index
                    .get(&row[si])
                    .ok_or("Undeclared sample in path map")?;
                if map.insert(row[pi].to_string(), *sample).is_some() {
                    return Err("Duplicate explicit sample path".into());
                }
            }
            Some(map)
        } else {
            None
        };
        let mut nodes = vec![Node {
            name: String::new(),
            sequence: String::new(),
        }];
        let mut names = HashMap::new();
        let mut namespace = HashSet::new();
        let mut edges = HashSet::new();
        rows(gfa, |_, f| {
            match f[0] {
                "H" => {
                    for tag in &f[1..] {
                        if let Some(v) = tag.strip_prefix("VN:Z:") {
                            if !["1.0", "1.1", "1.2"].contains(&v) {
                                return Err("Unsupported GFA version".into());
                            }
                        }
                    }
                }
                "S" | "P" => {
                    if !valid_name(f[1]) || !namespace.insert(f[1].to_string()) {
                        return Err("Invalid or duplicate segment/path name".into());
                    }
                }
                _ => {}
            }
            if f[0] == "S" {
                let id = handle(f[1], &mut names, &mut nodes)? as usize;
                let sequence = f[2].to_ascii_uppercase();
                if sequence.is_empty() || !io::iupac(&sequence) {
                    return Err("Missing or invalid node sequence".into());
                }
                for tag in &f[3..] {
                    if let Some(v) = tag.strip_prefix("LN:i:") {
                        if v.parse::<usize>().ok() != Some(sequence.len()) {
                            return Err("LN differs from segment sequence length".into());
                        }
                    }
                }
                nodes[id].sequence = sequence;
            } else if f[0] == "L" {
                if f[5] != "0M" {
                    return Err("Link overlap must be explicit 0M".into());
                }
                if !valid_name(f[1])
                    || !valid_name(f[3])
                    || !["+", "-"].contains(&f[2])
                    || !["+", "-"].contains(&f[4])
                {
                    return Err("Invalid link names/orientation".into());
                }
                let a = handle(f[1], &mut names, &mut nodes)? * if f[2] == "+" { 1 } else { -1 };
                let b = handle(f[3], &mut names, &mut nodes)? * if f[4] == "+" { 1 } else { -1 };
                edges.insert((a, b));
                edges.insert((-b, -a));
            }
            Ok(())
        })?;
        if nodes[1..].iter().any(|n| n.sequence.is_empty()) {
            return Err("GFA links reference undefined segments".into());
        }
        let mut tokens = HashMap::new();
        for (name, &id) in &names {
            tokens.insert(format!("{name}+"), id);
            tokens.insert(format!("{name}-"), -id);
        }
        let mut paths = Vec::new();
        let mut sample_paths = vec![Vec::new(); samples.len()];
        let mut inventory = Vec::new();
        let mut candidates = Vec::new();
        let mut haps = vec![HashSet::new(); samples.len()];
        let mut selected_names = HashSet::new();
        rows(gfa, |_, f| {
            if f[0] != "P" {
                return Ok(());
            }
            let name = f[1];
            let sample = if let Some(map) = &explicit {
                map.get(name).copied()
            } else {
                sample_index.get(name.split('#').next().unwrap()).copied()
            };
            let is_ref = if let Some(r) = reference_name {
                name == r
            } else {
                name == chrom || name.starts_with(&format!("{chrom}#"))
            };
            let selected = sample.is_some() || is_ref;
            if f[2].contains(';') {
                return Err("Jump-connected paths are unsupported".into());
            }
            let mut walk = Vec::new();
            let mut previous = None;
            let mut length = 0usize;
            let mut count = 0;
            for token in f[2].split(',') {
                let &id = tokens
                    .get(token)
                    .ok_or_else(|| format!("Unknown node or malformed path step {token}"))?;
                if let Some(prior) = previous {
                    if !edges.contains(&(prior, id)) {
                        return Err(format!("Path edge absent from graph: {name}"));
                    }
                }
                length = length
                    .checked_add(nodes[id.unsigned_abs() as usize].sequence.len())
                    .ok_or("Path length overflow")?;
                count += 1;
                previous = Some(id);
                if selected {
                    walk.push(id);
                }
            }
            if f[3] != "*" {
                let overlap: Vec<_> = f[3].split(',').collect();
                if overlap.len() != count - 1 || overlap.iter().any(|&s| s != "0M") {
                    return Err("P overlaps need exactly one 0M per edge".into());
                }
            }
            let sample_name = sample.map(|si| samples[si].clone()).unwrap_or_else(|| {
                if explicit.is_some() {
                    ".".into()
                } else {
                    name.split('#').next().unwrap().into()
                }
            });
            inventory.push((
                name.into(),
                sample_name,
                length,
                if is_ref {
                    "reference"
                } else if selected {
                    "selected"
                } else {
                    "not_selected"
                }
                .into(),
            ));
            if !selected {
                return Ok(());
            }
            selected_names.insert(name.to_string());
            let id = paths.len();
            if is_ref {
                if explicit.as_ref().is_some_and(|map| map.contains_key(name)) {
                    return Err("Reference cannot be assigned to a genotyped sample".into());
                }
                candidates.push(id);
            } else {
                let si = sample.unwrap();
                sample_paths[si].push(id);
                if explicit.is_none() {
                    if let Some(hap) = name.split('#').nth(1) {
                        haps[si].insert(hap.to_string());
                    }
                }
            }
            paths.push(Walk {
                name: name.into(),
                sample: if is_ref { None } else { sample },
                nodes: walk,
                length,
            });
            Ok(())
        })?;
        if candidates.len() != 1 {
            return Err("Select one complete GFA reference path with --reference-path".into());
        }
        if sample_paths.iter().any(Vec::is_empty) {
            return Err("Samples have no selected graph paths".into());
        }
        if haps.iter().any(|h| h.len() > 1) {
            return Err("Multiple haplotypes per sample are unsupported".into());
        }
        if explicit
            .as_ref()
            .is_some_and(|m| m.keys().any(|n| !selected_names.contains(n)))
        {
            return Err("Explicit sample paths are absent from graph".into());
        }
        let reference_path = candidates[0];
        let ref_offsets = offsets(&nodes, &paths[reference_path].nodes);
        let ref_occurrences = paths[reference_path]
            .nodes
            .iter()
            .enumerate()
            .map(|(i, &n)| ((n, ref_offsets[i]), i))
            .collect();
        let mut ref_places: HashMap<NodeId, Vec<usize>> = HashMap::new();
        for (i, &n) in paths[reference_path].nodes.iter().enumerate() {
            ref_places.entry(n.abs()).or_default().push(i);
        }
        let graph = Self {
            chrom,
            reference,
            samples,
            nodes,
            paths,
            sample_paths,
            reference_path,
            ref_offsets,
            ref_occurrences,
            ref_places,
            inventory,
        };
        if graph.spell(reference_path) != graph.reference {
            return Err("Spelled graph reference differs from FASTA".into());
        }
        Ok(graph)
    }
    pub fn refwalk(&self) -> &[NodeId] {
        &self.paths[self.reference_path].nodes
    }
    pub fn spell(&self, path: usize) -> String {
        spell(&self.nodes, &self.paths[path].nodes)
    }
    pub fn write_inventory(&self, output: &Path) -> Result<()> {
        let mut out = io::tsv(&output.join("Graph.nodes.tsv.gz"))?;
        out.write_record(["internal_id", "segment_id", "length", "sequence_sha256"])
            .map_err(|e| e.to_string())?;
        for (i, n) in self.nodes.iter().enumerate().skip(1) {
            out.write_record([
                i.to_string(),
                n.name.clone(),
                n.sequence.len().to_string(),
                io::hash_bytes(n.sequence.as_bytes()),
            ])
            .map_err(|e| e.to_string())?;
        }
        io::finish_tsv(out)?;
        let mut out = io::tsv(&output.join("Graph.input_paths.tsv.gz"))?;
        out.write_record(["path", "sample", "length", "selection"])
            .map_err(|e| e.to_string())?;
        for (name, sample, len, status) in &self.inventory {
            out.write_record([
                name.clone(),
                sample.clone(),
                len.to_string(),
                status.clone(),
            ])
            .map_err(|e| e.to_string())?;
        }
        io::finish_tsv(out)
    }
}
pub fn offsets(nodes: &[Node], walk: &[NodeId]) -> Vec<usize> {
    let mut result = Vec::with_capacity(walk.len() + 1);
    result.push(0);
    for &n in walk {
        result.push(result.last().unwrap() + nodes[n.unsigned_abs() as usize].sequence.len());
    }
    result
}
pub fn spell(nodes: &[Node], walk: &[NodeId]) -> String {
    let mut result = String::with_capacity(
        walk.iter()
            .map(|n| nodes[n.unsigned_abs() as usize].sequence.len())
            .sum(),
    );
    for &n in walk {
        let s = &nodes[n.unsigned_abs() as usize].sequence;
        if n > 0 {
            result.push_str(s)
        } else {
            result.push_str(&io::reverse_complement(s))
        }
    }
    result
}
