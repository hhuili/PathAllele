//! Physical occurrence controls retained across the standalone Rust migration.
use pathallele::{
    anchors::{MatchIndex, SeedCounts},
    catalogue::Temporary,
    context::{self, Alignment},
    graph::Graph,
    io,
};
use std::{collections::BTreeSet, fs};
fn fixture(reference: &[i32], query: &[i32]) -> (Temporary, Graph) {
    let t = Temporary::new(&std::env::temp_dir(), "pathallele-context").unwrap();
    let mut state = 401u64;
    let mut sequences: Vec<String> = [100, 2000, 100, 100, 1]
        .into_iter()
        .map(|n| {
            (0..n)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize] as char
                })
                .collect()
        })
        .collect();
    sequences.push("NNNN".into());
    let spell = |walk: &[i32]| {
        walk.iter()
            .map(|&n| {
                if n > 0 {
                    sequences[n as usize - 1].clone()
                } else {
                    io::reverse_complement(&sequences[-n as usize - 1])
                }
            })
            .collect::<String>()
    };
    let mut text = String::from("H\tVN:Z:1.0\n");
    for (i, seq) in sequences.iter().enumerate() {
        text += &format!("S\tn{}\t{seq}\n", i + 1);
    }
    let token = |n: i32| format!("n{}{}", n.abs(), if n > 0 { '+' } else { '-' });
    let edges: BTreeSet<_> = [reference, query]
        .into_iter()
        .flat_map(|w| w.windows(2).map(|p| (token(p[0]), token(p[1]))))
        .collect();
    for (a, b) in edges {
        text += &format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            &a[..a.len() - 1],
            &a[a.len() - 1..],
            &b[..b.len() - 1],
            &b[b.len() - 1..]
        );
    }
    for (name, walk) in [("ref#0#chr", reference), ("sample#1#ctg", query)] {
        text += &format!(
            "P\t{name}\t{}\t*\n",
            walk.iter().map(|&n| token(n)).collect::<Vec<_>>().join(",")
        );
    }
    fs::write(t.path.join("graph.gfa"), text).unwrap();
    fs::write(
        t.path.join("reference.fa"),
        format!(">ref\n{}\n", spell(reference)),
    )
    .unwrap();
    fs::write(t.path.join("samples.txt"), "sample\n").unwrap();
    let g = Graph::read(
        &t.path.join("graph.gfa"),
        &t.path.join("reference.fa"),
        Some(&t.path.join("samples.txt")),
        None,
        None,
    )
    .unwrap();
    (t, g)
}
fn alignment(id: &str, runs: &[(usize, usize, usize)]) -> Alignment {
    let (a, b) = (runs[0], runs[runs.len() - 1]);
    Alignment {
        identity: id.into(),
        reverse: false,
        qs: a.0,
        qe: b.1,
        rs: a.2,
        re: b.2 + b.1 - b.0,
        native_qs: a.0,
        native_qe: b.1,
        index: MatchIndex::new(runs.to_vec()),
    }
}
fn plan(g: &Graph, records: &[Alignment]) -> context::Plan {
    let counts = [
        SeedCounts::new(&g.reference, 19),
        SeedCounts::new(&g.reference, 31),
    ];
    context::plan(g, g.sample_paths[0][0], records, 31, 19, 1000, &counts).unwrap()
}
#[test]
fn separate_primary_records_cannot_manufacture_a_long_deletion() {
    let (_t, g) = fixture(&[1, 2, 3], &[1, 3]);
    let p = plan(
        &g,
        &[
            alignment("A", &[(0, 100, 0)]),
            alignment("B", &[(100, 200, 2100)]),
        ],
    );
    assert_eq!(p.blocks.len(), 2);
    assert!(p.missing.is_empty());
    assert!(p.blocks.iter().all(|b| b.runs.len() == 1));
    let (rows, mask) = context::structure_evidence(&p.blocks, 200, 2200).unwrap();
    assert_eq!(rows[0].4, "UNRESOLVED_ADJACENCY");
    assert_eq!(mask, vec![(100, 2100)]);
    let supported = plan(&g, &[alignment("A", &[(0, 100, 0), (100, 200, 2100)])]);
    assert_eq!(supported.blocks.len(), 1);
    assert_eq!(supported.blocks[0].runs.len(), 2);
    assert!(supported.missing.is_empty());
}
#[test]
fn graph_tail_recovery_requires_known_two_sided_sequence() {
    let (_t, g) = fixture(&[1, 3, 4], &[1, 3, 5, 4]);
    let p = plan(&g, &[alignment("A", &[(0, 100, 0)])]);
    assert!(p.missing.is_empty());
    assert!(p
        .blocks
        .iter()
        .any(|b| b.kind == "graph" && b.runs.len() == 2));
    let (_t, g) = fixture(&[1, 3], &[1, 5]);
    assert_eq!(plan(&g, &[]).missing, vec![(100, 101)]);
    let (_t, g) = fixture(&[1, 3], &[1, 6, 3]);
    let p = plan(&g, &[]);
    assert_eq!(p.missing, vec![(100, 104)]);
    assert!(p.blocks.iter().all(|b| b.runs.len() == 1));
}
#[test]
fn competing_primary_ownership_is_order_independent() {
    let (_t, g) = fixture(&[1, 3], &[1, 3]);
    let a = alignment("A", &[(0, 200, 0)]);
    let b = alignment("B", &[(50, 150, 50)]);
    let x = plan(&g, &[a.clone(), b.clone()]);
    let y = plan(&g, &[b, a]);
    assert_eq!(x.missing, y.missing);
    assert!(!x.missing.is_empty());
}
#[test]
fn repeated_reference_copy_does_not_establish_inversion_ownership() {
    let (_t, g) = fixture(&[1, 2, 2, 3], &[1, -2, 3]);
    assert!(plan(&g, &[]).blocks.iter().all(|b| b.kind != "inversion"));
}
#[test]
fn competing_offsets_remain_signed() {
    let index = MatchIndex::new(vec![(0, 100, 0), (10, 110, 0)]);
    assert_eq!(
        index.locations(0, 100).into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([-10, 0])
    );
}
