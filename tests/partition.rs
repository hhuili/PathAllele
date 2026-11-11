use pathallele::{
    alleles, anchors::Anchor, classify::Classifier, io, model, partition, source::Parent,
};

fn dna(n: usize, mut seed: u64) -> String {
    (0..n)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize] as char
        })
        .collect()
}

// Independently construct a collinear forward match chain, including the short
// coincidental matches that can occur inside a fragmented inversion graph.
fn matches(r: &str, a: &str) -> Vec<Anchor> {
    let (r, a) = (r.as_bytes(), a.as_bytes());
    let mut dp = vec![vec![0; a.len() + 1]; r.len() + 1];
    for i in (0..r.len()).rev() {
        for j in (0..a.len()).rev() {
            dp[i][j] = if r[i] == a[j] {
                dp[i + 1][j + 1] + 1
            } else {
                dp[i + 1][j].max(dp[i][j + 1])
            };
        }
    }
    let (mut i, mut j) = (0, 0);
    let mut out: Vec<Anchor> = vec![];
    while i < r.len() && j < a.len() {
        if r[i] == a[j] {
            if i > 0 && j > 0 && i + 1 < r.len() && j + 1 < a.len() {
                out.push(Anchor {
                    r: i,
                    q: j,
                    n: 1,
                    node: Some(1),
                });
            }
            i += 1;
            j += 1;
        } else if dp[i + 1][j] >= dp[i][j + 1] {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}
fn parent(reference: &str, alternate: &str, fragmented: bool) -> Parent {
    Parent {
        id: "synthetic".into(),
        start: 0,
        alleles: vec![reference.into(), alternate.into()],
        carriers: vec![vec![], vec![0]],
        anchors: vec![
            vec![],
            if fragmented {
                matches(reference, alternate)
            } else {
                vec![]
            },
        ],
        kinds: vec!["REFERENCE".into(), "SEQUENCE".into()],
    }
}
fn replay(plan: &partition::Plan) -> Vec<String> {
    let edits = alleles::decompose(plan).unwrap();
    plan.parents
        .iter()
        .zip(edits)
        .map(|(p, aa)| {
            let rebuilt = model::replay(
                &p.alleles[0],
                0,
                p.alleles[0].len(),
                &aa[1].iter().collect::<Vec<_>>(),
            )
            .unwrap();
            assert_eq!(rebuilt, p.alleles[1]);
            rebuilt
        })
        .collect()
}

#[test]
fn local_reversals_define_units_before_short_forward_matches() {
    for seed in [23, 1099, 1209, 98331] {
        let reference = dna(400, seed);
        let mut alt = reference.clone();
        alt.replace_range(100..200, &io::reverse_complement(&reference[100..200]));
        let snv = if reference.as_bytes()[330] == b'A' {
            "C"
        } else {
            "A"
        };
        alt.replace_range(330..331, snv);
        for fragmented in [false, true] {
            let p = parent(&reference, &alt, fragmented);
            let (full, partial) = partition::build(&[p], &reference, 31, 31).unwrap();
            for plan in [&full, &partial] {
                replay(plan);
                replay(&plan.reversed());
                let protected: Vec<_> = plan.units[0][1].iter().filter(|u| u.preserve).collect();
                assert_eq!(protected.len(), 1);
                let u = protected[0];
                assert!(u.r <= 100 && u.r + u.rn >= 200);
            }
            // Full projection retains an independent distant SNV; the repair
            // does not convert the complete parent into a giant replacement.
            let edits = alleles::decompose(&full).unwrap();
            assert!(edits[0][1].iter().any(|e| e.start == 330 && e.end == 331));
            let c = Classifier::new(&reference, 31).unwrap();
            assert!(edits[0][1].iter().any(|e| c
                .annotate(e.start, e.end, &e.alt, "SEQUENCE")
                .unwrap()
                .kind
                == "INV"));
            let twice = full.reversed().reversed();
            assert_eq!(replay(&twice), vec![alt.clone()]);
        }
    }
}

#[test]
fn separated_inversions_keep_the_forward_sequence_between_them() {
    let reference = dna(600, 23921);
    let mut alt = reference.clone();
    for (s, e) in [(100, 200), (350, 450)] {
        alt.replace_range(s..e, &io::reverse_complement(&reference[s..e]));
    }
    let (plan, _) =
        partition::build(&[parent(&reference, &alt, true)], &reference, 31, 31).unwrap();
    let units: Vec<_> = plan.units[0][1].iter().filter(|u| u.preserve).collect();
    assert_eq!(units.len(), 2);
    assert!(units[0].r + units[0].rn < units[1].r);
    replay(&plan);
}

#[test]
fn remote_copies_short_reversals_and_nonunique_evidence_do_not_protect_units() {
    let reference = dna(500, 811);
    let mut remote = reference.clone();
    remote.replace_range(100..200, &io::reverse_complement(&reference[300..400]));
    let mut short = reference.clone();
    short.replace_range(100..120, &io::reverse_complement(&reference[100..120]));
    let duplicate = format!("{reference}{reference}");
    let mut alt_duplicate = duplicate.clone();
    alt_duplicate.replace_range(100..200, &io::reverse_complement(&reference[100..200]));
    for (r, a) in [
        (&reference, &remote),
        (&reference, &short),
        (&duplicate, &alt_duplicate),
    ] {
        let (full, partial) = partition::build(&[parent(r, a, true)], r, 31, 31).unwrap();
        assert!(full.units[0][1].iter().all(|u| !u.preserve));
        assert!(partial.units[0][1].iter().all(|u| !u.preserve));
        replay(&full);
        replay(&partial);
    }
}

#[test]
fn sequence_anchor_setting_does_not_lower_structural_evidence_threshold() {
    let reference = dna(300, 4545);
    let mut alt = reference.clone();
    alt.replace_range(100..140, &io::reverse_complement(&reference[100..140]));
    let (plan, _) =
        partition::build(&[parent(&reference, &alt, true)], &reference, 19, 51).unwrap();
    assert!(plan.units[0][1].iter().all(|u| !u.preserve));
    replay(&plan);
}

#[test]
fn forward_matches_crossing_inversion_boundaries_keep_their_flanks() {
    for seed in 1..40 {
        let reference = dna(400, seed);
        let mut alt = reference.clone();
        alt.replace_range(100..200, &io::reverse_complement(&reference[100..200]));
        alt.replace_range(
            330..331,
            if &reference[330..331] == "A" {
                "C"
            } else {
                "A"
            },
        );
        // Maximal diagonal matches can extend by chance into either boundary.
        let mut p = parent(&reference, &alt, false);
        let mut i = 1;
        while i < reference.len() - 1 {
            if reference.as_bytes()[i] != alt.as_bytes()[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < reference.len() - 1 && reference.as_bytes()[i] == alt.as_bytes()[i] {
                i += 1;
            }
            p.anchors[1].push(Anchor {
                r: start,
                q: start,
                n: i - start,
                node: Some(1),
            });
        }
        let (full, partial) = partition::build(&[p], &reference, 31, 31).unwrap();
        for plan in [&full, &partial] {
            replay(plan);
            replay(&plan.reversed());
            let edits = alleles::decompose(plan).unwrap();
            assert!(
                edits[0][1].iter().any(|e| e.start == 330 && e.end == 331),
                "seed {seed}"
            );
            let u = plan.units[0][1].iter().find(|u| u.preserve).unwrap();
            assert!(u.r >= 99 && u.r + u.rn <= 201, "seed {seed}");
        }
    }
}
