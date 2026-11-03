use pathallele::{classify::Classifier, io};
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
#[test]
fn exact_and_shifted_reverse_homology_have_different_contracts() {
    let reference = dna(600, 1209);
    let c = Classifier::new(&reference, 31).unwrap();
    let a = c
        .annotate(
            180,
            280,
            &io::reverse_complement(&reference[180..280]),
            "SEQUENCE",
        )
        .unwrap();
    assert_eq!(
        (&*a.kind, &*a.sequence, &*a.basis),
        ("INV", "MNV", "EXACT_RC")
    );
    let b = c
        .annotate(
            180,
            280,
            &io::reverse_complement(&reference[184..284]),
            "SEQUENCE",
        )
        .unwrap();
    assert_eq!(
        (&*b.kind, &*b.sequence, &*b.basis),
        ("INV_COMPLEX", "MNV", "LOCAL_RC")
    );
    assert!(b
        .segments
        .iter()
        .all(|s| reference[s.start..s.end] == s.template));
    assert!(b.segments.iter().any(|s| s.start >= 184 && s.end <= 284));
    // Matching padding is irrelevant to the class and template coordinates.
    let padded = format!(
        "{}{}{}",
        &reference[160..180],
        io::reverse_complement(&reference[180..280]),
        &reference[280..300]
    );
    let p = c.annotate(160, 300, &padded, "SEQUENCE").unwrap();
    assert_eq!(p.kind, "INV");
    assert_eq!(p.segments[0].start, a.segments[0].start);
    assert_eq!(p.segments[0].query, a.segments[0].query + 20);
}
#[test]
fn local_inverted_component_survives_junction_changes_and_substitution() {
    let reference = dna(800, 98331);
    let c = Classifier::new(&reference, 31).unwrap();
    let mut alt = io::reverse_complement(&reference[200..400]).into_bytes();
    alt[100] = if alt[100] == b'A' { b'C' } else { b'A' };
    let alt = String::from_utf8(alt).unwrap();
    let a = c.annotate(200, 400, &alt, "SEQUENCE").unwrap();
    assert_eq!(a.kind, "INV_COMPLEX");
    assert!(a.segments.len() >= 2);
    let alt = format!(
        "GATTACAGATTACA{}",
        io::reverse_complement(&reference[204..404])
    );
    let a = c.annotate(200, 400, &alt, "SEQUENCE").unwrap();
    assert_eq!((&*a.kind, &*a.sequence), ("INV_COMPLEX", "REPLACEMENT"));
}
#[test]
fn reverse_copy_insertions_remote_templates_and_repeats_are_not_pure_inversions() {
    let reference = dna(800, 1099);
    let c = Classifier::new(&reference, 31).unwrap();
    let remote = io::reverse_complement(&reference[600..700]);
    assert_eq!(
        c.annotate(200, 300, &remote, "SEQUENCE").unwrap().kind,
        "MNV"
    );
    assert_eq!(
        c.annotate(200, 200, &remote, "SEQUENCE").unwrap().kind,
        "INS"
    );
    let duplicated = format!("{reference}{}", &reference[200..300]);
    let d = Classifier::new(&duplicated, 31).unwrap();
    let alt = io::reverse_complement(&reference[200..300]);
    assert_eq!(d.annotate(200, 300, &alt, "SEQUENCE").unwrap().kind, "MNV");
    let tandem = alt.repeat(2);
    assert_eq!(
        c.annotate(200, 300, &tandem, "SEQUENCE").unwrap().kind,
        "REPLACEMENT"
    );
    assert_eq!(
        c.annotate(
            200,
            220,
            &io::reverse_complement(&reference[200..220]),
            "SEQUENCE"
        )
        .unwrap()
        .kind,
        "MNV"
    );
}
#[test]
fn ordinary_edits_and_native_structure_use_the_same_classifier() {
    let reference = dna(500, 771);
    let c = Classifier::new(&reference, 31).unwrap();
    let alt = if &reference[200..201] == "A" {
        "C"
    } else {
        "A"
    };
    assert_eq!(c.annotate(200, 201, alt, "SEQUENCE").unwrap().kind, "SNV");
    assert_eq!(c.annotate(200, 280, "", "SEQUENCE").unwrap().kind, "DEL");
    assert_eq!(
        c.annotate(200, 280, &dna(80, 23), "COMPLEX").unwrap().kind,
        "COMPLEX"
    );
    assert!(c
        .annotate(200, 280, &reference[200..280], "SEQUENCE")
        .is_err());
    assert!(c.annotate(200, 280, "NNNN", "SEQUENCE").is_err());
    // Palindromes and length changes alone cannot establish inversion.
    assert!(c.annotate(200, 280, &dna(80, 23), "BAD").is_err());
}

#[test]
fn classification_is_strand_invariant_and_respects_the_evidence_threshold() {
    for seed in 1..=40 {
        let reference = dna(1000, seed * 1231);
        let reverse = io::reverse_complement(&reference);
        let c = Classifier::new(&reference, 31).unwrap();
        let d = Classifier::new(&reverse, 31).unwrap();
        for shift in [0, 4] {
            let alt = io::reverse_complement(&reference[300 + shift..500 + shift]);
            let a = c.annotate(300, 500, &alt, "SEQUENCE").unwrap();
            let b = d
                .annotate(500, 700, &io::reverse_complement(&alt), "SEQUENCE")
                .unwrap();
            assert_eq!(a.kind, b.kind);
            // Shifted templates can occasionally give a smaller exact INV if
            // the junction bases coincide; classify the observed change.
            let r = &reference[300..500];
            let left = r
                .bytes()
                .zip(alt.bytes())
                .take_while(|(x, y)| x == y)
                .count();
            let right = r[left..]
                .bytes()
                .rev()
                .zip(alt[left..].bytes().rev())
                .take_while(|(x, y)| x == y)
                .count();
            let exact =
                alt[left..alt.len() - right] == io::reverse_complement(&r[left..r.len() - right]);
            assert_eq!(a.kind, if exact { "INV" } else { "INV_COMPLEX" });
        }
        let c = Classifier::new(&reference, 101).unwrap();
        let a = c
            .annotate(
                300,
                400,
                &io::reverse_complement(&reference[300..400]),
                "SEQUENCE",
            )
            .unwrap();
        assert_eq!(a.kind, "MNV");
    }
}

#[test]
fn palindrome_flanks_cannot_inflate_short_reverse_evidence_in_compound_loci() {
    for n in [17, 45] {
        let arm = dna(20, 539);
        let template = format!("{arm}{}{}", dna(n, 9871), io::reverse_complement(&arm));
        let reference = format!("{}{template}{}", dna(100, 397), dna(100, 887));
        let end = 100 + template.len() + 20;
        let mut alt = format!(
            "{}{}{}",
            &reference[80..100],
            io::reverse_complement(&template),
            &reference[100 + template.len()..end]
        )
        .into_bytes();
        let last = alt.len() - 1;
        for pos in [0, last] {
            alt[pos] = if alt[pos] == b'A' { b'C' } else { b'A' };
        }
        let alt = String::from_utf8(alt).unwrap();
        let result = Classifier::new(&reference, 31)
            .unwrap()
            .annotate(80, end, &alt, "SEQUENCE")
            .unwrap();
        assert_eq!(result.kind, if n < 31 { "MNV" } else { "INV_COMPLEX" });
    }
}
