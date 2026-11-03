use pathallele::{
    model::Edit,
    structure::{cluster, compatible, resolve, within_distance, Allele},
};
fn allele(
    start: usize,
    kind: &str,
    reference: &str,
    alt: &str,
    reverse: Option<(usize, usize)>,
) -> Allele {
    Allele {
        edit: Edit {
            start,
            end: start + reference.len(),
            alt: alt.into(),
            complex: false,
        },
        reference: reference.into(),
        kind: kind.into(),
        reverse,
    }
}
fn distance(a: &[u8], b: &[u8]) -> usize {
    let mut d: Vec<_> = (0..=b.len()).collect();
    for (i, x) in a.iter().enumerate() {
        let mut left = i + 1;
        let mut diag = d[0];
        d[0] = left;
        for (j, y) in b.iter().enumerate() {
            let old = d[j + 1];
            left = (old + 1).min(left + 1).min(diag + usize::from(x != y));
            diag = old;
            d[j + 1] = left;
        }
    }
    d[b.len()]
}
#[test]
fn wavefront_distance_agrees_with_full_dynamic_programming() {
    let strings: Vec<Vec<u8>> = (0usize..128)
        .map(|x| {
            let n = usize::BITS as usize - x.leading_zeros() as usize;
            (0..n)
                .map(|i| if x & (1 << i) == 0 { b'A' } else { b'C' })
                .collect()
        })
        .collect();
    for a in &strings {
        for b in &strings {
            for k in 0..=7 {
                assert_eq!(
                    within_distance(a, b, k),
                    distance(a, b) <= k,
                    "{a:?} {b:?} {k}"
                );
            }
        }
    }
}
#[test]
fn inversion_subtypes_merge_but_distinct_axes_or_classes_do_not() {
    let r = "A".repeat(84);
    let a = "C".repeat(84);
    let x = allele(100, "INV_COMPLEX", &r, &a, Some((104, 188)));
    let mut b = a.clone();
    b.replace_range(50..51, "T");
    let y = allele(100, "INV_COMPLEX", &r, &b, Some((104, 188)));
    let z = allele(100, "INV_COMPLEX", &r, &(a.clone() + "G"), Some((103, 188)));
    assert!(compatible(&x, &y));
    assert!(compatible(&x, &z));
    assert!(compatible(&y, &z));
    let mut bad = y.clone();
    bad.reverse = Some((200, 284));
    assert!(!compatible(&x, &bad));
    bad = y.clone();
    bad.kind = "MNV".into();
    assert!(!compatible(&x, &bad));
    bad = y.clone();
    bad.reverse = None;
    assert!(!compatible(&x, &bad));
    bad = y;
    bad.edit.start += 100;
    bad.edit.end += 100;
    assert!(!compatible(&x, &bad));
}
#[test]
fn complete_linkage_blocks_transitive_bridges_and_is_order_invariant() {
    let r = "A".repeat(100);
    let a = "C".repeat(100);
    let b = "T".repeat(5) + &"C".repeat(95);
    let c = "T".repeat(10) + &"C".repeat(90);
    let items: Vec<_> = [a, b, c]
        .iter()
        .map(|s| allele(100, "REPLACEMENT", &r, s, None))
        .collect();
    assert!(compatible(&items[0], &items[1]));
    assert!(compatible(&items[1], &items[2]));
    assert!(!compatible(&items[0], &items[2]));
    let signature = |x: &[Allele]| {
        let mut r: Vec<_> = cluster(x)
            .iter()
            .map(|g| {
                let mut s: Vec<_> = g.iter().map(|&i| x[i].clone()).collect();
                s.sort();
                s
            })
            .collect();
        r.sort();
        r
    };
    let mut reversed = items.clone();
    reversed.reverse();
    assert_eq!(signature(&items), signature(&reversed));
    assert_eq!(cluster(&items).len(), 2);
}
#[test]
fn deletions_and_insertions_have_appropriate_sequence_comparisons() {
    let x = allele(100, "DEL", &"A".repeat(100), "", None);
    let y = allele(100, "DEL", &"A".repeat(99), "", None);
    assert!(compatible(&x, &y));
    let z = allele(100, "DEL", &"C".repeat(99), "", None);
    assert!(!compatible(&x, &z));
    let x = allele(100, "INS", "", &"A".repeat(100), None);
    let y = allele(100, "INS", "", &("A".repeat(99) + "C"), None);
    assert!(compatible(&x, &y));
    let mut shifted = y;
    shifted.edit.start += 1;
    shifted.edit.end += 1;
    assert!(!compatible(&x, &shifted));
}
#[test]
fn complete_locus_evidence_preserves_missingness_and_conflicts() {
    assert_eq!(resolve(false, &[0], true), (0, "PASS"));
    assert_eq!(resolve(false, &[-1], true), (-1, "NO_DATA"));
    assert_eq!(resolve(false, &[0], false), (-1, "NO_DATA"));
    assert_eq!(resolve(true, &[-1], true), (1, "PASS"));
    assert_eq!(resolve(true, &[0, 1], true), (-1, "AMBIGUOUS_ALIGNMENT"));
    assert_eq!(resolve(false, &[0, -1], true), (-1, "NO_DATA"));
    assert_eq!(resolve(false, &[1], true), (1, "PASS"));
}

#[test]
fn shared_reference_background_does_not_define_a_structural_group() {
    let r = "A".repeat(100);
    let mut a = r.clone();
    a.replace_range(0..1, "C");
    a.replace_range(99..100, "C");
    let mut b = a.clone();
    b.replace_range(50..51, "G");
    assert!(!compatible(
        &allele(100, "MNV", &r, &a, None),
        &allele(100, "MNV", &r, &b, None)
    ));
}
