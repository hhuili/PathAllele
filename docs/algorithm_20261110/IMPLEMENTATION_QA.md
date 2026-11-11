# Checks for the documentation revision

2026-11-11; production source baseline `057d45a9c51aba959ecd88eece279a047a02681d`.
This revision changes documentation and artwork only. Caller code, parameters and
production outputs are unchanged. Historical adapter/build reports are archived
locally and are not current deployment claims.

- Reviewed graph/source/context/anchors, complete-allele projection and the general sequence classifier against the source.
- Local and deployed graph/source/context/anchors source hashes matched at inspection.
- Existing `classification` integration tests: 6 passed. Run from `/tmp` with a temporary official-registry `CARGO_HOME` to bypass an unavailable inherited mirror; no permanent Cargo configuration changed.
- The test suite covers exact and boundary-shifted inversions, additional changes, remote reversed insertion, repeat and palindrome controls using generated sequences. It does not establish biological accuracy.
- Read-only oriP audit: PathAllele (1 carried record), current vg/wave (2), and Linear (21) reconstruct the same 104-bp ALT. The raw vg snarl already contains that ALT. See `oriP_case_audit.json` and its script.
- Figure bases independently checked for reverse-complement identity. They are a short explanatory simulation, not a default-parameter caller test.
- Final PDF text, geometry and collision reports accompany the figure.

Reproduce the focused code check:

```sh
cargo test --locked --test classification
```

Reproduce the existing-data case audit from a checkout of the EBV analysis project:

```sh
python3 audit_orip_example.py --project-root /path/to/04_pangenome
```

The case audit needs the existing local reference, graph and VCF outputs; those
cohort data are not bundled with these documentation files.
