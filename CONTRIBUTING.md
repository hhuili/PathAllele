# Contributing

Report a minimal synthetic graph, reference and sample list, the command, version, dependency versions and relevant error log. Avoid uploading private cohort data. Describe expected sequence behavior rather than an expected catalogue count alone.

A calling change must have a general rule and a failing synthetic test that becomes passing. No production path may consume benchmark target sites. Preserve unknown states; do not use missing-to-reference imputation. Check reference orientation, repeated occurrences, overlapping edits and compound representation before proposing a change. Run the full integration test suite and validate a new output directory before changing a scientific release.

Use Python 3.10+ and format Python with Ruff. Pull requests should explain the sequence or ownership problem, the new behavior, validation performed and remaining limitations. Original contributions use the project's MIT license.
