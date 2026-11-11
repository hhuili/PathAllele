# Figure QA

2026-11-11. Source: `build_source_resolved.py`; official VectorCraft CLI 0.23.0,
safe-mode export. Output hashes are in `export_manifest.json`.

- One 180 × 58.5 mm page; 102 vector paths, no raster images.
- 34 text spans; embedded Helvetica regular/bold; minimum 5.53 pt.
- Collision audit: PASS, zero failures and warnings. Fourteen contained overlays are intentional node labels.
- Render inspected: compact three-part information flow, short synthetic bases, no actual oriP coordinates, sample names, 84-bp sequence or ellipses.
- Sequence check: reverse complement of AGTCCCGA is TCGGGACT; TT→AC adds a neighboring change. The 8-bp drawing is below the default 31-bp evidence threshold, not a full-caller fixture.
- Two forward graph branches do not imply a path-native inversion. Crossed links and ALT arrows in the complete-allele view show sequence relationships relative to REF.
- Occurrence evidence and mapping converge on source resolution; shared-node identity alone does not establish a coordinate or eligible anchor.
- Main art omits bounded graph-only supplementation and conditional WFA2; the legend and method describe them.
- The existing VectorCraft illustration route is retained. Matplotlib/R source-export checks do not apply; no claim is made that they passed.

Automated QA does not establish aesthetic approval, method novelty or biological accuracy.
