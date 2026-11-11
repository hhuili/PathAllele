"""Read-only, example-specific audit of frozen caller outputs; not caller logic."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[4])
p.add_argument('--output', type=Path, default=Path(__file__).with_name('oriP_case_audit.json'))
a = p.parse_args()
root = a.project_root
sample = 'E5_sLCL'
start, end = 8403, 8507  # Published local illustration, 0-based half-open.
refpath = root/'variants/variants_PathAllele/input/reference.fa'
reference = ''.join(l.strip() for l in refpath.read_text().splitlines() if not l.startswith('>'))
gfapath = root/'graph_inversion/output/cpc_paths/graph.gfa'
nodes, walks = {}, {}
for line in gfapath.read_text().splitlines():
    f = line.split('\t')
    if f[0] == 'S': nodes[f[1]] = f[2]
    if f[0] == 'P': walks[f[1]] = f[2].split(',')
def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def spell(name):
    return ''.join(nodes[t[:-1]] if t[-1] == '+' else rc(nodes[t[:-1]]) for t in walks[name])
ref_window, expected = spell('Type1A'), spell('Type1B')
assert reference[start:end] == ref_window
assert len(expected) == len(ref_window) == 104
assert expected[:8] == ref_window[:8] and expected[-10:] == ref_window[-10:]
assert ref_window[8:10] == 'TT' and expected[8:10] == 'AC'
assert expected[10:94] == rc(ref_window[10:94])
# Cross-check the illustrative reconstruction against the actual PathAllele source.
sourcepath = root/'variants/variants_PathAllele/output/Graph.alleles.tsv.gz'
with gzip.open(sourcepath, 'rt') as f:
    source = next(d for d in csv.DictReader(f, delimiter='\t')
                  if d['block_id'] == 'PC009732' and d['allele_index'] == '1')
assert source['sequence'].count(expected) == 1
report = {'sample': sample, 'reference_interval_1based': [start+1, end],
          'selected_inversion_1based': [8414, 8497],
          'reference_sequence': ref_window, 'sample_sequence': expected,
          'display_source': 'validated-variant-informed reconstruction, not native PGGB topology',
          'source_allele_id': 'PC009732', 'source_crosscheck': True,
          'callers': {}, 'inputs': {}}
paths = {
    'PathAllele': 'variants/variants_PathAllele/output/Graph.vcf.gz',
    'VG': 'variants/vg/output/Graph.vcf.gz',
    'Linear': 'variants/linear/output/Linear.vcf.gz',
}
for name, rel in paths.items():
    file = root/rel
    edits = []
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                ix = line.rstrip().split('\t').index(sample)
                continue
            if line.startswith('#'): continue
            t = line.rstrip().split('\t')
            s, e = int(t[1])-1, int(t[1])-1+len(t[3])
            if s >= end: break
            if e <= start: continue
            fields = dict(zip(t[8].split(':'), t[ix].split(':')))
            gt = fields.get('GT', '.')
            # Missing stays missing; only carried ALT edits are replayed.
            if not gt.isdigit() or int(gt) == 0: continue
            assert s >= start and e <= end
            alt = t[4].split(',')[int(gt)-1]
            assert set(alt) <= set('ACGT') and reference[s:e] == t[3]
            info = dict(x.split('=', 1) for x in t[7].split(';') if '=' in x)
            selected = {}
            for key in ['TYPE', 'SEQTYPE', 'NATIVE', 'INVBASIS']:
                if key in info:
                    v = info[key].split(',')
                    selected[key] = v[int(gt)-1] if len(v) > 1 else v[0]
            edits.append({'start0': s, 'end0': e, 'ref': t[3], 'alt': alt,
                          'id': t[2], 'gt': int(gt), 'annotation': selected})
    edits.sort(key=lambda d: (d['start0'], d['end0']))
    replay, cursor = '', start
    for d in edits:
        assert d['start0'] >= cursor
        replay += reference[cursor:d['start0']] + d['alt']
        cursor = d['end0']
    replay += reference[cursor:end]
    assert replay == expected, name
    report['callers'][name] = {'n_carried_records': len(edits), 'replay_matches': True,
                               'reconstructed_length': len(replay), 'edits': edits,
                               'missing_policy': 'No missing genotype was recoded to REF; replay is of carried ALT records, not a genome-wide callability assertion.'}
    report['inputs'][rel] = hashlib.sha256(file.read_bytes()).hexdigest()
# Check raw vg representation before decomposition.
rawpath = root/'variants/vg/output/Graph.raw.vcf.gz'
with gzip.open(rawpath, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'): ix = line.rstrip().split('\t').index(sample); continue
        if line.startswith('#'): continue
        t = line.rstrip().split('\t')
        if t[1] != '6823' or t[2] != '>56518>65227': continue
        gt = int(t[ix].split(':')[0]); alt = t[4].split(',')[gt-1]
        assert expected in alt
        report['vg_raw'] = {'id': t[2], 'pos1': int(t[1]), 'ref_length': len(t[3]),
                            'selected_alt_length': len(alt), 'gt': gt,
                            'contains_exact_104bp_sample_sequence': True}
        break
assert 'vg_raw' in report
for file in [refpath, gfapath, sourcepath, rawpath]:
    report['inputs'][str(file.relative_to(root))] = hashlib.sha256(file.read_bytes()).hexdigest()
report['interpretation'] = 'All three frozen representations preserve this sample sequence. Structural classification, not sequence discovery alone, distinguishes the outputs. This case is not an independent accuracy benchmark.'
a.output.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({k: {'records':v['n_carried_records'], 'replay':v['replay_matches']}
                  for k,v in report['callers'].items()}))
