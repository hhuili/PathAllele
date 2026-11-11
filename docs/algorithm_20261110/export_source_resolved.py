"""Rebuild the editable SVG, native VectorCraft document and publication PDF."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--cli', default=shutil.which('vectorcraft-cli'))
args = p.parse_args()
if not args.cli:
    p.error('Supply --cli /path/to/vectorcraft-cli (official VectorCraft release).')
cli = str(Path(args.cli).resolve())
subprocess.run([sys.executable, str(ROOT/'build_source_resolved.py')], check=True)
cmd = [cli, '--safe-mode', 'run', '--in', str(ROOT/'PathAllele_panel_a.svg'),
       '--export', str(ROOT/'PathAllele_panel_a.vectorcraft'),
       '--export', str(ROOT/'PathAllele_panel_a.pdf')]
result = subprocess.run(cmd, check=True, text=True, capture_output=True)
steps = [json.loads(x) for x in result.stdout.splitlines() if x.strip()]
for step in steps:
    if step.get('result', {}).get('warnings'):
        raise RuntimeError(f'Export warning: {step}')
for step in steps:
    if 'path' in step: step['path'] = Path(step['path']).name
    if 'path' in step.get('result', {}): step['result']['path'] = Path(step['result']['path']).name
files = {}
for ext in ['svg','vectorcraft','pdf']:
    f=ROOT/f'PathAllele_panel_a.{ext}'
    files[f.name]={'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
manifest={'backend':'Official VectorCraft CLI','tested_version':'0.23.0',
          'executable':Path(cli).name,'executable_sha256':hashlib.sha256(Path(cli).read_bytes()).hexdigest(),
          'command':[Path(v).name if v.startswith('/') else v for v in cmd],'steps':steps,'files':files,'figure_mm':[180,58.5],
          'native_master':'PathAllele_panel_a.vectorcraft','publication':'PathAllele_panel_a.pdf'}
(ROOT/'export_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(result.stdout,end='')
