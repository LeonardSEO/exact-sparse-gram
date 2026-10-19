#!/usr/bin/env python3
"""Recompute archived CSV from source/data-only copies, no executable or timing."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run import HERE,ROOT,SOURCES,BINARY


def main():
    original_data=ROOT/'data/specialized_2026-10-19'
    with tempfile.TemporaryDirectory(prefix='specialized-portable-') as td:
        root=Path(td)/'release'
        for source in SOURCES:
            destination=root/source.relative_to(ROOT)
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,destination)
        shutil.copytree(HERE/'frozen_sources',root/HERE.relative_to(ROOT)/'frozen_sources')
        data=root/'data/specialized_2026-10-19';shutil.copytree(original_data,data)
        local_binary=root/BINARY.relative_to(ROOT)
        assert not local_binary.exists()
        for mode in ('absent','different'):
            if mode=='different':
                local_binary.parent.mkdir(parents=True,exist_ok=True)
                local_binary.write_text('Deliberately different local build fixture; never executed.\n')
            output=root/f'{mode}.csv';report=root/f'{mode}.json'
            subprocess.run([sys.executable,str(root/HERE.relative_to(ROOT)/'summarize.py'),
                            str(data/'specialized_raw.jsonl'),'--manifest',str(data/'frozen.json'),
                            '--output',str(output),'--integrity-report',str(report)],check=True)
            assert output.read_bytes()==(original_data/'summary.csv').read_bytes()
            integrity=json.loads(report.read_text())
            assert integrity['validated_job_receipts']==20
            if mode=='absent':assert integrity['local_binary_observation']=='absent'
            else:assert integrity['local_binary_observation'].startswith('differs from frozen binary')
            assert len(integrity['verified_sources'])==len(SOURCES)
        print('portable recorded analysis: no-binary and different-binary source/data-only copies reproduce all 20 CSV rows byte for byte')

if __name__=='__main__':main()
