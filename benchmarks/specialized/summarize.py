#!/usr/bin/env python3
"""Summarize separate H+E measurements; paired block bootstrap, no timing."""
import argparse
import csv
import json
import random
import statistics
import sys
from pathlib import Path
from audit import validate_complete_records
from run import BINARY, HERE, ROOT, SOURCES, digest, validate_manifest


def interval(values, rng):
    blocks=[values[i:i+4] for i in range(0,len(values),4)]
    res=sorted(statistics.median(v for _ in blocks for v in rng.choice(blocks)) for _ in range(5000))
    return res[124],res[4874]


def validate_historical_sources(manifest):
    """Require frozen source bytes, never a reader's native build artifact."""
    frozen=manifest['source_binary_sha256']
    expected_paths={str(p.relative_to(ROOT)) for p in SOURCES+[BINARY]}
    if set(frozen)!=expected_paths:raise ValueError('Unexpected frozen source/binary path set')
    if any(len(value)!=64 or any(c not in '0123456789abcdef' for c in value) for value in frozen.values()):raise ValueError('Invalid source/binary SHA256')
    verified=[]
    for source in SOURCES:
        relative=str(source.relative_to(ROOT));expected=frozen[relative]
        historical=HERE/'frozen_sources'/expected/source.name
        for candidate in (source,historical):
            if candidate.is_file() and digest(candidate)==expected:
                verified.append(dict(manifest_path=relative,sha256=expected,verified_from=str(candidate.relative_to(ROOT))))
                break
        else:raise ValueError(f'Frozen source bytes unavailable or hash mismatch: {relative}')
    binary_hash=frozen[str(BINARY.relative_to(ROOT))]
    local_hash=digest(BINARY) if BINARY.is_file() else None
    return dict(historical_source_integrity='all manifest source hashes independently verified from available bytes',
                verified_sources=verified,historical_binary_sha256=binary_hash,
                historical_binary_integrity='hash consistency across frozen manifest, completed environment and every per-job receipt; historical execution is trusted from these receipts',
                local_binary_observation=('absent' if local_hash is None else 'matches frozen binary' if local_hash==binary_hash else 'differs from frozen binary; irrelevant to recorded-data analysis'),
                local_binary_sha256=local_hash)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('raw',type=Path)
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--environment',type=Path)
    ap.add_argument('--integrity-report',type=Path,help='write explicit historical source verification and binary receipt trust report')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    manifest=json.loads(args.manifest.read_text());validate_manifest(manifest)
    manifest_sha=digest(args.manifest)
    integrity=validate_historical_sources(manifest)
    env_path=args.environment or args.raw.with_name('specialized_environment.json')
    environment=json.loads(env_path.read_text())
    if environment['manifest']!=manifest or environment['manifest_sha256']!=manifest_sha:raise ValueError('Environment/manifest mismatch')
    if not environment.get('completed_utc'):raise ValueError('Run has no completion receipt')
    records=[json.loads(line) for line in args.raw.read_text().splitlines()]
    validate_complete_records(records,manifest,manifest_sha)
    integrity.update(manifest_sha256=manifest_sha,environment_manifest_matches=True,completed_environment_receipt=True,validated_job_receipts=len(records))
    if args.integrity_report:
        with args.integrity_report.open('x') as f:json.dump(integrity,f,indent=2);f.write('\n')
    print(f"Integrity: all frozen source hashes verified; {len(records)} complete matching job receipts; local binary {integrity['local_binary_observation']}. Historical execution/binary provenance relies on the recorded receipts.",file=sys.stderr)
    rows=[]
    for index,(r,job) in enumerate(zip(records,manifest['jobs'],strict=True)):
        samples=r['samples']
        row={k:r[k] for k in ('family','m','n','q_cap','scale_exp','seed','dist','output_full','actual_q','max_observed_q','actual_e')}
        for k in ('0','4','13','14'):row[f'method{k}_median_ms']=statistics.median(s[k] for s in samples)*1000
        for index,(numerator,denominator) in enumerate((('0','13'),('4','13'),('14','13'),('0','14'),('4','14'))):
            values=[s[numerator]/s[denominator] for s in samples]
            name=f'{numerator}_over_{denominator}'
            row[name]=statistics.median(values)
            row[name+'_low'],row[name+'_high']=interval(values,random.Random(2026101903+r['job_index']*7+index))
        rows.append(row)
    if not rows:raise SystemExit('No measurements')
    with args.output.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(args.output)

if __name__=='__main__':main()
