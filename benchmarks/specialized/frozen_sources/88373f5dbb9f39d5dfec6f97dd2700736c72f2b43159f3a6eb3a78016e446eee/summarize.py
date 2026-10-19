#!/usr/bin/env python3
"""Summarize separate H+E measurements; paired block bootstrap, no timing."""
import argparse
import csv
import json
import random
import statistics
from pathlib import Path
from audit import validate_complete_records
from run import digest, hashes, validate_manifest


def interval(values, rng):
    blocks=[values[i:i+4] for i in range(0,len(values),4)]
    res=sorted(statistics.median(v for _ in blocks for v in rng.choice(blocks)) for _ in range(5000))
    return res[124],res[4874]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('raw',type=Path)
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--environment',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    manifest=json.loads(args.manifest.read_text());validate_manifest(manifest)
    manifest_sha=digest(args.manifest)
    if hashes()!=manifest['source_binary_sha256']:raise ValueError('Current sources/binary differ from frozen hashes')
    env_path=args.environment or args.raw.with_name('specialized_environment.json')
    environment=json.loads(env_path.read_text())
    if environment['manifest']!=manifest or environment['manifest_sha256']!=manifest_sha:raise ValueError('Environment/manifest mismatch')
    if not environment.get('completed_utc'):raise ValueError('Run has no completion receipt')
    records=[json.loads(line) for line in args.raw.read_text().splitlines()]
    validate_complete_records(records,manifest,manifest_sha)
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
