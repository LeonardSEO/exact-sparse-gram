#!/usr/bin/env python3
"""Prepare/freeze or run the separate H+E comparison. Never runs by default."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
from audit import validate_input, validate_record

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BINARY = HERE / 'build/specialized_benchmark'
SOURCES = [HERE/'src/benchmark.c', HERE/'src/hankel_sparse.c',
           HERE/'src/hankel_sparse.h', HERE/'tests/test_specialized.c',
           HERE/'Makefile', HERE/'run.py', HERE/'summarize.py', HERE/'audit.py',
           HERE/'tests/test_harness.py', HERE/'tests/test_receipts.py', ROOT/'src/lattice_sequential.c']
METHODS = {'0':'Accelerate SYRK', '4':'Accelerate GEMM',
           '13':'preselected sequential lattice', '14':'reconstructed Hankel plus sparse corrections'}

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def capture(*cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return dict(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)

def planned():
    jobs = [dict(family=f, m=512, n=4096, q=2, e=-15, seed=s, dist=d, full=1)
            for f in ('windows', 'innov') for d in (0, 1)
            for s in range(20261201, 20261206)]
    random.Random(2026101903).shuffle(jobs)
    return dict(phase='specialized_comparator_separate_from_original_230',
                jobs=jobs, method_names=METHODS, samples=32, blocks=8,
                balanced_block_size=4, batch_ms=10, warmups=5,
                candidate_preselected='13',
                distributions={'0':'uniform signed16', '1':'nearconstant integer'},
                output='full; same blocked32 mirror for triangular methods',
                input='same dense generator and seeds as original Round2; no latent H/edit metadata',
                timing='complete calls, allocations, validation, H reconstruction, E detection/fill and reflection included',
                accuracy='poison+allfinite; independent int64 sample oracle; exact equality for 13 and 14; norm-scaled 2e-12 for vendors',
                flags='-O3 -std=c11 -fno-fast-math -ffp-contract=off -Wall -Wextra -Werror')

def hashes():
    return {str(p.relative_to(ROOT)):digest(p) for p in SOURCES+[BINARY]}

def command(job,samples=32,batch_ms=10,diagnostic=False):
    result=[str(BINARY),job['family'],str(job['m']),str(job['n']),str(job['q']),str(job['e']),
            str(job['seed']),str(job['dist']),str(samples),str(job['full']),str(batch_ms)]
    if diagnostic:result.append('--input-diagnostic')
    return result

def validate_manifest(manifest):
    for key,value in planned().items():
        if manifest.get(key)!=value:raise ValueError(f'Unexpected frozen parameter {key}')
    if len(manifest['input_diagnostics'])!=20:raise ValueError('Expected 20 frozen input diagnostics')
    fingerprints={}
    for job,record in zip(manifest['jobs'],manifest['input_diagnostics'],strict=True):
        validate_input(record,job,manifest['samples'])
        group=(job['family'],job['dist'])
        fingerprints.setdefault(group,set()).add(record['input_fingerprint_fnv1a64'])
    if any(len(values)!=5 for values in fingerprints.values()):raise ValueError('Seeds did not generate five distinct inputs in each family/distribution')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    action=ap.add_mutually_exclusive_group()
    action.add_argument('--freeze', type=Path, help='write immutable manifest; no timing')
    action.add_argument('--run', type=Path, help='run an existing frozen manifest')
    ap.add_argument('--output-dir', type=Path, default=HERE/'results')
    args=ap.parse_args()
    plan=planned()
    if args.freeze:
        plan.update(source_binary_sha256=hashes(), created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                    input_diagnostics=[json.loads(subprocess.check_output(command(job,diagnostic=True),text=True)) for job in plan['jobs']])
        validate_manifest(plan)
        if hashes()!=plan['source_binary_sha256']:raise RuntimeError('Source changed during untimed input diagnostics')
        with args.freeze.open('x') as f:json.dump(plan,f,indent=2);f.write('\n')
        print(args.freeze)
        return
    if not args.run:
        print(json.dumps(plan,indent=2));return
    manifest=json.loads(args.run.read_text())
    expected=manifest['source_binary_sha256']
    validate_manifest(manifest)
    manifest_sha=digest(args.run)
    def check_frozen():
        if hashes()!=expected:raise RuntimeError('Frozen source or binary changed')
    check_frozen()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    raw_path=args.output_dir/'specialized_raw.jsonl'
    env_path=args.output_dir/'specialized_environment.json'
    if raw_path.exists() or env_path.exists():raise SystemExit('Refusing to overwrite measurement receipts')
    env=os.environ.copy();env['VECLIB_MAXIMUM_THREADS']='1'
    receipt=dict(manifest=manifest, manifest_sha256=manifest_sha,
                 started_utc=dt.datetime.now(dt.timezone.utc).isoformat(), platform=platform.platform(),
                 compiler=capture('clang','--version'), os=capture('sw_vers'),
                 hardware=capture('sysctl','-n','machdep.cpu.brand_string','hw.ncpu','hw.memsize'),
                 power=capture('pmset','-g','batt'), libraries=capture('otool','-L',str(BINARY)),
                 threading='BLASSetThreading single, BLASGetThreading checked, VECLIB_MAXIMUM_THREADS=1',
                 scheduling='standard macOS scheduler, no affinity; operator ensures quiet serial timing')
    env_path.write_text(json.dumps(receipt,indent=2)+'\n')
    with raw_path.open('x') as out:
        for index, job in enumerate(manifest['jobs']):
            check_frozen()
            cmd=command(job,manifest['samples'],manifest['batch_ms'])
            r=subprocess.run(cmd,env=env,capture_output=True,text=True,check=True)
            check_frozen()
            if digest(args.run)!=manifest_sha:raise RuntimeError('Manifest changed during run')
            raw=json.loads(r.stdout);raw.update(job_index=index,phase=manifest['phase'],checks=r.stderr,
                                               manifest_sha256=manifest_sha,source_binary_sha256=expected)
            validate_record(raw,job,manifest,index,manifest_sha)
            out.write(json.dumps(raw)+'\n');out.flush()
            ratio={k:statistics.median(s[k]/s['13'] for s in raw['samples']) for k in ('0','4','14')}
            print(f'{index+1}/{len(manifest["jobs"])} {job} ratios_to_seq={ratio} nnzE={raw["actual_e"]}',flush=True)
    check_frozen()
    receipt.update(completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                   final_power=capture('pmset','-g','batt'),thermal=capture('pmset','-g','therm'))
    env_path.write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
