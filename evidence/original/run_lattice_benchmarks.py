#!/usr/bin/env python3
"""Frozen round-2 evaluation. Standard library; run without competing compute."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess

ROOT = Path(__file__).resolve().parent
METHODS = {'0': 'Accelerate SYRK', '4': 'Accelerate GEMM',
           '9': 'certified floating', '11': 'grouped lattice',
           '12': 'merged substitution', '13': 'sequential lattice'}
SOURCES = [ROOT/'src'/f for f in (
    'lattice_benchmark.c', 'benchmark.c', 'lattice_exact.c',
    'lattice_substitute.c', 'lattice_sequential.c',
    'communication_guard_fast.c', 'numerical_sparse_fast.c')]

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def capture(*cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    return {'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}

def jobs_for(grid):
    jobs = []
    if grid in ('heldout', 'all'):
        jobs += [dict(family='innov', m=m, n=n, q=q, e=e, seed=seed,
                      dist=0, full=full, set='heldout')
                 for m, n in ((192, 768), (384, 2048), (768, 8192))
                 for q in (1, 2, 4) for e in (0, -15)
                 for seed in range(20261201, 20261206) for full in (0, 1)]
    if grid in ('distribution', 'all'):
        jobs += [dict(family=family, m=512, n=4096, q=2, e=-15,
                      seed=seed, dist=dist, full=1, set='distribution')
                 for family in ('windows', 'innov') for dist in range(5)
                 for seed in range(20261201, 20261206)]
    random.Random(2026101902).shuffle(jobs)
    return jobs

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--grid', choices=('heldout', 'distribution', 'all'), default='all')
    p.add_argument('--dry-run', action='store_true')
    args = p.parse_args()
    jobs = jobs_for(args.grid)
    manifest = {'jobs': jobs, 'samples': 30, 'batch_ms': 10,
                'method_names': METHODS, 'candidate_preselected': '13',
                'distribution_names': ['uniform signed16', 'nearconstant integer',
                                       'alternating saturated', 'impulse saturated',
                                       'ternary'],
                'deterministic_inputs': 'windows dist2/3 repeat across seeds; innovation edit locations remain seeded',
                'source_sha256': {str(s.relative_to(ROOT)): digest(s) for s in SOURCES},
                'binary_sha256': digest(ROOT/'lattice_benchmark'),
                'runner_sha256': digest(Path(__file__)),
                'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    if args.dry_run:
        print(json.dumps(manifest, indent=2))
        return
    out = ROOT/'results'/f'lattice_final_{args.grid}.jsonl'
    meta = ROOT/'results'/f'lattice_environment_{args.grid}.json'
    if out.exists() or meta.exists():
        raise SystemExit(f'Refusing to overwrite {out} or {meta}; preserve prior receipts.')
    env = os.environ.copy()
    env['VECLIB_MAXIMUM_THREADS'] = '1'
    manifest.update(platform=platform.platform(), python=platform.python_version(),
                    compiler=capture('clang', '--version'),
                    os=capture('sw_vers'), hardware=capture('sysctl', '-n',
                    'machdep.cpu.brand_string', 'hw.ncpu', 'hw.memsize', 'hw.cachelinesize'),
                    power=capture('pmset', '-g', 'batt'),
                    linked_libraries=capture('otool', '-L', str(ROOT/'lattice_benchmark')),
                    flags='-O3 -std=c11 -fno-fast-math -ffp-contract=off -Wall -Wextra',
                    threads='BLASSetThreading single; BLASGetThreading checked; VECLIB_MAXIMUM_THREADS=1',
                    scheduling='standard macOS scheduler, no affinity',
                    cache='warm repeated operands; per-call validation, allocation and conversion included',
                    accuracy='fixed norm-scaled2e-12; exact lattice output theorem; independent int64 preflight',
                    output='lower or full; blocked32mirror counted for triangular methods, native GEMM already full')
    meta.write_text(json.dumps(manifest, indent=2)+'\n')
    with out.open('x') as f:
        for index, job in enumerate(jobs):
            for source in SOURCES:
                if digest(source) != manifest['source_sha256'][str(source.relative_to(ROOT))]:
                    raise RuntimeError(f'Frozen source changed: {source}')
            if digest(ROOT/'lattice_benchmark') != manifest['binary_sha256']:
                raise RuntimeError('Frozen binary changed')
            cmd = [str(ROOT/'lattice_benchmark'), job['family'], str(job['m']),
                   str(job['n']), str(job['q']), str(job['e']), str(job['seed']),
                   str(job['dist']), '30', str(job['full']), '10']
            r = subprocess.run(cmd, env=env, text=True, capture_output=True, check=True)
            raw = json.loads(r.stdout)
            raw.update(phase='round2_final', evaluation_set=job['set'],
                       job_index=index, checks=r.stderr.strip(),
                       source_manifest=meta.name, candidate_preselected='13')
            f.write(json.dumps(raw)+'\n')
            f.flush()
            med = {k: statistics.median(s[k] for s in raw['samples']) for k in METHODS}
            ratios = [s['0']/s['13'] for s in raw['samples']]
            print(f"{index+1}/{len(jobs)} {job} SYRK_ms={med['0']*1e3:.4f} "
                  f"sequential_ms={med['13']*1e3:.4f} paired={statistics.median(ratios):.3f}", flush=True)
    manifest.update(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    final_power=capture('pmset', '-g', 'batt'),
                    thermal_status=capture('pmset', '-g', 'therm'))
    meta.write_text(json.dumps(manifest, indent=2)+'\n')
    print(out)

if __name__ == '__main__':
    main()
