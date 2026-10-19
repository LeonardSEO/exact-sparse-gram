#!/usr/bin/env python3
"""Prespecified method13; paired medians and six-trial block bootstraps."""
import csv
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

ROOT = Path(__file__).resolve().parent

def interval(values):
    assert len(values) == 30
    blocks = [values[k:k+6] for k in range(0, 30, 6)]
    rng = random.Random(2026101903)
    draws = sorted(statistics.median(x for block in rng.choices(blocks, k=5) for x in block)
                   for _ in range(5000))
    return statistics.median(values), draws[124], draws[4874]

def summarize(d):
    samples = d['samples']
    med = {k: statistics.median(s[k] for s in samples) for k in samples[0]}
    sr, lo, hi = interval([s['0']/s['13'] for s in samples])
    gr, glo, ghi = interval([s['4']/s['13'] for s in samples])
    # Baseline chosen once by its per-run median, not a per-sample oracle.
    baseline = min(('0', '4'), key=med.get)
    br, blo, bhi = interval([s[baseline]/s['13'] for s in samples])
    r = {k: d[k] for k in ('evaluation_set', 'family', 'm', 'n', 'q_cap',
                           'scale_exp', 'seed', 'dist', 'output_full')}
    r.update({f'method_{k}_ms': v*1000 for k, v in med.items()})
    r.update(speedup_syrk=sr, syrk_ci_lo=lo, syrk_ci_hi=hi,
             speedup_gemm=gr, gemm_ci_lo=glo, gemm_ci_hi=ghi,
             best_vendor=baseline, speedup_best_vendor=br,
             best_ci_lo=blo, best_ci_hi=bhi)
    return r

def write_selected(rows):
    examples = [('innov', 768, 8192, 2, -15), ('innov', 768, 8192, 4, 0),
                ('windows', 512, 4096, 2, -15), ('innov', 512, 4096, 2, -15),
                ('innov', 384, 2048, 4, -15), ('innov', 192, 768, 2, -15)]
    lines = ['# Frozen Round2 comparisons', '',
             'All example rows: uniform signed16-bit coefficients, seed20261201, full output.',
             'Times are per-method medians. Speedups are paired-ratio medians;',
             'the 95% intervals resample five contiguous six-trial balanced blocks.',
             'Intervals describe within-run variability, not cross-hardware coverage.', '',
             '| Family | Shape | q cap | e | SYRK ms | GEMM ms | Earlier certified ms | Sequential ms | SYRK/seq (95% CI) | GEMM/seq (95% CI) |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for family, m, n, q, e in examples:
        r = next(r for r in rows if (r['family'], r['m'], r['n'], r['q_cap'],
                 r['scale_exp'], r['seed'], r['dist'], r['output_full']) ==
                 (family, m, n, q, e, 20261201, 0, 1))
        lines.append(f"|{family}|{m}×{n}|{q}|{e}|{r['method_0_ms']:.4f}|"
                     f"{r['method_4_ms']:.4f}|{r['method_9_ms']:.4f}|{r['method_13_ms']:.4f}|"
                     f"{r['speedup_syrk']:.3f} [{r['syrk_ci_lo']:.3f}, {r['syrk_ci_hi']:.3f}]|"
                     f"{r['speedup_gemm']:.3f} [{r['gemm_ci_lo']:.3f}, {r['gemm_ci_hi']:.3f}]|")
    groups = defaultdict(list)
    for r in rows:
        key = (r['evaluation_set'], r['family'], r['m'], r['n'], r['q_cap'],
               r['dist'] if r['evaluation_set'] == 'distribution' else -1,
               r['output_full'])
        groups[key].append(r)
    lines += ['', '## Every group, including losses', '',
              'Heldout groups combine two scales and five seeds (10 runs); distribution groups have five runs.',
              'These are ranges and medians of per-run paired speedups, not pooled confidence intervals.',
              'Lower/full, scaled versions and deterministic patterns are repeated inputs, not independent draws.', '',
              '| Set | Family | Shape | q | dist | output | runs | SYRK/seq min, median, max | Best-vendor/seq min, median, max |',
              '|---|---|---:|---:|---:|---|---:|---:|---:|']
    stats = []
    for key, rr in sorted(groups.items()):
        eval_set, family, m, n, q, dist, full = key
        s = [r['speedup_syrk'] for r in rr]
        b = [r['speedup_best_vendor'] for r in rr]
        v = [min(s), statistics.median(s), max(s)]
        w = [min(b), statistics.median(b), max(b)]
        stats.append(dict(evaluation_set=eval_set, family=family, m=m, n=n,
                          q_cap=q, dist=dist, output_full=full, runs=len(rr),
                          syrk_min_median_max=v, best_vendor_min_median_max=w))
        lines.append(f"|{eval_set}|{family}|{m}×{n}|{q}|{dist if dist>=0 else '0'}|"
                     f"{'full' if full else 'lower'}|{len(rr)}|"
                     + ', '.join(f'{x:.3f}' for x in v) + '|'
                     + ', '.join(f'{x:.3f}' for x in w) + '|')
    wins = sum(r['speedup_syrk'] > 1 and r['speedup_gemm'] > 1 for r in rows)
    ci_wins = sum(r['syrk_ci_lo'] > 1 and r['gemm_ci_lo'] > 1 for r in rows)
    lines += ['', f'{len(rows)} manifest-checked runs; {wins} paired-median wins against both vendors.',
              f'{ci_wins} runs also have both conditional interval lower endpoints above1.',
              'All methods11/12/13, individual times, vendor choices and intervals remain in lattice_summary.csv.',
              'No per-case candidate selection is used: method13 is fixed throughout.', '']
    (ROOT/'results'/'lattice_selected_measurements.md').write_text('\n'.join(lines))
    (ROOT/'results'/'lattice_statistics.json').write_text(json.dumps(
        dict(runs=len(rows), paired_median_wins_both=wins,
             conditional_interval_wins_both=ci_wins, groups=stats,
             summarizer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')

def main():
    rows = []
    seen = set()
    for p in sorted((ROOT/'results').glob('lattice_final_*.jsonl')):
        raw_rows = [json.loads(line) for line in p.read_text().splitlines()]
        if not raw_rows:
            raise SystemExit(f'Empty comparison file: {p}')
        manifest_name = raw_rows[0]['source_manifest']
        manifest = json.loads((ROOT/'results'/manifest_name).read_text())
        if not manifest.get('completed_utc') or len(raw_rows) != len(manifest['jobs']):
            raise SystemExit(f'Incomplete evaluation: {p}')
        indices = set()
        for d in raw_rows:
            assert d['source_manifest'] == manifest_name
            index = d['job_index']
            assert index not in indices
            indices.add(index)
            job = manifest['jobs'][index]
            got = (d['family'], d['m'], d['n'], d['q_cap'], d['scale_exp'],
                   d['seed'], d['dist'], d['output_full'], d['evaluation_set'])
            expected = tuple(job[k] for k in ('family', 'm', 'n', 'q', 'e',
                                               'seed', 'dist', 'full', 'set'))
            assert got == expected and got not in seen
            seen.add(got)
            assert len(d['samples']) == 30 and d['candidate_preselected'] == '13'
            for sample in d['samples']:
                assert set(sample) == {'0', '4', '9', '11', '12', '13'}
                assert all(math.isfinite(v) and v > 0 for v in sample.values())
            rows.append(summarize(d))
        assert indices == set(range(len(manifest['jobs'])))
    if not rows:
        raise SystemExit('No round2 final receipts')
    rows.sort(key=lambda r: tuple(r[k] for k in ('evaluation_set', 'family', 'm',
              'n', 'q_cap', 'scale_exp', 'dist', 'seed', 'output_full')))
    out = ROOT/'results'/'lattice_summary.csv'
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    write_selected(rows)
    print(out, len(rows), 'runs; intervals conditional within-run, not cross-hardware guarantees')

if __name__ == '__main__':
    main()
