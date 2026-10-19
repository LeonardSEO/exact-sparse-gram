#!/usr/bin/env python3
"""Generate manuscript tables and vector picture from retained statistics."""
import csv
import statistics
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((ROOT/'analysis/generated/lattice_summary.csv').open()))
assert len(rows) == 230
out = ROOT/'figures'
out.mkdir(exist_ok=True)
keys = [('innov',768,8192,2,-15),('innov',768,8192,4,0),
        ('windows',512,4096,2,-15),('innov',512,4096,2,-15),
        ('innov',384,2048,4,-15),('innov',192,768,2,-15)]
lines=[]
for family,m,n,q,e in keys:
 r=next(r for r in rows if tuple(r[k] for k in ('family','m','n','q_cap','scale_exp','seed','dist','output_full'))
        ==(family,str(m),str(n),str(q),str(e),'20261201','0','1'))
 label='Innovations' if family=='innov' else 'Windows'
 lines.append(f"{label} & ${m}\\times{n}$ & {q} & {e} & "
              + ' & '.join(f"{float(r[f'method_{k}_ms']):.4f}" for k in ('0','4','13'))
              + f" & {float(r['speedup_syrk']):.3f} [{float(r['syrk_ci_lo']):.3f}, {float(r['syrk_ci_hi']):.3f}] " + chr(92)*2)
(out/'selected_results.tex').write_text('\n'.join(lines)+'\n'+chr(92)+'bottomrule\n')
p=[r'\setlength{\unitlength}{1cm}',r'\begin{picture}(15,5.15)',r'\small']
for panel,(m,n) in enumerate(((192,768),(384,2048),(768,8192))):
 base=panel*5+.7
 p += [f'\\put({base:.2f},.75){{\\line(1,0){{3.65}}}}',
       f'\\put({base:.2f},.75){{\\line(0,1){{3.5}}}}',
       f'\\put({base+1.8:.2f},4.65){{\\makebox(0,0){{$ {m}\\times{n} $}}}}']
 for y in range(4):
  yy=.75+y*1.08
  p.append(f'\\put({base-.08:.2f},{yy:.3f}){{\\line(1,0){{.08}}}}')
  if panel==0:p.append(f'\\put({base-.24:.2f},{yy:.3f}){{\\makebox(0,0){{{y}}}}}')
 p.append(f'\\put({base:.2f},1.83){{\\color{{gray}}\\line(1,0){{3.65}}}}')
 for index,q in enumerate((1,2,4)):
  values=[float(r['speedup_syrk']) for r in rows if (r['evaluation_set'],r['m'],r['n'],r['q_cap'],r['output_full'])
          ==('heldout',str(m),str(n),str(q),'1')]
  assert len(values)==10
  x=base+.65+index*1.2
  lo,med,hi=(.75+1.08*v for v in (min(values),statistics.median(values),max(values)))
  p += [f'\\put({x:.3f},{lo:.3f}){{\\line(0,1){{{hi-lo:.3f}}}}}',
        f'\\put({x-.1:.3f},{lo:.3f}){{\\line(1,0){{.2}}}}',
        f'\\put({x-.1:.3f},{hi:.3f}){{\\line(1,0){{.2}}}}',
        f'\\put({x-.05:.3f},{med-.05:.3f}){{\\rule{{.10cm}}{{.10cm}}}}',
        f'\\put({x:.3f},.43){{\\makebox(0,0){{{q}}}}}']
 p.append(f'\\put({base+1.8:.2f},.10){{\\makebox(0,0){{defect cap $q$}}}}')
p += [r'\put(.70,5.05){\makebox(0,0)[l]{Paired SYRK/sequential speedup; full output}}',r'\end{picture}']
(out/'speedup.tex').write_text('\n'.join(p)+'\n')
print('generated selected_results.tex and speedup.tex from all230 retained rows')

# GitHub-readable visual evidence uses exactly the same retained observations.
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="390" viewBox="0 0 1080 390" role="img" aria-label="Measured SYRK over sequential speedup including losses">',
     '<rect width="1080" height="390" rx="12" fill="#f7fafc"/>',
     '<g font-family="Arial, sans-serif" fill="#142b40">',
     '<text x="35" y="33" font-size="21" font-weight="700">Exact sparse Gram: measured speedup over Accelerate SYRK</text>',
     '<text x="35" y="58" font-size="14">Full output · original frozen data · markers: median; bars: observed range across 10 runs</text>']
for panel,(m,n) in enumerate(((192,768),(384,2048),(768,8192))):
 base=65+panel*350
 svg.append(f'<text x="{base+120}" y="95" text-anchor="middle" font-size="18" font-weight="700">{m} × {n}</text>')
 for value in range(4):
  yy=310-60*value
  svg.append(f'<path d="M{base} {yy}h260" stroke="'+('#64748b" stroke-dasharray="5 4' if value==1 else '#d9e2ec')+'"/>')
  svg.append(f'<text x="{base-12}" y="{yy+5}" text-anchor="end" font-size="13">{value}×</text>')
 for index,q in enumerate((1,2,4)):
  values=[float(r['speedup_syrk']) for r in rows if (r['evaluation_set'],r['m'],r['n'],r['q_cap'],r['output_full'])==('heldout',str(m),str(n),str(q),'1')]
  x=base+40+index*85
  low,mid,high=(310-60*v for v in (min(values),statistics.median(values),max(values)))
  svg.append(f'<path d="M{x} {low}V{high}M{x-7} {low}h14M{x-7} {high}h14" stroke="#086b8a" stroke-width="2"/>')
  svg.append(f'<rect x="{x-5}" y="{mid-5}" width="10" height="10" fill="#086b8a"/><text x="{x}" y="{mid-13}" text-anchor="middle" font-size="13">{statistics.median(values):.2f}×</text>')
  svg.append(f'<text x="{x}" y="333" text-anchor="middle" font-size="14">q = {q}</text>')
svg += ['<text x="35" y="371" font-size="14">Apple M4 · binary64 · 1 thread · complete calls including admission, allocation and symmetric output</text>','</g></svg>']
(out/'speedup.svg').write_text('\n'.join(svg)+'\n')

# New specialized campaign is never pooled with the original 230 runs.
from collections import defaultdict
special=list(csv.DictReader((ROOT/'data/specialized_2026-10-19/summary.csv').open()))
assert len(special)==20
groups=defaultdict(list)
for r in special: groups[(r['family'],r['dist'])].append(r)
lines=[r'\begin{table}[h]',r'\centering\footnotesize\setlength{\tabcolsep}{3pt}',r'\begin{tabular}{llrll}',r'\toprule',
       r'Family & Coefficients & $R$ range & SYRK/seq. & H+E/seq. '+chr(92)*2,r'\midrule']
for (family,dist), rr in sorted(groups.items()):
 assert len(rr)==5
 label='Innovations' if family=='innov' else 'Windows'
 distribution='Uniform' if dist=='0' else 'Near constant'
 R=[int(r['actual_e']) for r in rr]
 ranges=[min(R),max(R)]
 rs=str(ranges[0]) if ranges[0]==ranges[1] else str(ranges[0])+'--'+str(ranges[1])
 line=f'{label} & {distribution} & {rs}'
 for key in ['0_over_13','14_over_13']:
  values=[float(r[key]) for r in rr]
  line+=f' & {statistics.median(values):.3f} [{min(values):.3f},{max(values):.3f}]'
 lines.append(line+' '+chr(92)*2)
lines += [r'\bottomrule',r'\end{tabular}',
 r'\caption{Separate new structured comparison: medians of five run-level paired ratios with observed minimum/maximum in brackets, not confidence intervals. All jobs use $512\times4096$, cap two, $e=-15$, full output. $R$ counts actual reconstructed Hankel corrections.}',
 r'\label{tab:special}',r'\end{table}']
lines += ['The sequential member beats both dense vendors and this reconstructed-Hankel implementation in all 20 new runs. It is about 1.7--1.8 times faster than the latter on these corrected windows, and about 42--46 times faster on accumulated innovations. The large latter ratio reflects the growth of reconstructed corrections, not a claim against every structured method. An optimized alternative Hankel reconstruction or a different published structured implementation remains an outstanding comparison. No structured state-of-the-art claim follows.']
(out/'specialized_results.tex').write_text('\n'.join(lines)+'\n')
