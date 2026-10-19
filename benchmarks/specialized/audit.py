"""Shared untimed generator/order and receipt validation for the frozen protocol."""
import math
import struct

METHODS=('0','4','13','14')
MASK=(1<<64)-1
FIELDS={'family':'family','m':'m','n':'n','q':'q_cap','e':'scale_exp','seed':'seed','dist':'dist','full':'output_full'}

class RNG:
    def __init__(self,seed):self.state=seed or 1
    def rnd(self):
        x=self.state;x^=(x<<13)&MASK;x^=x>>7;x^=(x<<17)&MASK;self.state=x&MASK
        return self.state

def reference_input(job):
    rng=RNG(job['seed']);m,n,q=job['m'],job['n'],job['q']
    def coeff(k):
        if job['dist']==0:return (rng.rnd()&65535)-32768
        if job['dist']==1:return 16384+(rng.rnd()%17)-8
        if job['dist']==2:return -32768 if k&1 else 32767
        if job['dist']==3:return 0 if k%17 else (-32768 if k&1 else 32767)
        return (rng.rnd()%3)-1
    if job['family']=='innov':
        rows=[[coeff(k) for k in range(n)]]
        for i in range(1,m):
            row=rows[-1][1:]+[coeff(i+n-1)]
            for _ in range(q):
                k=rng.rnd()%(n-1);row[k]=coeff(i*n+k)
            rows.append(row)
    else:
        h=[coeff(k) for k in range(m+n-1)]
        rows=[h[i:i+n] for i in range(m)]
        if q:
            for i in range(m):rows[i][(17*i+7)%n]=0
    hash_value=14695981039346656037
    for row in rows:
        for x in row:
            for byte in struct.pack('<d',math.ldexp(float(x),job['e'])):
                hash_value=((hash_value^byte)*1099511628211)&MASK
    return f'{hash_value:016x}',rng.state

def reference_orders(seed,samples):
    rng=RNG((seed+713)&MASK);orders=[]
    # C state assignment can be zero (unlike generator's seed-or-one).
    rng.state=(seed+713)&MASK
    for s in range(samples):
        if s%4==0:
            base=list(range(4))
            for z in range(3,0,-1):
                j=rng.rnd()%(z+1);base[z],base[j]=base[j],base[z]
        orders.append([int(METHODS[base[(z+s%4)%4]]) for z in range(4)])
    return orders

def validate_input(record,job,samples):
    for key,field in FIELDS.items():
        if record[field]!=job[key]:raise ValueError(f'CLI/input mismatch {field}: {record[field]} != {job[key]}')
    fp=record['input_fingerprint_fnv1a64']
    if len(fp)!=16 or any(c not in '0123456789abcdef' for c in fp):raise ValueError('Invalid input fingerprint')
    if record['method_orders']!=reference_orders(job['seed'],samples):raise ValueError('Seeded method order mismatch')
    for block in range(0,samples,4):
        for position in range(4):
            if sorted(record['method_orders'][s][position] for s in range(block,block+4))!=[0,4,13,14]:raise ValueError('Unbalanced method order')

def validate_record(record,job,manifest,index,manifest_sha):
    validate_input(record,job,manifest['samples'])
    if record['job_index']!=index or record['phase']!=manifest['phase']:raise ValueError('Index/phase mismatch')
    if record['manifest_sha256']!=manifest_sha or record['source_binary_sha256']!=manifest['source_binary_sha256']:raise ValueError('Receipt hash mismatch')
    if record['input_fingerprint_fnv1a64']!=manifest['input_diagnostics'][index]['input_fingerprint_fnv1a64']:raise ValueError('Frozen input mismatch')
    if record['blas_threads']!=1:raise ValueError('BLAS threads mismatch')
    if len(record['samples'])!=manifest['samples']:raise ValueError('Sample count mismatch')
    for sample in record['samples']:
        if set(sample)!=set(METHODS) or any(not isinstance(v,(int,float)) or not math.isfinite(v) or v<=0 for v in sample.values()):raise ValueError('Invalid measurement sample')
    if set(record['batches'])!=set(METHODS) or any(not isinstance(v,int) or v<1 for v in record['batches'].values()):raise ValueError('Invalid calibration batches')
    if record['max_observed_q']>job['q'] or not 0<=record['actual_q']<=(job['m']-1)*job['q']:raise ValueError('Defect cap mismatch')
    if not 0<=record['actual_e']<=job['m']*job['n']:raise ValueError('Invalid E count')

def validate_complete_records(records,manifest,manifest_sha):
    if len(records)!=len(manifest['jobs']):raise ValueError('Expected complete 20-job grid')
    if len({r['job_index'] for r in records})!=len(manifest['jobs']):raise ValueError('Duplicate/missing job index')
    for index,(record,job) in enumerate(zip(records,manifest['jobs'],strict=True)):
        validate_record(record,job,manifest,index,manifest_sha)
