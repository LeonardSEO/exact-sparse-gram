#!/usr/bin/env python3
"""Synthetic in-memory schema fixtures only; never produce measurement results."""
import copy
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit import reference_orders,validate_complete_records,FIELDS
from run import planned,validate_manifest

def rejected(callback):
    try:callback()
    except (ValueError,KeyError):return
    raise AssertionError('Invalid receipt accepted')

def main():
    manifest=planned();manifest['source_binary_sha256']={'synthetic':'a'*64}
    manifest['input_diagnostics']=[];records=[];manifest_sha='b'*64
    for index,job in enumerate(manifest['jobs']):
        raw={field:job[key] for key,field in FIELDS.items()}
        raw.update(input_fingerprint_fnv1a64=f'{index:016x}',method_orders=reference_orders(job['seed'],32))
        manifest['input_diagnostics'].append(copy.deepcopy(raw))
        raw.update(job_index=index,phase=manifest['phase'],manifest_sha256=manifest_sha,
                   source_binary_sha256=manifest['source_binary_sha256'],blas_threads=1,
                   samples=[{key:1. for key in ('0','4','13','14')} for _ in range(32)],
                   batches={key:1 for key in ('0','4','13','14')},max_observed_q=2,actual_q=1022,actual_e=1022)
        records.append(raw)
    validate_manifest(manifest);validate_complete_records(records,manifest,manifest_sha)
    rejected(lambda:validate_complete_records(records[:-1],manifest,manifest_sha))
    mutations=[('job_index',1),('seed',2),('manifest_sha256','c'*64),('source_binary_sha256',{}),
               ('input_fingerprint_fnv1a64','f'*16),('blas_threads',2),('samples',[]),('method_orders',[[0,4,13,14]]*32)]
    for key,value in mutations:
        bad=copy.deepcopy(records);bad[0][key]=value
        rejected(lambda:validate_complete_records(bad,manifest,manifest_sha))
    bad=copy.deepcopy(records);bad[0]['samples'][0]['13']=float('nan')
    rejected(lambda:validate_complete_records(bad,manifest,manifest_sha))
    bad=copy.deepcopy(records);bad[1]=copy.deepcopy(bad[0])
    rejected(lambda:validate_complete_records(bad,manifest,manifest_sha))
    print('receipt validation: complete synthetic grid accepted; missing/duplicate/seed/hash/input/order/thread/sample mutations rejected')

if __name__=='__main__':main()
