#!/usr/bin/env python3
"""No performance calls: diagnostic branch exits before warmup/calibration/clock."""
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit import reference_input,validate_input
from run import BINARY

def cmd(job):
    return [str(BINARY),job['family'],str(job['m']),str(job['n']),str(job['q']),str(job['e']),str(job['seed']),str(job['dist']),'32',str(job['full']),'10','--input-diagnostic']

def main():
    records=[]
    for family in ('windows','innov'):
        for dist in (0,1):
            seen=set()
            for seed in range(20261201,20261206):
                job=dict(family=family,m=17,n=37,q=2,e=-15,seed=seed,dist=dist,full=1)
                record=json.loads(subprocess.check_output(cmd(job),text=True))
                validate_input(record,job,32)
                expected,state=reference_input(job)
                assert record['input_fingerprint_fnv1a64']==expected
                assert record['generator_final_state']==state
                assert expected not in seen;seen.add(expected)
                records.append({k:record[k] for k in ('family','dist','seed','input_fingerprint_fnv1a64','generator_final_state')})
            assert len(seen)==5
    base=cmd(job)
    for index,value in ((2,'-1'),(3,'37junk'),(4,'2x'),(5,'-15junk'),(6,'20261201x'),(7,'1x'),(8,'32x'),(9,'1x'),(10,'10x')):
        invalid=base.copy();invalid[index]=value
        assert subprocess.run(invalid,capture_output=True).returncode==2
    print(json.dumps(dict(status='pass',mode='input diagnostic; no clock or kernel calls',records=records),indent=2))

if __name__=='__main__':main()
