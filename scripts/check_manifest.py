#!/usr/bin/env python3
"""Verify the frozen source release SHA256SUMS; no dependency caches required."""
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
n=0
for line in (root/'SHA256SUMS').read_text().splitlines():
 digest,name=line.split('  ',1)
 p=root/name
 assert p.is_file(),name
 assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
 n+=1
print(f'PASS: {n} release files match SHA256SUMS')
