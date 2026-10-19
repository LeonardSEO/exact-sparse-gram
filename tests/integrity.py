#!/usr/bin/env python3
"""Verify byte-preserved sources/data/receipts against the release mapping."""
import hashlib
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
EXPECTED = '7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5'
assert hashlib.sha256((ROOT/'src/lattice_sequential.c').read_bytes()).hexdigest() == EXPECTED
mapping = json.loads((ROOT/'evidence/original/source_mapping.json').read_text())
for item in mapping:
    assert hashlib.sha256((ROOT/item['release']).read_bytes()).hexdigest() == item['sha256'], item
assert (ROOT/'src/lattice_sequential.c').read_bytes() == (ROOT/'benchmarks/round2/src/lattice_sequential.c').read_bytes()
print(f'PASS: {len(mapping)} byte-preserved sources, data, and historical receipts; frozen kernel hash matches.')
