#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
python3 tests/integrity.py
python3 tests/sequential_verify.py
