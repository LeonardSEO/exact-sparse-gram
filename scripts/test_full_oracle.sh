#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
if [ "$(uname -s)" != Darwin ]; then
  printf '%s\n' 'SKIP: full-entry vendor oracle requires macOS Accelerate; portable exact tests are scripts/test.sh.'
  exit 0
fi
python3 tests/algebra_lattice_full_entry_verify.py --run
