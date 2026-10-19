#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
mkdir -p build
CC=${CC:-clang}
CFLAGS='-O3 -std=c11 -fno-fast-math -ffp-contract=off -Wall -Wextra -Wconversion -Werror'
$CC $CFLAGS -Iinclude -c src/lattice_sequential.c -o build/lattice_sequential.o
ar rcs build/libgram_exact.a build/lattice_sequential.o
printf '%s\n' 'Built build/libgram_exact.a (portable core; no benchmark timing).'
if [ "${1:-}" = '--benchmarks' ]; then
  [ "$(uname -s)" = Darwin ] || { printf '%s\n' 'Six-method benchmark requires macOS Accelerate.' >&2; exit 1; }
  S=benchmarks/round2/src
  $CC -O3 -std=c11 -fno-fast-math -ffp-contract=off -Wall -Wextra \
    "$S/lattice_benchmark.c" "$S/joint_kernel.c" "$S/communication_skinny.c" \
    "$S/adapters.c" "$S/algebra_hankel.c" "$S/numerical_guard.c" \
    "$S/algebra_sparse.c" "$S/numerical_sparse.c" "$S/communication_guard.c" \
    "$S/communication_guard_fast.c" "$S/numerical_sparse_fast.c" \
    "$S/lattice_exact.c" "$S/lattice_substitute.c" "$S/lattice_sequential.c" \
    -framework Accelerate -o build/lattice_benchmark
  printf '%s\n' 'Built build/lattice_benchmark. Run timing only without competing compute.'
fi
