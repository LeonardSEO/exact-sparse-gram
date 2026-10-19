CC = clang
CFLAGS = -O3 -std=c11 -fno-fast-math -ffp-contract=off -Wall -Wextra -Werror
ROOT = ../..
.PHONY: all check clean
all: build/specialized_benchmark build/test_specialized
build:
	mkdir -p build
build/specialized_benchmark: src/benchmark.c src/hankel_sparse.c src/hankel_sparse.h $(ROOT)/src/lattice_sequential.c | build
	$(CC) $(CFLAGS) -Isrc src/benchmark.c src/hankel_sparse.c $(ROOT)/src/lattice_sequential.c -framework Accelerate -o $@
build/test_specialized: tests/test_specialized.c src/hankel_sparse.c src/hankel_sparse.h $(ROOT)/src/lattice_sequential.c | build
	$(CC) $(CFLAGS) -Isrc tests/test_specialized.c src/hankel_sparse.c $(ROOT)/src/lattice_sequential.c -o $@
check: all
	./build/test_specialized
	python3 tests/test_harness.py
	python3 tests/test_receipts.py
clean:
	rm -rf build
