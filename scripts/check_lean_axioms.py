#!/usr/bin/env python3
"""Fail closed on incomplete, missing, duplicate or unsupported Lean axiom reports."""
import argparse
from pathlib import Path
import re

ALLOWED = {"propext", "Classical.choice", "Quot.sound"}
EXPECTED = {"GramExact." + name for name in """
mixed
firstState
secondState
support
dot
boundary
update
reconstruct
first_prefix
second_prefix
first_step
second_step
groups_join
groups_complete
boundary_identity
recurrence
reconstruction
product_subset_bound
tail_sum_bound
first_prefix_bound
second_prefix_bound
first_prefix_cap
second_prefix_cap
residual_bound
support_work
total_work
coefficientCap
GridRepresentable
RoundingModel
Binary64Bridge
Binary64SafeDyadic
roundedFMA_coefficients
roundedSub_coefficients
dyadicScale_square
accumulator_exponent_safe
roundedFirst
roundedSecond
rounded_first_exact
rounded_second_exact
rounded_groups_exact
rounded_boundary_exact
rounded_anchor_exact
rounded_anchor_merge
rounded_residual_exact
roundedUpdate
rounded_update_exact
roundedReconstruct
rounded_reconstruction
computed_anchor_exact
roundedGram
rounded_gram_exact
scaled_dot
binary64_gram_exact
""".split()}

def check(path):
    content = path.read_text()
    reports = re.findall(r"^'([^']+)' depends on axioms: \[([^\]]*)\]$", content, re.MULTILINE)
    found = {}
    for name, raw in reports:
        if name in found:
            raise ValueError(f"Duplicate axiom report: {name}")
        axioms = set(raw.split(", ")) if raw else set()
        if not axioms <= ALLOWED:
            raise ValueError(f"Unsupported axioms for {name}: {sorted(axioms - ALLOWED)}")
        found[name] = axioms
    if set(found) != EXPECTED:
        raise ValueError(f"Axiom report mismatch: missing={sorted(EXPECTED-set(found))}; unexpected={sorted(set(found)-EXPECTED)}")
    if content.count("depends on axioms:") != len(reports):
        raise ValueError("Malformed axiom report")
    if re.search(r"(?im)^.*(?:error:|warning:.*sorry|\bsorryAx\b).*$", content):
        raise ValueError("Lean error or incomplete proof diagnostic")
    print(f"PASS: all {len(found)} expected principal axiom closures contain only standard logical axioms.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        check(args.output)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"FAIL: {exc}\n")
