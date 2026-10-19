# Frozen Round2 comparisons

All example rows: uniform signed16-bit coefficients, seed20261201, full output.
Times are per-method medians. Speedups are paired-ratio medians;
the 95% intervals resample five contiguous six-trial balanced blocks.
Intervals describe within-run variability, not cross-hardware coverage.

| Family | Shape | q cap | e | SYRK ms | GEMM ms | Earlier certified ms | Sequential ms | SYRK/seq (95% CI) | GEMM/seq (95% CI) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|innov|768×8192|2|-15|19.6336|20.5861|13.2805|9.0460|2.169 [2.164, 2.175]|2.277 [2.272, 2.282]|
|innov|768×8192|4|0|30.4593|31.1561|23.4610|15.2613|2.045 [1.889, 2.229]|2.075 [1.971, 2.290]|
|windows|512×4096|2|-15|3.2515|3.6209|2.9428|2.0573|1.596 [1.554, 1.612]|1.766 [1.732, 1.780]|
|innov|512×4096|2|-15|2.9476|3.3151|2.8378|1.9707|1.541 [1.438, 1.588]|1.690 [1.584, 1.817]|
|innov|384×2048|4|-15|1.2021|1.3621|1.9897|1.2632|0.952 [0.943, 0.982]|1.069 [1.063, 1.104]|
|innov|192×768|2|-15|0.0820|0.0858|0.2061|0.1353|0.612 [0.606, 0.614]|0.623 [0.598, 0.634]|

## Every group, including losses

Heldout groups combine two scales and five seeds (10 runs); distribution groups have five runs.
These are ranges and medians of per-run paired speedups, not pooled confidence intervals.
Lower/full, scaled versions and deterministic patterns are repeated inputs, not independent draws.

| Set | Family | Shape | q | dist | output | runs | SYRK/seq min, median, max | Best-vendor/seq min, median, max |
|---|---|---:|---:|---:|---|---:|---:|---:|
|distribution|innov|512×4096|2|0|full|5|1.533, 1.548, 1.567|1.533, 1.548, 1.567|
|distribution|innov|512×4096|2|1|full|5|1.455, 1.554, 1.598|1.455, 1.554, 1.598|
|distribution|innov|512×4096|2|2|full|5|1.635, 1.654, 1.749|1.635, 1.654, 1.749|
|distribution|innov|512×4096|2|3|full|5|1.837, 1.909, 2.002|1.837, 1.909, 2.002|
|distribution|innov|512×4096|2|4|full|5|1.471, 1.615, 1.641|1.471, 1.615, 1.641|
|distribution|windows|512×4096|2|0|full|5|1.520, 1.596, 1.629|1.520, 1.596, 1.629|
|distribution|windows|512×4096|2|1|full|5|1.610, 1.660, 1.694|1.610, 1.660, 1.694|
|distribution|windows|512×4096|2|2|full|5|1.588, 1.644, 1.657|1.588, 1.644, 1.657|
|distribution|windows|512×4096|2|3|full|5|1.886, 1.930, 1.949|1.886, 1.930, 1.949|
|distribution|windows|512×4096|2|4|full|5|1.654, 1.717, 1.877|1.654, 1.717, 1.877|
|heldout|innov|192×768|1|0|lower|10|0.634, 0.641, 0.665|0.634, 0.641, 0.665|
|heldout|innov|192×768|1|0|full|10|0.674, 0.684, 0.691|0.674, 0.684, 0.694|
|heldout|innov|192×768|2|0|lower|10|0.570, 0.575, 0.654|0.570, 0.575, 0.654|
|heldout|innov|192×768|2|0|full|10|0.592, 0.608, 0.612|0.592, 0.608, 0.612|
|heldout|innov|192×768|4|0|lower|10|0.455, 0.464, 0.548|0.455, 0.464, 0.548|
|heldout|innov|192×768|4|0|full|10|0.489, 0.497, 0.507|0.489, 0.496, 0.507|
|heldout|innov|384×2048|1|0|lower|10|1.198, 1.212, 1.222|1.198, 1.212, 1.222|
|heldout|innov|384×2048|1|0|full|10|1.197, 1.204, 1.277|1.197, 1.204, 1.277|
|heldout|innov|384×2048|2|0|lower|10|1.073, 1.112, 1.246|1.073, 1.112, 1.246|
|heldout|innov|384×2048|2|0|full|10|1.102, 1.116, 1.145|1.102, 1.116, 1.145|
|heldout|innov|384×2048|4|0|lower|10|0.924, 0.951, 0.998|0.924, 0.951, 0.998|
|heldout|innov|384×2048|4|0|full|10|0.952, 0.958, 0.971|0.952, 0.958, 0.971|
|heldout|innov|768×8192|1|0|lower|10|2.303, 2.496, 2.573|2.303, 2.496, 2.573|
|heldout|innov|768×8192|1|0|full|10|2.418, 2.477, 2.840|2.418, 2.477, 3.083|
|heldout|innov|768×8192|2|0|lower|10|1.958, 2.204, 2.296|1.958, 2.204, 2.296|
|heldout|innov|768×8192|2|0|full|10|1.917, 2.166, 2.348|1.917, 2.166, 2.348|
|heldout|innov|768×8192|4|0|lower|10|1.610, 1.826, 1.867|1.610, 1.826, 1.867|
|heldout|innov|768×8192|4|0|full|10|1.613, 1.811, 2.589|1.613, 1.811, 2.551|

230 manifest-checked runs; 150 paired-median wins against both vendors.
150 runs also have both conditional interval lower endpoints above1.
All methods11/12/13, individual times, vendor choices and intervals remain in lattice_summary.csv.
No per-case candidate selection is used: method13 is fixed throughout.
