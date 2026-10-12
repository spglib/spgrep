# v0.8.0 performance comparison

Table 1 compares v0.7.0 with the v0.8.0 candidate at `4f7f871`. Both used Apple M1, macOS 26.5, Python 3.13.5, NumPy 2.5.2 with Accelerate, spglib 2.7.0, and single-thread limits. Each timing is the median of nine batch means, with 20 warm calls per batch and alternating revision order.

Table 1. Execution time in milliseconds per call. Speedup is the ratio of the v0.7.0 and candidate medians.

| Workload | v0.7.0, ms | v0.8.0 candidate, ms | Speedup |
| --- | ---: | ---: | ---: |
| Cubic crystal, Neto | 16.261 | 8.063 | 2.02× |
| Cubic crystal, random | 25.480 | 9.337 | 2.73× |
| Cubic spinor, Neto | 58.700 | 6.437 | 9.12× |
| P4₂/mnm at Z, Neto | 2.795 | 2.040 | 1.37× |
| P4₂/mnm at Z, random | 2.835 | 1.356 | 2.09× |
| Cayley table, C₃ᵥ, order 6 | 0.096 | 0.017 | 5.63× |
| Cayley table, Oₕ, order 48 | 6.098 | 0.541 | 11.27× |
| Cayley table, gray Oₕ, order 96 | 29.786 | 2.105 | 14.15× |

Cubic-crystal calls use a one-atom cubic cell at Γ and include symmetry search. Spinor calls use the 48 operations of Oₕ. The nonsymmorphic cases use supplied P4₂/mnm operations at k = (0, 0, ½). Imports and input preparation are untimed. Cayley-table timings measure only table construction.

These results apply to the selected workloads and environment. Peak memory was not measured; batching uses additional temporary arrays. [Recorded samples and metadata](https://github.com/spglib/spgrep/blob/838a95ce79eb8cf875cb056244fce969935fce90/benchmarks/results/v0.8.0.json) retain the measurement details.
