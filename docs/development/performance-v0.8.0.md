# Execution-time comparison of spgrep v0.7.0 and v0.8.0

## Abstract

Execution times were compared between spgrep v0.7.0 and the v0.8.0 release candidate using identical Python dependencies and single-thread limits. Eight workloads covered Cayley-table construction, ordinary and spinor irreps, and a nonsymmorphic Brillouin-zone boundary. The default cubic-crystal irrep call decreased from 16.261 to 8.063 ms, while the cubic spinor call decreased from 58.700 to 6.437 ms. All measured workloads improved. These measurements quantify warm-call performance for the selected inputs on an Apple M1 processor.

## Introduction

The v0.8.0 candidate includes five performance changes to irrep extraction, Cayley-table construction, spinor factors, nonsymmorphic factors, and regular-representation lookup [1]. The comparison measures their combined effect on public irrep calls and distinguishes these timings from the cost of Cayley-table construction alone. Both Neto and random enumeration methods are included because they exercise different decomposition paths.

## Methods

The baseline was v0.7.0 at commit `207aaf6e2a1057357a929819839760b238f1680b`. The candidate was commit `4f7f871c2de5fa86e64c1393c8640c47e8a481b0`, which contains all five merged optimizations. Source archives from these commits were imported by the same benchmark script [2]. This compares source revisions with fixed dependencies rather than independently resolved release environments.

Measurements used an Apple M1 processor running macOS 26.5 on ARM64, Python 3.13.5, NumPy 2.5.2 with Apple Accelerate BLAS and LAPACK, and spglib 2.7.0. `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and `VECLIB_MAXIMUM_THREADS` were set to `1`. Nine paired batches alternated which revision ran first. Each batch started a fresh Python process and measured the mean time of 20 calls after validation and one warmup call. Tables report the median and interquartile range of the nine batch means. Speedup is the ratio of baseline and candidate medians.

The cubic-crystal workload uses a unit cubic lattice with one atom at the origin and a zero k-point. Its timed `get_spacegroup_irreps` call includes the spglib symmetry search. The spinor workload calls `get_crystallographic_pointgroup_spinor_irreps_from_symmetry` for the 48 operations of $O_h$, including spinor preparation. The nonsymmorphic workload uses the 16 operations of $P4_2/mnm$, Hall number 419, at $\mathbf{k}=(0,0,1/2)$. Its `get_spacegroup_irreps_from_primitive_symmetry` call starts from supplied operations and includes little-group selection, factor construction, and irrep enumeration.

Cayley-table workloads use $C_{3v}$, $O_h$, and the gray group $O_h1'$, with orders 6, 48, and 96. Imports, input construction, database lookup for the supplied nonsymmorphic operations, validation, and warmup are outside every timed region. Each process checks Cayley-table rows as permutations, or irrep completeness and character orthogonality. The nonsymmorphic cases also require four two-dimensional irreps, matching the existing regression fixture.

## Results

### Public irrep calls

The combined optimizations reduce execution time for all five public-call workloads in Table 1. The largest measured ratio occurs in the cubic spinor case. The nonsymmorphic boundary-point calls also improve, with a smaller reduction for the Neto method than for the random method.

Table 1. Public-call execution times in milliseconds. Brackets give the interquartile range of batch means. Speedup is the v0.7.0 median divided by the candidate median. Cubic-crystal calls include symmetry search; the other calls start from supplied symmetry operations.

| Workload | v0.7.0, ms | v0.8.0 candidate, ms | Speedup |
| --- | ---: | ---: | ---: |
| Cubic crystal, Neto | 16.261 [16.170, 16.456] | 8.063 [8.008, 8.173] | 2.02× |
| Cubic crystal, random | 25.480 [25.351, 25.715] | 9.337 [9.299, 9.413] | 2.73× |
| Cubic spinor, Neto | 58.700 [58.564, 58.846] | 6.437 [6.428, 6.463] | 9.12× |
| $P4_2/mnm$ at Z, Neto | 2.795 [2.772, 2.808] | 2.040 [2.039, 2.063] | 1.37× |
| $P4_2/mnm$ at Z, random | 2.835 [2.829, 2.860] | 1.356 [1.340, 1.361] | 2.09× |

### Cayley-table construction

Table 2 isolates table construction, whose implementation batches matrix products and uses indexed operation lookup. These kernel timings measure one component of the public-call workloads. They do not represent complete irrep calculations.

Table 2. Cayley-table construction times in milliseconds, with interquartile ranges of batch means in brackets. The order includes antiunitary operations in the gray group. Speedup is the ratio of medians.

| Group | Order | v0.7.0, ms | v0.8.0 candidate, ms | Speedup |
| --- | ---: | ---: | ---: | ---: |
| $C_{3v}$ | 6 | 0.096 [0.095, 0.097] | 0.017 [0.017, 0.017] | 5.63× |
| $O_h$ | 48 | 6.098 [6.051, 6.152] | 0.541 [0.538, 0.541] | 11.27× |
| $O_h1'$ | 96 | 29.786 [29.704, 29.934] | 2.105 [2.099, 2.111] | 14.15× |

## Discussion

The public-call measurements compare the combined release changes, with each workload exercising the applicable subset. They do not assign a separate contribution to each change. The ratios depend on the group, k-point, enumeration method, hardware, and numerical libraries. The interquartile ranges describe variation within this run, not uncertainty across machines or input distributions. Startup cost, peak memory, magnetic co-representation enumeration, and exhaustive coverage of space groups were not measured.

Batching increases temporary memory. Cayley-table and spinor intermediates scale quadratically with group order, while regular-representation lookup adds a boolean array with the same cubic shape as the existing dense representation. At order 48, that boolean array occupies 108 KiB. Public function signatures and dependency requirements are unchanged. Invalid empty regular representations now raise `ValueError` instead of returning an empty list.

## Conclusions

The candidate reduces warm execution time across the selected ordinary, spinor, and nonsymmorphic workloads. The results support a performance release while leaving workload-specific benchmarking necessary for applications outside this set.

## References

1. [Source changes from v0.7.0 to the measured candidate](https://github.com/spglib/spgrep/compare/v0.7.0...4f7f871c2de5fa86e64c1393c8640c47e8a481b0).
2. {download}`Benchmark worker <../../benchmarks/benchmark_cayley_table.py>`, {download}`comparison runner <../../benchmarks/compare_releases.py>`, and {download}`raw measurements and metadata <../../benchmarks/results/v0.8.0.json>`.

## Appendix: reproduction

The comparison runner extracts both source revisions into temporary directories and verifies the imported package paths. It uses the active interpreter and dependencies for both revisions. The JSON artifact records exact source commits, benchmark and runner SHA-256 hashes, thread limits, environment versions, and every batch mean. The release-preparation changes add documentation and benchmarks without changing the measured library source.

From the repository root, with the dependencies listed in Methods installed, the measurement command is:

```shell
python benchmarks/compare_releases.py \
    --baseline v0.7.0 \
    --candidate 4f7f871c2de5fa86e64c1393c8640c47e8a481b0 \
    --repeat 9 --number 20 \
    --output .cache/v0.8.0-comparison.json
```

Processor identification used `sysctl -n machdep.cpu.brand_string`. The numerical-library backend was checked with `numpy.show_config()`. The comparison excludes dependency upgrades in the release lockfile by using the same installed versions for both revisions.
