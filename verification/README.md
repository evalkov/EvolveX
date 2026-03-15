# Verification Suite

Tools to verify that the compute optimizations produce correct results and to measure the performance improvement.

## Quick Start

```bash
./verification/benchmark.sh /path/to/foldx 4
```

This runs the full benchmark: old code vs optimized code on the Vsig4 example, compares outputs, validates the Raw_ file parsing, and prints a timing summary.

## What Was Optimized

The GA search hot loop was the bottleneck — each of hundreds of thousands of MC steps invoked 4–6 FoldX subprocesses sequentially, performed redundant file I/O, and compressed PDB files synchronously. The optimizations are:

| Optimization | Description | Estimated Impact |
|---|---|---|
| Eliminate Stability calls | Read `Raw_`/`Dif_` files from BuildModel instead of running the Stability command | 33–50% of FoldX time |
| Consolidate file reads | `parse_interaction_file()` reads once instead of 3× per file | ~10ms/step, millions of reads saved |
| Early rejection | Check stability criteria before reading mutant interaction files | 0–50% of remaining I/O |
| Fix random sampling | Pre-filter to allowed AAs instead of rejection-sampling all 20 | Eliminates worst-case loops |
| Parallelize AnalyseComplex | Run 2–4 independent FoldX calls concurrently via `ThreadPoolExecutor` | 1.5–2× on remaining FoldX time |
| Defer PDB compression | Copy during search, batch compress after | 5–15ms/step (40–125 min total) |
| Cap combinatorial product | Random-sample if >1M combinations instead of materializing | Prevents OOM |
| Dict-based residue IDs | `dict[str, str]` instead of `list[str]` with O(n) `.remove()` | O(1) updates |

**Expected net effect:** ~50–70% wall time reduction on the GA search phase.

## Scripts

### `benchmark.sh`

Full end-to-end comparison of old vs new code.

```bash
./verification/benchmark.sh <path_to_foldx_dir> [n_cores]
```

**Arguments:**
- `path_to_foldx_dir` — directory containing the `foldx` binary (required)
- `n_cores` — number of CPU cores for Dask (default: 4)

**What it does:**
1. Creates a git worktree of the pre-optimization commit
2. Generates two identical Vsig4 configs pointing to isolated working directories
3. Runs the OLD code with wall-time measurement
4. Runs the NEW code with wall-time measurement
5. Runs `compare_outputs.py` to check correctness
6. Runs `validate_raw_files.py` to verify the Raw_ file parsing
7. Prints a speedup summary
8. Cleans up the worktree

**Output:** All logs and reports are saved to `verification/benchmark_run_<timestamp>/`.

### `compare_outputs.py`

Compares two `generated_models_info.csv` files from old and new runs.

```bash
python verification/compare_outputs.py <old_csv> <new_csv>
```

**Checks performed:**
1. **Schema** — same columns in both files
2. **Row count** — same number of iterations
3. **Column types** — energy columns are numeric
4. **Value ranges** — energies are within physically reasonable bounds
5. **Metadata structure** — step types (`MC`, `recombination`), backbone names, residue ID format
6. **Distribution comparison** — side-by-side mean/std of all energy columns

Since the GA search is stochastic, exact value matches are not expected. The check confirms that both runs produce structurally correct output with energy values in the same ballpark.

### `validate_raw_files.py`

Validates that the Raw_/Dif_ file approach produces the same complex stability values as the explicit FoldX Stability command. This is the critical correctness proof for the largest optimization.

```bash
python verification/validate_raw_files.py <evolvex_working_dir> <foldx_dir>
```

**What it does:**
1. Finds a PDB backbone in the EvolveX working directory
2. Runs BuildModel to generate `Raw_` and `Dif_` files
3. Reads mutant total energy from `Raw_`, ddG from `Dif_`, derives wildtype energy
4. Runs FoldX Stability on both mutant and wildtype PDBs
5. Compares values with a 0.01 kcal/mol tolerance
6. Cleans up temporary files

**Expected result:** PASS — the values should match because both sources derive from the same FoldX force field evaluation. The `Raw_` file contains the total energy computed during BuildModel, which is the same quantity that the Stability command recomputes.

## Interpreting Results

### Timing

The benchmark prints:
```
OLD (pre-optimization): 1200s
NEW (optimized):         500s
Speedup:                2.40x
Time reduction:         58.3%
```

The actual speedup depends on:
- Number of backbones and models (more = better amortization of parallelism overhead)
- FoldX binary performance on your hardware
- Fraction of mutations rejected on stability (more rejections = bigger benefit from early rejection)
- Whether `calculate_binding_dG_with_water` is enabled (more AnalyseComplex calls = bigger benefit from parallelization)

### Output Comparison

The comparison report shows per-column statistics. What to look for:
- **All checks PASS** — the optimized code produces structurally valid output
- **Energy means/stds are in the same range** — the search is exploring the same energy landscape
- **Exact matches are NOT expected** — the GA is stochastic, and internal parallelism ordering can vary

### Raw_ Validation

This should always PASS. If it fails, it means the Raw_ file format differs from what was assumed, and the optimization in `keep_mutant_decision()` would produce incorrect complex stability values. In that case, revert Phase 1A (re-enable the Stability calls in `run_foldx_commands()`).
