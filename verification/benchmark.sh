#!/usr/bin/env bash
#
# EvolveX Optimization Verification Benchmark
#
# Runs the old (pre-optimization) and new (optimized) code on the same Vsig4 config,
# then compares wall time and output correctness.
#
# Usage:
#   ./verification/benchmark.sh <path_to_foldx_dir> [n_cores]
#
# Prerequisites:
#   - FoldX binary installed and accessible
#   - Python environment with EvolveX dependencies (dask, biopython, pyyaml, pandas)
#   - Git (for worktree checkout of old code)
#
# What it does:
#   1. Creates a git worktree of the pre-optimization commit (ffacd2d)
#   2. Generates two separate config files pointing to isolated working directories
#   3. Runs OLD code with wall-time measurement
#   4. Runs NEW code with wall-time measurement
#   5. Runs compare_outputs.py to validate correctness
#   6. Prints timing summary
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OLD_COMMIT="ffacd2d"  # Last commit before optimization
NEW_COMMIT="6c6c95c"  # Optimization commit

# --- Arguments ---
if [ $# -lt 1 ]; then
    echo "Usage: $0 <path_to_foldx_dir> [n_cores]"
    echo "  path_to_foldx_dir: directory containing the 'foldx' binary"
    echo "  n_cores: number of cores for Dask (default: 4)"
    exit 1
fi

FOLDX_DIR="$(realpath "$1")"
N_CORES="${2:-4}"

if [ ! -x "$FOLDX_DIR/foldx" ]; then
    echo "ERROR: $FOLDX_DIR/foldx not found or not executable"
    exit 1
fi

# --- Setup directories ---
BENCH_DIR="$REPO_ROOT/verification/benchmark_run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BENCH_DIR"

OLD_WORKDIR="$BENCH_DIR/old_run"
NEW_WORKDIR="$BENCH_DIR/new_run"
OLD_WORKTREE="$BENCH_DIR/old_code"
mkdir -p "$OLD_WORKDIR" "$NEW_WORKDIR"

echo "=== EvolveX Optimization Benchmark ==="
echo "Benchmark directory: $BENCH_DIR"
echo "FoldX directory:     $FOLDX_DIR"
echo "Cores:               $N_CORES"
echo ""

# --- Create config files ---
# Both configs are identical except for working_dir, so they produce comparable runs
create_config() {
    local workdir="$1"
    local config_path="$2"
    cat > "$config_path" <<YAML
# Auto-generated benchmark config
antibody_chains: "N"
antigen_chains: "L"

working_dir: "$workdir"
foldx_dir: "$FOLDX_DIR"
Backbones_dir: "$REPO_ROOT/Vsig4_example/Backbones"
PositionsToExplore_file_path: "$REPO_ROOT/PositionsToExplore_Vsig4.tsv"

search_algorithm: "GA"
max_iterations: 10
population_size: 2
recombine_every_nth_iteration: 5

compute_env: "local"
n_cores: $N_CORES

residues_to_ignore: "GMHC"
vdwDesign: 2
print_stdout: false
calculate_binding_dG_with_water: true
YAML
}

create_config "$OLD_WORKDIR" "$BENCH_DIR/config_old.yaml"
create_config "$NEW_WORKDIR" "$BENCH_DIR/config_new.yaml"

# --- Checkout old code ---
echo "--- Setting up old code worktree at $OLD_COMMIT ---"
git -C "$REPO_ROOT" worktree add "$OLD_WORKTREE" "$OLD_COMMIT" 2>/dev/null || {
    echo "Worktree already exists or error; cleaning up and retrying..."
    git -C "$REPO_ROOT" worktree remove "$OLD_WORKTREE" --force 2>/dev/null || true
    git -C "$REPO_ROOT" worktree add "$OLD_WORKTREE" "$OLD_COMMIT"
}
echo ""

# --- Detect Python and site-packages ---
# The worktree only contains EvolveX source code; third-party dependencies (dask,
# biopython, pandas, etc.) come from the active virtual environment / system install.
# We prepend the OLD or NEW src/ to PYTHONPATH so that `import evolvex` resolves to
# the correct version, while all other packages resolve normally.
SITE_PACKAGES="$(python3 -c 'import site; print(":".join(site.getsitepackages()))')"
BASE_PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$SITE_PACKAGES"

# --- Run OLD code ---
echo "--- Running OLD code (pre-optimization) ---"
OLD_START=$(date +%s)
PYTHONPATH="$OLD_WORKTREE/src:$BASE_PYTHONPATH" python3 -m evolvex.main "$BENCH_DIR/config_old.yaml" 2>&1 | tee "$BENCH_DIR/old_stdout.log"
OLD_END=$(date +%s)
OLD_ELAPSED=$((OLD_END - OLD_START))
echo "OLD code finished in ${OLD_ELAPSED}s"
echo ""

# --- Run NEW code ---
echo "--- Running NEW code (optimized) ---"
NEW_START=$(date +%s)
PYTHONPATH="$REPO_ROOT/src:$BASE_PYTHONPATH" python3 -m evolvex.main "$BENCH_DIR/config_new.yaml" 2>&1 | tee "$BENCH_DIR/new_stdout.log"
NEW_END=$(date +%s)
NEW_ELAPSED=$((NEW_END - NEW_START))
echo "NEW code finished in ${NEW_ELAPSED}s"
echo ""

# --- Compare outputs ---
echo "--- Comparing outputs ---"
python3 "$REPO_ROOT/verification/compare_outputs.py" \
    "$OLD_WORKDIR/generated_models_info.csv" \
    "$NEW_WORKDIR/generated_models_info.csv" \
    2>&1 | tee "$BENCH_DIR/comparison_report.txt"
echo ""

# --- Run Raw_ vs Stability validation (on new run's intermediate files, if available) ---
echo "--- Validating Raw_ file parsing ---"
python3 "$REPO_ROOT/verification/validate_raw_files.py" \
    "$NEW_WORKDIR" "$FOLDX_DIR" \
    2>&1 | tee "$BENCH_DIR/raw_validation_report.txt"
echo ""

# --- Timing summary ---
echo "=========================================="
echo "         TIMING SUMMARY"
echo "=========================================="
echo "OLD (pre-optimization): ${OLD_ELAPSED}s"
echo "NEW (optimized):        ${NEW_ELAPSED}s"
if [ "$OLD_ELAPSED" -gt 0 ]; then
    SPEEDUP=$(echo "scale=2; $OLD_ELAPSED / $NEW_ELAPSED" | bc 2>/dev/null || echo "N/A")
    REDUCTION=$(echo "scale=1; (1 - $NEW_ELAPSED / $OLD_ELAPSED) * 100" | bc 2>/dev/null || echo "N/A")
    echo "Speedup:                ${SPEEDUP}x"
    echo "Time reduction:         ${REDUCTION}%"
fi
echo "=========================================="
echo ""
echo "Full results saved to: $BENCH_DIR"

# --- Cleanup worktree ---
echo "Cleaning up old code worktree..."
git -C "$REPO_ROOT" worktree remove "$OLD_WORKTREE" --force 2>/dev/null || true

echo "Done."
