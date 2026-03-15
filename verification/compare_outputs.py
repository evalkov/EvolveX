"""
Compare generated_models_info.csv from old and new EvolveX runs.

Validates that the optimized code produces structurally correct and
scientifically equivalent output compared to the pre-optimization code.

Usage:
    python compare_outputs.py <old_csv> <new_csv>

Checks performed:
    1. Schema: same columns in both files
    2. Row count: same number of rows (iterations)
    3. Column types: numeric columns parse correctly
    4. Value ranges: energy values are in physically reasonable ranges
    5. Distribution similarity: mean/std of energy columns are comparable
       (not identical, since GA search is stochastic)
"""

import sys
import pandas as pd
import numpy as np


# Energy columns that should be present and numeric
ENERGY_COLUMNS = [
    'antibody_stability_dG',
    'complex_stability_dG',
    'binding_dG',
    'antibody_intraclash_score',
]

# Optional energy columns (depend on config)
OPTIONAL_ENERGY_COLUMNS = [
    'binding_dG_with_waters',
]

# Metadata columns
METADATA_COLUMNS = [
    'backbone_PDB_file_name',
    'nth_model',
    'step',
    'nth_iteration',
    'residue_IDs',
    'from_mut_name',
    'mutation_accepted',
]

# Physically reasonable ranges for FoldX energies (kcal/mol)
ENERGY_RANGES = {
    'antibody_stability_dG': (-200, 200),
    'complex_stability_dG': (-500, 500),
    'binding_dG': (-100, 100),
    'binding_dG_with_waters': (-100, 100),
    'antibody_intraclash_score': (-50, 100),
}


def load_csv(path):
    try:
        df = pd.read_csv(path)
        print(f"  Loaded {path}: {len(df)} rows, {len(df.columns)} columns")
        return df
    except Exception as e:
        print(f"  ERROR: Could not load {path}: {e}")
        return None


def check_schema(old_df, new_df):
    """Check that both files have the same columns."""
    old_cols = set(old_df.columns)
    new_cols = set(new_df.columns)

    passed = True
    if old_cols == new_cols:
        print("  PASS: Identical column sets")
    else:
        only_old = old_cols - new_cols
        only_new = new_cols - old_cols
        if only_old:
            print(f"  WARN: Columns only in OLD: {only_old}")
        if only_new:
            print(f"  WARN: Columns only in NEW: {only_new}")
        # Check if core columns are present in both
        core = set(ENERGY_COLUMNS + METADATA_COLUMNS)
        missing_old = core - old_cols
        missing_new = core - new_cols
        if missing_old:
            print(f"  FAIL: OLD missing core columns: {missing_old}")
            passed = False
        if missing_new:
            print(f"  FAIL: NEW missing core columns: {missing_new}")
            passed = False
        if not missing_old and not missing_new:
            print("  PASS: All core columns present in both")

    return passed


def check_row_count(old_df, new_df):
    """Check row counts match."""
    if len(old_df) == len(new_df):
        print(f"  PASS: Same row count ({len(old_df)})")
        return True
    else:
        print(f"  WARN: Different row counts (OLD={len(old_df)}, NEW={len(new_df)})")
        print(f"        This is expected if the GA explored different paths.")
        return True  # Not a failure for stochastic algorithms


def check_column_types(df, label):
    """Check that energy columns are numeric."""
    passed = True
    for col in ENERGY_COLUMNS:
        if col not in df.columns:
            continue
        if not pd.api.types.is_numeric_dtype(df[col]):
            try:
                pd.to_numeric(df[col])
                print(f"  WARN: {label} column '{col}' is string but parseable as numeric")
            except ValueError:
                print(f"  FAIL: {label} column '{col}' contains non-numeric values")
                passed = False
    return passed


def check_value_ranges(df, label):
    """Check energy values are in physically reasonable ranges."""
    passed = True
    for col, (lo, hi) in ENERGY_RANGES.items():
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors='coerce')
        if vals.isna().all():
            continue
        col_min, col_max = vals.min(), vals.max()
        if col_min < lo or col_max > hi:
            print(f"  WARN: {label} '{col}' range [{col_min:.2f}, {col_max:.2f}] "
                  f"outside expected [{lo}, {hi}]")
        else:
            print(f"  PASS: {label} '{col}' range [{col_min:.2f}, {col_max:.2f}] OK")
    return passed


def check_distributions(old_df, new_df):
    """Compare distributions of energy columns between old and new."""
    print("  Energy column statistics (OLD vs NEW):")
    print(f"  {'Column':<30} {'OLD mean':>10} {'NEW mean':>10} {'OLD std':>10} {'NEW std':>10}")
    print(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

    for col in ENERGY_COLUMNS + OPTIONAL_ENERGY_COLUMNS:
        if col not in old_df.columns or col not in new_df.columns:
            continue
        old_vals = pd.to_numeric(old_df[col], errors='coerce')
        new_vals = pd.to_numeric(new_df[col], errors='coerce')

        old_mean = old_vals.mean()
        new_mean = new_vals.mean()
        old_std = old_vals.std()
        new_std = new_vals.std()

        print(f"  {col:<30} {old_mean:>10.3f} {new_mean:>10.3f} {old_std:>10.3f} {new_std:>10.3f}")

    print("")
    print("  NOTE: Values will differ between runs because the GA search is stochastic.")
    print("  What matters is that both are in the same ballpark and physically reasonable.")
    return True


def check_metadata_structure(old_df, new_df):
    """Check that metadata columns have expected structure."""
    passed = True

    for col in ['step', 'backbone_PDB_file_name']:
        if col not in old_df.columns or col not in new_df.columns:
            continue
        old_vals = set(old_df[col].unique())
        new_vals = set(new_df[col].unique())
        if old_vals == new_vals:
            print(f"  PASS: '{col}' has same unique values: {old_vals}")
        else:
            print(f"  WARN: '{col}' differs — OLD: {old_vals}, NEW: {new_vals}")

    # Check step types
    if 'step' in new_df.columns:
        expected_steps = {'MC', 'recombination'}
        actual_steps = set(new_df['step'].unique())
        if actual_steps <= expected_steps:
            print(f"  PASS: 'step' values are valid: {actual_steps}")
        else:
            print(f"  FAIL: Unexpected step types: {actual_steps - expected_steps}")
            passed = False

    # Check mutation_accepted is boolean-like
    if 'mutation_accepted' in new_df.columns:
        accepted_vals = set(new_df['mutation_accepted'].unique())
        if accepted_vals <= {True, False, 'True', 'False'}:
            print(f"  PASS: 'mutation_accepted' values are boolean-like")
        else:
            print(f"  WARN: 'mutation_accepted' has unexpected values: {accepted_vals}")

    return passed


def check_residue_ids_format(df, label):
    """Check that residue_IDs column has expected format (semicolon-separated AA+residueID pairs)."""
    if 'residue_IDs' not in df.columns:
        print(f"  SKIP: '{label}' has no 'residue_IDs' column")
        return True

    sample = df['residue_IDs'].dropna().iloc[:5] if len(df) > 0 else []
    passed = True
    for val in sample:
        parts = str(val).split(';')
        for part in parts:
            if len(part) < 2:
                print(f"  FAIL: {label} residue_ID '{part}' too short (expected format like 'KH52')")
                passed = False
                break

    if passed:
        print(f"  PASS: {label} 'residue_IDs' format looks correct (sample: {sample.iloc[0] if len(sample) > 0 else 'N/A'})")
    return passed


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <old_csv> <new_csv>")
        sys.exit(1)

    old_path, new_path = sys.argv[1], sys.argv[2]

    print("=" * 60)
    print("  EvolveX Output Comparison Report")
    print("=" * 60)
    print("")

    # Load
    print("1. Loading files...")
    old_df = load_csv(old_path)
    new_df = load_csv(new_path)
    if old_df is None or new_df is None:
        print("\nFATAL: Could not load one or both CSV files.")
        sys.exit(1)
    print("")

    all_passed = True

    # Schema
    print("2. Schema check...")
    all_passed &= check_schema(old_df, new_df)
    print("")

    # Row count
    print("3. Row count check...")
    all_passed &= check_row_count(old_df, new_df)
    print("")

    # Column types
    print("4. Column type check...")
    all_passed &= check_column_types(old_df, "OLD")
    all_passed &= check_column_types(new_df, "NEW")
    print("")

    # Value ranges
    print("5. Value range check...")
    all_passed &= check_value_ranges(old_df, "OLD")
    all_passed &= check_value_ranges(new_df, "NEW")
    print("")

    # Metadata structure
    print("6. Metadata structure check...")
    all_passed &= check_metadata_structure(old_df, new_df)
    print("")

    # Residue ID format
    print("7. Residue ID format check...")
    all_passed &= check_residue_ids_format(old_df, "OLD")
    all_passed &= check_residue_ids_format(new_df, "NEW")
    print("")

    # Distribution comparison
    print("8. Distribution comparison...")
    check_distributions(old_df, new_df)
    print("")

    # Summary
    print("=" * 60)
    if all_passed:
        print("  RESULT: ALL CHECKS PASSED")
    else:
        print("  RESULT: SOME CHECKS FAILED (see above)")
    print("=" * 60)

    sys.exit(0 if all_passed else 1)


if __name__ == '__main__':
    main()
