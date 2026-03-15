"""
Validate that Raw_/Dif_ file parsing produces the same complex stability values
as the explicit FoldX Stability command.

This is the critical correctness check for Phase 1A: we eliminated the Stability
command and now read total energy from BuildModel's Raw_ file instead.

Usage:
    python validate_raw_files.py <evolvex_working_dir> <foldx_dir>

What it does:
    1. Finds a model directory that still has BuildModel output files (Raw_, Dif_)
       If none exist (clean_up_model_dir removes them), it creates a fresh model
       from the Vsig4 example to generate the files.
    2. Reads the mutant total energy from Raw_ file
    3. Reads the ddG from Dif_ file
    4. Computes wildtype energy as (mutant - ddG)
    5. Runs FoldX Stability on both mutant and wildtype PDBs
    6. Compares the Stability _ST.fxout values with the Raw_/Dif_ derived values
    7. Reports pass/fail with tolerance
"""

import subprocess
import sys
from pathlib import Path


TOLERANCE = 0.01  # kcal/mol tolerance for floating point comparison
NDIGIS_ROUNDING = 4


def find_or_create_test_model(evolvex_working_dir, foldx_dir):
    """Find an existing model dir with Raw_/Dif_ files, or create one."""
    evolvex_dir = Path(evolvex_working_dir) / 'EvolveX'

    if not evolvex_dir.exists():
        print(f"  ERROR: {evolvex_dir} does not exist. Run EvolveX first.")
        return None

    # Look for any PDB backbone directory
    for pdb_dir in evolvex_dir.iterdir():
        if not pdb_dir.is_dir():
            continue

        alanine_pdb = list(pdb_dir.glob('*_1_Alanine_mutant.pdb'))
        if not alanine_pdb:
            continue

        alanine_pdb = alanine_pdb[0]
        pdb_name = pdb_dir.name

        # Create a temporary test model using BuildModel
        test_dir = pdb_dir / '_raw_validation_test'
        test_dir.mkdir(exist_ok=True)

        import shutil
        input_pdb = shutil.copy(alanine_pdb, test_dir)
        input_pdb = Path(input_pdb)

        # Get a valid mutation from the positions to explore
        # Use a simple self-mutation (A->A at first Ala position) to keep it simple
        # Actually, let's read the sequence and make a real mutation
        search_results = pdb_dir / 'search_results'
        mutations_summary = pdb_dir / 'hotspot_mutants' / 'all_mutations_summary.csv'

        if mutations_summary.exists():
            import pandas as pd
            df = pd.read_csv(mutations_summary, header=0, index_col=0)
            # Pick the first mutation
            mut_name = df.index[0]
            print(f"  Using test mutation: {mut_name}")
        else:
            print(f"  ERROR: No mutations summary found at {mutations_summary}")
            shutil.rmtree(test_dir)
            continue

        # Write individual_list file
        mutations_file = test_dir / 'individual_list_foldx_mutations_file.txt'
        with open(mutations_file, 'w') as f:
            f.write(f'{mut_name};\n')

        # Run BuildModel
        print(f"  Running BuildModel with mutation {mut_name}...")
        cmd = [
            str(Path(foldx_dir) / 'foldx'),
            '--command', 'BuildModel',
            '--pdb-dir', str(test_dir),
            '--pdb', input_pdb.name,
            '--mutant-file', str(mutations_file),
            '--pdbHydrogens', 'true',
            '--output-dir', str(test_dir),
            '--moveNeighbours', 'true',
            '--vdwDesign', '2',
            '--screen', 'false'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: BuildModel failed: {result.stderr[:500]}")
            shutil.rmtree(test_dir)
            continue

        # Check that Raw_ and Dif_ files were generated
        raw_file = test_dir / f'Raw_{input_pdb.stem}.fxout'
        dif_file = test_dir / f'Dif_{input_pdb.stem}.fxout'
        mutant_pdb = test_dir / f'{input_pdb.stem}_1.pdb'
        wildtype_pdb = test_dir / f'WT_{input_pdb.stem}_1.pdb'

        if not raw_file.exists():
            print(f"  ERROR: Raw_ file not generated at {raw_file}")
            shutil.rmtree(test_dir)
            continue

        if not mutant_pdb.exists():
            # wildtype copy
            shutil.copy(input_pdb, wildtype_pdb)

        return {
            'test_dir': test_dir,
            'raw_file': raw_file,
            'dif_file': dif_file,
            'mutant_pdb': mutant_pdb,
            'wildtype_pdb': wildtype_pdb,
            'pdb_stem': input_pdb.stem,
        }

    print("  ERROR: Could not find any suitable PDB backbone directory")
    return None


def read_raw_energy(raw_file_path):
    """Read total energy from BuildModel Raw_ file."""
    with open(raw_file_path, 'rt') as f:
        lines = f.readlines()
    return round(float(lines[9].split('\t')[1]), NDIGIS_ROUNDING)


def read_dif_ddg(dif_file_path):
    """Read stability ddG from BuildModel Dif_ file."""
    with open(dif_file_path, 'rt') as f:
        lines = f.readlines()
    return round(float(lines[9].split('\t')[1]), NDIGIS_ROUNDING)


def run_stability_and_read(foldx_dir, pdb_dir, pdb_stem, output_dir):
    """Run FoldX Stability command and read the _ST.fxout result."""
    cmd = [
        str(Path(foldx_dir) / 'foldx'),
        '--command', 'Stability',
        '--pdb-dir', str(pdb_dir),
        '--pdb', f'{pdb_stem}.pdb',
        '--vdwDesign', '2',
        '--output-dir', str(output_dir),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: Stability command failed for {pdb_stem}: {result.stderr[:500]}")
        return None

    st_file = output_dir / f'{pdb_stem}_0_ST.fxout'
    if not st_file.exists():
        print(f"  ERROR: ST file not found at {st_file}")
        return None

    with open(st_file, 'rt') as f:
        lines = f.readlines()
    return round(float(lines[0].split('\t')[1]), NDIGIS_ROUNDING)


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <evolvex_working_dir> <foldx_dir>")
        sys.exit(1)

    evolvex_working_dir = sys.argv[1]
    foldx_dir = sys.argv[2]

    print("=" * 60)
    print("  Raw_ File Validation Report")
    print("=" * 60)
    print("")

    # Step 1: Find or create test model
    print("1. Setting up test model...")
    test_info = find_or_create_test_model(evolvex_working_dir, foldx_dir)
    if test_info is None:
        print("\nFATAL: Could not set up test model. Ensure EvolveX has been run first.")
        sys.exit(1)
    print("")

    test_dir = test_info['test_dir']
    raw_file = test_info['raw_file']
    dif_file = test_info['dif_file']
    mutant_pdb = test_info['mutant_pdb']
    wildtype_pdb = test_info['wildtype_pdb']

    # Step 2: Read Raw_ and Dif_ values
    print("2. Reading BuildModel output files...")
    mutant_energy_from_raw = read_raw_energy(raw_file)
    ddg_from_dif = read_dif_ddg(dif_file)
    wildtype_energy_from_raw = round(mutant_energy_from_raw - ddg_from_dif, NDIGIS_ROUNDING)
    print(f"  Raw_ mutant energy:   {mutant_energy_from_raw}")
    print(f"  Dif_ ddG:             {ddg_from_dif}")
    print(f"  Derived WT energy:    {wildtype_energy_from_raw}")
    print("")

    # Step 3: Run Stability command
    print("3. Running FoldX Stability for comparison...")
    mutant_energy_from_stability = run_stability_and_read(
        foldx_dir, mutant_pdb.parent, mutant_pdb.stem, test_dir
    )
    wildtype_energy_from_stability = run_stability_and_read(
        foldx_dir, wildtype_pdb.parent, wildtype_pdb.stem, test_dir
    )

    if mutant_energy_from_stability is None or wildtype_energy_from_stability is None:
        print("\nFATAL: Stability command failed.")
        cleanup(test_dir)
        sys.exit(1)

    print(f"  Stability mutant energy:  {mutant_energy_from_stability}")
    print(f"  Stability wildtype energy: {wildtype_energy_from_stability}")
    print("")

    # Step 4: Compare
    print("4. Comparing values...")
    mutant_diff = abs(mutant_energy_from_raw - mutant_energy_from_stability)
    wildtype_diff = abs(wildtype_energy_from_raw - wildtype_energy_from_stability)

    print(f"  Mutant energy difference:   {mutant_diff:.4f} kcal/mol (tolerance: {TOLERANCE})")
    print(f"  Wildtype energy difference: {wildtype_diff:.4f} kcal/mol (tolerance: {TOLERANCE})")
    print("")

    passed = True
    if mutant_diff <= TOLERANCE:
        print(f"  PASS: Mutant energy matches (Raw_={mutant_energy_from_raw}, "
              f"Stability={mutant_energy_from_stability})")
    else:
        print(f"  FAIL: Mutant energy mismatch (Raw_={mutant_energy_from_raw}, "
              f"Stability={mutant_energy_from_stability}, diff={mutant_diff})")
        passed = False

    if wildtype_diff <= TOLERANCE:
        print(f"  PASS: Wildtype energy matches (derived={wildtype_energy_from_raw}, "
              f"Stability={wildtype_energy_from_stability})")
    else:
        print(f"  FAIL: Wildtype energy mismatch (derived={wildtype_energy_from_raw}, "
              f"Stability={wildtype_energy_from_stability}, diff={wildtype_diff})")
        passed = False

    print("")

    # Cleanup
    print("5. Cleaning up test files...")
    cleanup(test_dir)
    print("")

    # Summary
    print("=" * 60)
    if passed:
        print("  RESULT: Raw_ FILE PARSING VALIDATED SUCCESSFULLY")
        print("  The Raw_/Dif_ approach produces equivalent results to Stability.")
    else:
        print("  RESULT: VALIDATION FAILED")
        print("  The Raw_/Dif_ values do not match Stability output.")
        print("  Investigate before using the optimized code in production.")
    print("=" * 60)

    sys.exit(0 if passed else 1)


def cleanup(test_dir):
    """Remove the temporary test directory."""
    import shutil
    try:
        shutil.rmtree(test_dir)
        print(f"  Removed {test_dir}")
    except Exception as e:
        print(f"  WARN: Could not clean up {test_dir}: {e}")


if __name__ == '__main__':
    main()
