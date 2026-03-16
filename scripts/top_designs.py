"""
Extract the top N designs from a generated_models_info.csv, ranked by binding_dG_with_waters.

Usage:
    python scripts/top_designs.py <generated_models_info.csv> <model_PDB_files_dir> [N]

    N defaults to 10.

Output:
    Writes top_designs.tsv to the same directory as the input CSV.
"""

import sys
from pathlib import Path

import pandas as pd


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <generated_models_info.csv> <model_PDB_files_dir> [N]")
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    pdb_dir = Path(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 10

    df = pd.read_csv(csv_path)

    # Only keep accepted mutations
    if 'mutation_accepted' in df.columns:
        df = df[df['mutation_accepted'].astype(str).str.strip().isin(['True', 'true', '1'])]

    # Sort by binding_dG_with_waters (most negative = best binding)
    if 'binding_dG_with_waters' in df.columns:
        sort_col = 'binding_dG_with_waters'
    else:
        sort_col = 'binding_dG'

    df = df.sort_values(sort_col, ascending=True).head(n)

    # Build output columns
    output_cols = [
        'backbone_PDB_file_name',
        'nth_model',
        'nth_iteration',
    ]
    if 'binding_dG_with_waters' in df.columns:
        output_cols.append('binding_dG_with_waters')
    output_cols.extend([
        'antibody_stability_dG',
        'complex_stability_dG',
        'binding_dG',
        'from_mut_name',
        'residue_IDs',
    ])

    result = df[output_cols].copy()

    # Add pdb_file column
    result['pdb_file'] = result.apply(
        lambda row: str(pdb_dir / f"{row['backbone_PDB_file_name']}_{row['nth_model']}_{row['nth_iteration']}.pdb"),
        axis=1,
    )

    output_path = csv_path.parent / 'top_designs.tsv'
    result.to_csv(output_path, sep='\t', index=False)
    print(f"Wrote {len(result)} designs to {output_path}")


if __name__ == '__main__':
    main()
