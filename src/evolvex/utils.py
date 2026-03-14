
import shutil
import tarfile

NDIGIS_ROUNDING = 4

def save_compressed_PDB_file(PDB_file_path, output_name, output_dir):
    with tarfile.open(output_dir / f'{output_name}.tar.gz', 'w:gz') as tar_file_handle:
        tar_file_handle.add(PDB_file_path, arcname=output_name)
    return

def save_PDB_file_copy(PDB_file_path, output_name, output_dir):
    """Copy a PDB file to output_dir without compression. Use batch_compress_PDB_files() after search completes."""
    shutil.copy(PDB_file_path, output_dir / output_name)

def batch_compress_PDB_files(pdb_dir):
    """Compress all PDB files in a directory to tar.gz and remove the originals."""
    for pdb_file in pdb_dir.glob('*.pdb'):
        with tarfile.open(pdb_dir / f'{pdb_file.name}.tar.gz', 'w:gz') as tar:
            tar.add(pdb_file, arcname=pdb_file.name)
        pdb_file.unlink()
