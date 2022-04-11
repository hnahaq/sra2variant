import os

from sra2variant.artifacts.base_file import _FileArtifacts


def init_fasta_file(file_path: str):
    fasta_file = _FileArtifacts(
        file_path,
        cwd=os.path.dirname(file_path),
        working_id=os.path.basename(file_path)
    )
    if not fasta_file.exist():
        raise FileNotFoundError(f"Reference genome {file_path} not found")
    return fasta_file
