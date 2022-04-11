import os

from sra2variant.artifacts.base_file import _FileArtifacts


def init_tsv_file(file_path: str):
    tsv_file = _FileArtifacts(
        file_path,
        cwd=os.path.dirname(file_path)
    )
    if not tsv_file.exist():
        raise FileNotFoundError(f"Primer file {file_path} not found")
    return tsv_file
