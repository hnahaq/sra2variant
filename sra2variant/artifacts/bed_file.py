import os

from sra2variant.artifacts.base_file import _FileArtifacts


def init_bed_file(file_path: str):
    bed_file = _FileArtifacts(
        file_path,
        cwd=os.path.dirname(file_path),
    )
    if not bed_file.exist():
        raise FileNotFoundError(f"Primer file {file_path} not found")
    return bed_file
