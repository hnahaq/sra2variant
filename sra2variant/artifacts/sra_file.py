import os

from sra2variant.artifacts.base_file import _FileArtifacts


def init_sra_file(
    sra_file: str,
    output_dir: str,
    csv_dir: str,
) -> _FileArtifacts:
    sra_id = os.path.basename(sra_file)
    if sra_file.endswith(".sra"):
        (sra_id, _) = os.path.splitext(sra_id)
    return _FileArtifacts(
        sra_file,
        working_id=sra_id,
        cwd=os.path.join(output_dir, sra_id),
        res_dir=csv_dir
    )
