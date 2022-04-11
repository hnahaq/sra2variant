import os
import pathlib

from sra2variant.artifacts.base_file import _FileArtifacts


def init_vcf_file(
    vcf_file: pathlib.Path,
    csv_dir: str,
) -> _FileArtifacts:
    vcf_file: str = str(vcf_file)
    vcf_id = os.path.basename(vcf_file)
    (vcf_id, _) = os.path.splitext(vcf_id)
    return _FileArtifacts(
        vcf_file,
        working_id=vcf_id,
        cwd=os.path.join(csv_dir, vcf_id),
        res_dir=csv_dir
    )
