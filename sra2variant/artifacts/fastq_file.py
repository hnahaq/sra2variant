import os
import pathlib

from sra2variant.artifacts.base_file import _FileArtifacts


def init_fastq_pair(
    fastq_forward: pathlib.Path,
    forward_ext: str,
    reverse_ext: str,
    output_dir: str,
    csv_dir: str
) -> _FileArtifacts:
    fastq_forward = str(fastq_forward)
    fastq_id = fastq_forward[:-len(forward_ext)]
    fastq_reverse = fastq_id + reverse_ext
    fastq_id = os.path.basename(fastq_id)
    return _FileArtifacts(
        fastq_forward,
        fastq_reverse,
        working_id=fastq_id,
        cwd=os.path.join(output_dir, fastq_id),
        res_dir=csv_dir
    )
