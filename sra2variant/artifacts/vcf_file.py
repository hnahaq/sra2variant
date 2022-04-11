# import os
# import pathlib

# from sra2variant.artifacts.base_file import _FileArtifacts


# def init_vcf_file(
#     vcf_file: pathlib.Path,
#     output_dir: str,
#     csv_dir: str,
# ) -> _FileArtifacts:
#     vcf_id = os.path.basename(sra_file)
#     return _FileArtifacts(
#         sra_file,
#         working_id=sra_id,
#         cwd=os.path.join(output_dir, sra_id),
#         res_dir=csv_dir
#     )
