import os

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class FastqDumpWrapper(CMDwrapperBase):

    exec_name: str = "fastq-dump"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (sra_file, ) = input_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "--skip-technical",
            "--split-3",
            "--outdir", input_files.cwd,
            sra_file
        )

    def _post_execution(self) -> None:
        paired_files = self.input_files.coupled_files(
            "_1.fastq",
            "_2.fastq",
            exec_name=""
        )
        unpaired_file = self.input_files.coupled_files(
            ".fastq",
            exec_name=""
        )
        if paired_files.exist():
            self.output_files = paired_files
        elif unpaired_file.exist():
            self.output_files = unpaired_file


class PrefetchWrapper(CMDwrapperBase):

    exec_name: str = "prefetch"

    def __init__(self, sra_id: str, working_dir: str = None, *args: str) -> None:
        if working_dir is None:
            working_dir = "."
        sra_file = os.path.join(working_dir, f"{sra_id}.sra")
        self.output_files = _FileArtifacts(
            sra_file,
            cwd=working_dir,
            working_id=sra_id,
        )
        input_files = _FileArtifacts(working_id=sra_id, cwd=working_dir)
        super().__init__(
            input_files,
            *args,
            "--output-file", sra_file,
            sra_id
        )

    def _post_execution(self) -> None:
        pass
