import os

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class FastqDumpWrapper(CMDwrapperBase):

    exec_name: str = "fastq-dump"

    def __init__(
        self,
        input_files: _FileArtifacts,
    ) -> None:
        (sra_file, ) = input_files.path_from_cwd()
        super().__init__(
            input_files,
            *self.exec_args,
            "--skip-technical",
            "--split-3",
            "--outdir", input_files.cwd,
            sra_file
        )

    def post_execution(self) -> None:
        paired_files = self.input_files.coupled_files(
            ("_1.fastq", "_2.fastq")
        )
        unpaired_file = self.input_files.coupled_files(
            (".fastq", )
        )
        if paired_files.exist():
            self.output_files = paired_files
        elif unpaired_file.exist():
            self.output_files = unpaired_file


class PrefetchWrapper(CMDwrapperBase):

    exec_name: str = "prefetch"

    def __init__(self, sra_id: str, working_dir: str = None) -> None:
        if working_dir is None:
            working_dir = "."
        sra_file = os.path.join(working_dir, f"{sra_id}.sra")
        self.output_files = _FileArtifacts(
            file_path=(sra_file, ),
            cwd=working_dir,
            working_id=sra_id,
        )
        input_files = _FileArtifacts(
            file_path=(".", ),
            cwd=".",
            working_id="placeholder"
        )
        super().__init__(
            input_files,
            *self.exec_args,
            "--output-file", sra_file,
            sra_id
        )

    def post_execution(self) -> None:
        pass
