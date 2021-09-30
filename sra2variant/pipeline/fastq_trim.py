from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class FastpWrapper(CMDwrapperBase):

    exec_name: str = "fastp"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        if len(input_files) == 2:
            self.output_files = input_files.coupled_files(
                "_1.fastq.gz",
                "_2.fastq.gz",
                exec_name=self.exec_name
            )
            (paired1_file, paired2_file) = input_files.path_from_cwd()
            (trimmed1_file, trimmed2_file) = self.output_files.path_from_cwd()
            io_files_flags = (
                "-i", paired1_file,
                "-I", paired2_file,
                "-o", trimmed1_file,
                "-O", trimmed2_file,
            )
        elif len(input_files) == 1:
            self.output_files = input_files.coupled_files(
                ".fastq.gz",
                exec_name=self.exec_name
            )
            (unpaired_file, ) = input_files.path_from_cwd()
            (trimmed_file, ) = self.output_files.path_from_cwd()
            io_files_flags = (
                "-i", unpaired_file,
                "-o", trimmed_file,
            )
        super().__init__(
            input_files,
            "--thread", self.threads,
            *args,
            *io_files_flags
        )

    def _post_execution(self) -> None:
        pass


class TrimmomaticWrapper(CMDwrapperBase):

    exec_name: str = "trimmomatic"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        if len(input_files) == 2:
            self.output_files = input_files.coupled_files(
                "_1.fastq.gz",
                "_2.fastq.gz",
                exec_name=self.exec_name
            )
            (paired1_file, paired2_file) = input_files.path_from_cwd()
            (trimmed1_file, trimmed2_file) = self.output_files.path_from_cwd()
            io_files_flags = (
                "PE",
                paired1_file,
                paired2_file,
                trimmed1_file,
                f"{input_files.workding_id}_1_unpaired.fastq.gz",
                trimmed2_file,
                f"{input_files.workding_id}_2_unpaired.fastq.gz",
            )
        elif len(input_files) == 1:
            self.output_files = input_files.coupled_files(
                ".fastq.gz",
                exec_name=self.exec_name
            )
            (unpaired_file, ) = input_files.path_from_cwd()
            (trimmed_file, ) = self.output_files.path_from_cwd()
            io_files_flags = (
                "SE",
                unpaired_file,
                trimmed_file,
            )
        super().__init__(
            input_files,
            *io_files_flags,
            "-summary", f"{input_files.workding_id}_trimmomatic.summary",
            "-threads", self.threads,
            "-phred33",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:3:15",
            "MINLEN:36",
            *args,
        )

    def _post_execution(self) -> None:
        pass
