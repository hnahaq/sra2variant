from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class SAMtoolsSortWrapper(CMDwrapperBase):

    exec_name: str = "samtools sort"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (sam_file, ) = input_files.path_from_cwd()
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name="sort"
        )
        (sort_bam_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            "--threads", self.threads,
            *args,
            "--output-fmt", "bam",
            "-o", sort_bam_file,
            sam_file
        )

    def _post_execution(self) -> None:
        pass


class SAMtoolsIndexWrapper(CMDwrapperBase):

    exec_name: str = "samtools index"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (sam_file, ) = input_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-b",
            sam_file
        )

    def _post_execution(self) -> None:
        self.output_files = self.input_files


class SAMtoolsViewWrapper(CMDwrapperBase):

    exec_name: str = "samtools view"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (sam_file, ) = input_files.path_from_cwd()
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name="view"
        )
        (trimmed_sam_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            "--threads", self.threads,
            *args,
            "-b",
            "-o", trimmed_sam_file,
            sam_file
        )

    def _post_execution(self) -> None:
        pass


class PicardMarkDuplicatedWrapper(CMDwrapperBase):

    exec_name: str = "picard MarkDuplicates"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (sort_bam_file, ) = input_files.path_from_cwd()
        (deduplicated_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            f"INPUT={sort_bam_file}",
            f"OUTPUT={deduplicated_file}",
            f"METRICS_FILE={input_files.workding_id}_picard_MarkDuplicates.metrics.txt",
        )

    def _post_execution(self) -> None:
        pass
