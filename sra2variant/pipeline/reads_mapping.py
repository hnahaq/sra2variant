import os

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class BWAindexWrapper(CMDwrapperBase):

    exec_name = "bwa index"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (refSeq_file, ) = input_files.file_path
        self.output_files = None
        super().__init__(
            input_files,
            *args,
            "-p", os.path.abspath(refSeq_file),
            os.path.abspath(refSeq_file)
        )

    def _post_execution(self) -> None:
        self.output_files = self.input_files


class BWAmemWrapper(CMDwrapperBase):
    
    exec_name = "bwa mem"
    bwa_index: BWAindexWrapper = None

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        self.output_files = input_files.coupled_files(
            ".sam",
            exec_name=self.exec_name
        )
        (sam_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = self.bwa_index.output_files.file_path
        input_flags = input_files.path_from_cwd()
        super().__init__(
            input_files,
            "-t", self.threads,
            *args,
            "-o", sam_file,
            os.path.abspath(refSeq_file),
            *input_flags,
        )

    @classmethod
    def set_base_index(cls, bwa_index: BWAindexWrapper) -> None:
        if bwa_index.output_files is None:
            bwa_index.execute_cmd()
        cls.bwa_index = bwa_index

    def _post_execution(self) -> None:
        pass
