from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class PicardMarkDuplicatedWrapper(CMDwrapperBase):

    exec_name: str = "picard MarkDuplicates"

    def __init__(self, input_files: _FileArtifacts) -> None:
        super().__init__(
            input_files,
            *self.exec_args
        )
        self.output_files = input_files.coupled_files(("_picard_MarkDuplicates.bam", ))
        (sort_bam_file, ) = input_files.path_from_cwd()
        (deduplicated_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            f"INPUT={sort_bam_file}",
            f"OUTPUT={deduplicated_file}",
            f"METRICS_FILE={input_files.workding_id}_picard_MarkDuplicates.metrics.txt",
            "REMOVE_DUPLICATES=true",
            "ASSUME_SORTED=true",
            "DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES",
            "OPTICAL_DUPLICATE_PIXEL_DISTANCE=100",
            "VALIDATION_STRINGENCY=LENIENT",
            *self.exec_args,
        )

    def post_execution(self) -> None:
        pass
