import os
import pysam

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class BWAindexWrapper(CMDwrapperBase):

    exec_name = "bwa index"

    def __init__(self, input_files: _FileArtifacts) -> None:
        (refSeq_file, ) = input_files.file_path
        self.output_files = None
        super().__init__(
            input_files,
            *self.exec_args,
            "-p", os.path.abspath(refSeq_file),
            os.path.abspath(refSeq_file)
        )

    def post_execution(self) -> None:
        self.output_files = self.input_files


class BWAmemWrapper(CMDwrapperBase):
    
    exec_name = "bwa mem"
    bwa_index: BWAindexWrapper = None

    def __init__(self, input_files: _FileArtifacts) -> None:
        self.output_files = input_files.coupled_files(("_bwa-mem.sam", ))
        (sam_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = self.bwa_index.output_files.file_path
        input_flags = input_files.path_from_cwd()
        super().__init__(
            input_files,
            "-t", self.threads,
            *self.exec_args,
            "-o", sam_file,
            os.path.abspath(refSeq_file),
            *input_flags,
        )

    @classmethod
    def set_base_index(cls, bwa_index: BWAindexWrapper) -> None:
        if bwa_index.output_files is None:
            bwa_index.execute_cmd()
        cls.bwa_index = bwa_index

    def post_execution(self) -> None:
        (sam_file, ) = self.output_files.file_path
        self.output_files = self.output_files.coupled_files(
            ("_bwa-mem.sort.bam", )
        )
        (sort_bam_file, ) = self.output_files.file_path
        print(f"Sort {sam_file}")
        pysam.sort("-@", self.threads, "-o", sort_bam_file, sam_file)
        pysam.index(sort_bam_file)
        if len(self.input_files) == 2:
            self.output_files = self.output_files.coupled_files(
                ("_bwa-mem.sort.trimmed.bam", )
            )
            (sort_trimmed_bam_file, ) = self.output_files.file_path
            pysam.view(
                "-b",
                "-@", self.threads,
                "-q", "20",
                "-f", "1",
                "-F", "268",
                "-o", sort_trimmed_bam_file,
                sort_bam_file,
                catch_stdout = False
            )


class Hisat2BuildWrapper(CMDwrapperBase):

    exec_name = "hisat2-build"

    __slots__ = ["ht2_index_base"]

    def __init__(self, input_files: _FileArtifacts) -> None:
        (refSeq_file, ) = input_files.file_path
        self.ht2_index_base = None
        self.output_files = None
        super().__init__(
            input_files,
            *self.exec_args,
            "-p", self.threads,
            os.path.abspath(refSeq_file),
            os.path.abspath(refSeq_file)
        )

    def post_execution(self) -> None:
        (refSeq_file, ) = self.input_files.file_path
        self.ht2_index_base = os.path.abspath(refSeq_file)


class Hisat2Wrapper(CMDwrapperBase):

    exec_name: str = "hisat2"
    ht2_index_base: str = None

    def __init__(
        self,
        input_files: _FileArtifacts,
    ) -> None:
        self.output_files = input_files.coupled_files(("_hisat2.sam", ))
        (sam_file, ) = self.output_files.path_from_cwd()
        if len(input_files) == 2:
            (paired1_file, paired2_file) = input_files.path_from_cwd()
            input_flags = ("-1", paired1_file, "-2", paired2_file)
        elif len(input_files) == 1:
            (unpaired_file, ) = input_files.path_from_cwd()
            input_flags = ("-U", unpaired_file, )
        super().__init__(
            input_files,
            "--threads", self.threads,
            *self.exec_args,
            "-x", input_files.relative_path(self.ht2_index_base),
            *input_flags,
            "-S", sam_file
        )

    @classmethod
    def set_base_index(cls, ht2_build: Hisat2BuildWrapper) -> None:
        if ht2_build.ht2_index_base is None:
            ht2_build.execute_cmd()
        cls.ht2_index_base = ht2_build.ht2_index_base

    def post_execution(self) -> None:
        (sam_file, ) = self.output_files.file_path
        self.output_files = self.output_files.coupled_files(
            ("_hisat2.sort.bam", )
        )
        (sort_bam_file, ) = self.output_files.file_path
        print(f"Sort {sam_file}")
        pysam.sort("-@", self.threads, "-o", sort_bam_file, sam_file)
        pysam.index(sort_bam_file)


class Bowtie2BuildWrapper(CMDwrapperBase):

    exec_name: str = "bowtie2-build"

    __slots__ = ["bt2_index_base"]

    def __init__(self, input_files: _FileArtifacts) -> None:
        (refSeq_file, ) = input_files.file_path
        self.bt2_index_base = None
        self.output_files = None
        super().__init__(
            input_files,
            "--threads", self.threads,
            *self.exec_args,
            os.path.abspath(refSeq_file),
            os.path.abspath(refSeq_file)
        )

    def post_execution(self) -> None:
        (refSeq_file, ) = self.input_files.file_path
        self.bt2_index_base = os.path.abspath(refSeq_file)


class Bowtie2Wrapper(CMDwrapperBase):

    exec_name: str = "bowtie2"
    bt2_index_base: str = None

    def __init__(
        self,
        input_files: _FileArtifacts,
    ) -> None:
        self.output_files = input_files.coupled_files(("_bowtie2.sam", ))
        (sam_file, ) = self.output_files.path_from_cwd()
        if len(input_files) == 2:
            (paired1_file, paired2_file) = input_files.path_from_cwd()
            input_flags = ("-1", paired1_file, "-2", paired2_file)
        elif len(input_files) == 1:
            (unpaired_file, ) = input_files.path_from_cwd()
            input_flags = ("-U", unpaired_file, )
        super().__init__(
            input_files,
            "--threads", self.threads,
            *self.exec_args,
            "-x", input_files.relative_path(self.bt2_index_base),
            *input_flags,
            "-S", sam_file
        )

    @classmethod
    def set_base_index(cls, bt2_build: Bowtie2BuildWrapper) -> None:
        if bt2_build.bt2_index_base is None:
            bt2_build.execute_cmd()
        cls.bt2_index_base = bt2_build.bt2_index_base

    def post_execution(self) -> None:
        (sam_file, ) = self.output_files.file_path
        
        # if len(self.input_files) == 2:
        #     self.output_files = self.output_files.coupled_files(
        #         ("_bowtie2.fixmate.bam", )
        #     )
        #     (fixmate_bam_file, ) = self.output_files.file_path
        #     print(f"Fixmate {fixmate_bam_file}")
        #     pysam.fixmate("-@", self.threads, "-o", fixmate_bam_file, sam_file)
        #     # The input for samtools sort changed
        #     sam_file = fixmate_bam_file
        self.output_files = self.output_files.coupled_files(
            ("_bowtie2.sort.bam", )
        )
        (sort_bam_file, ) = self.output_files.file_path
        print(f"Sort {sam_file}")
        pysam.sort("-@", self.threads, "-o", sort_bam_file, sam_file)
        pysam.index(sort_bam_file)
