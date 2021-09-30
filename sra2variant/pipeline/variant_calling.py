import os

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


"""
gatk CreateSequenceDictionary -R NC_045512.2.fasta
gatk AddOrReplaceReadGroups I=SRR15432222.sort.bam O=SRR15432222_with_RG.sort.bam RGID=4 RGLB=lib1 RGPL=illumina RGSM=1 RGPU=unit1 CREATE_INDEX=True
gatk HaplotypeCaller -R NC_045512.2.fasta -I SRR15432222.sort.bam -O SRR15432222.vcf -ERC GVCF -ploidy 1
"""


class LoFreqFaidxWrapper(CMDwrapperBase):
    
    exec_name: str = "lofreq faidx"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        self.output_files = None
        (refSeq_file, ) = input_files.file_path
        super().__init__(
            input_files,
            *args,
            refSeq_file
        )

    def _post_execution(self) -> None:
        self.output_files = self.input_files


class LoFreqViterbiWraper(CMDwrapperBase):

    exec_name: str = "lofreq viterbi"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        (sort_bam_file, ) = input_files.path_from_cwd()
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (viterbi_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            *args,
            "--ref", os.path.abspath(refSeq_file),
            "--out", viterbi_file,
            sort_bam_file,
        )

    def _post_execution(self) -> None:
        pass


class LoFreqIndelQualWrapper(CMDwrapperBase):

    exec_name: str = "lofreq indelqual"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (sort_bam_file, ) = input_files.path_from_cwd()
        (indelqual_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            *args,
            "--ref", os.path.abspath(refSeq_file),
            "--out", indelqual_file,
            sort_bam_file,
        )

    def _post_execution(self) -> None:
        pass


class LoFreqCallWrapper(CMDwrapperBase):

    exec_name: str = "lofreq call-parallel"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.output_files = input_files.coupled_files(
            ".vcf",
            exec_name="lofreq call"
        )
        (indelqual_file, ) = input_files.path_from_cwd()
        (vcf_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--pp-threads", self.threads,
            *args,
            "--verbose",
            "--ref", os.path.abspath(refSeq_file),
            "--out", vcf_file,
            indelqual_file
        )

    def _post_execution(self) -> None:
        pass


class LoFreqFilterWrapper(CMDwrapperBase):

    exec_name: str = "lofreq filter"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        self.output_files = input_files.coupled_files(
            ".vcf",
            exec_name=self.exec_name
        )
        (vcf_file, ) = input_files.path_from_cwd()
        (filter_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-i", vcf_file,
            "-o", filter_file,
        )

    def _post_execution(self) -> None:
        pass


class BCFtoolsMpileupWrapper(CMDwrapperBase):

    exec_name: str = "bcftools mpileup"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.output_files = input_files.coupled_files("_bcftools_mpileup")
        (sort_bam_file, ) = input_files.path_from_cwd()
        (mpileup_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--max-idepth", "1000000",
            "--max-depth", "1000000",
            "--threads", self.threads,
            "--fasta-ref", os.path.abspath(refSeq_file),
            "--output-type", "b",
            *args,
            "--annotate", "INFO/AD,INFO/ADF,INFO/ADR",
            "--output", mpileup_file,
            sort_bam_file
        )

    def _post_execution(self) -> None:
        pass


# class BCFtoolsCallWrapper(CMDwrapperBase):

#     exec_name: str = "bcftools call"

#     def __init__(
#         self,
#         input_files: _FileArtifacts,
#         *args: str
#     ) -> None:
#         self.output_files = input_files.coupled_files("_bcftools_call")
#         (mpileup_file, ) = input_files.path_from_cwd()
#         (call_bcf_file, ) = self.output_files.path_from_cwd()
#         super().__init__(
#             input_files,
#             "--threads", self.threads,
#             "--output-type", "b",
#             "--prior", "0",
#             "--ploidy", "1",
#             *args,
#             "--keep-alts",
#             "--variants-only",
#             "--multiallelic-caller",
#             "--output", call_bcf_file,
#             mpileup_file
#         )

#     def _post_execution(self) -> None:
#         pass


class BCFtoolsNormWrapper(CMDwrapperBase):

    exec_name: str = "bcftools norm"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.output_files = input_files.coupled_files("_bcftools_norm")
        (call_bcf_file, ) = input_files.path_from_cwd()
        (norm_vcf_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--threads", self.threads,
            "--multiallelics", "-any",
            "--output-type", "v",
            "--fasta-ref", os.path.abspath(refSeq_file),
            *args,
            "--output", norm_vcf_file,
            call_bcf_file
        )

    def _post_execution(self) -> None:
        pass
