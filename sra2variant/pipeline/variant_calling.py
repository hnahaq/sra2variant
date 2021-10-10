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
            exec_name="viterbi"
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
            exec_name="indelqual"
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


class IvarTrimWrapper(CMDwrapperBase):

    exec_name = "ivar trim"
    primer_bed_file: _FileArtifacts = None
    amplicon_info_file: _FileArtifacts = None

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (sort_bam_file, ) = input_files.path_from_cwd()
        bed_file_flag = ()
        if self.primer_bed_file:
            (bed_file, ) = self.primer_bed_file.file_path
            bed_file_flag = ("-b", os.path.abspath(bed_file))
        amplicon_info_flag = ()
        if self.amplicon_info_file:
            (tsv_file, ) = self.amplicon_info_file.file_path
            amplicon_info_flag = ("-f", os.path.abspath(tsv_file))
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (trimmed_bam_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            *bed_file_flag,
            *amplicon_info_flag,
            "-i", sort_bam_file,
            "-p", trimmed_bam_file[:-4]
        )

    @classmethod
    def set_primer_file(cls, primer_bed_file: _FileArtifacts) -> None:
        cls.primer_bed_file = primer_bed_file

    @classmethod
    def set_amplicon_info(cls, amplicon_info_file: _FileArtifacts) -> None:
        cls.amplicon_info_file = amplicon_info_file

    def _post_execution(self) -> None:
        pass


class IvarGetMaskedWrapper(CMDwrapperBase):

    exec_name = "ivar getmasked"
    primer_bed_file: _FileArtifacts = None
    amplicon_info_file: _FileArtifacts = None
    
    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (vcf_file, ) = input_files.path_from_cwd()
        (bed_file, ) = self.primer_bed_file.file_path
        (tsv_file, ) = self.amplicon_info_file.file_path
        self.output_files = input_files.coupled_files(
            ".txt",
            exec_name=self.exec_name
        )
        (masked_primers, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-b", os.path.abspath(bed_file),
            "-f", os.path.abspath(tsv_file),
            "-i", vcf_file,
            "-p", masked_primers[:-4]
        )

    @classmethod
    def set_primer_file(cls, primer_bed_file: _FileArtifacts) -> None:
        cls.primer_bed_file = primer_bed_file

    @classmethod
    def set_amplicon_info(cls, amplicon_info_file: _FileArtifacts) -> None:
        cls.amplicon_info_file = amplicon_info_file

    def _post_execution(self) -> None:
        pass


class IvarRemoveReadsWrapper(CMDwrapperBase):

    exec_name = "ivar removereads"
    primer_bed_file: _FileArtifacts = None

    def __init__(
        self,
        input_files: _FileArtifacts,
        masked_primers_file: _FileArtifacts,
        *args: str
    ) -> None:
        (sort_bam_file, ) = input_files.file_path
        (masked_primers, ) = masked_primers_file.file_path
        (bed_file, ) = self.primer_bed_file.file_path
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (trimmed_bam_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-t", os.path.abspath(masked_primers),
            "-b", os.path.abspath(bed_file),
            "-i", sort_bam_file,
            "-p", trimmed_bam_file[:-4]
        )

    @classmethod
    def set_primer_file(cls, primer_bed_file: _FileArtifacts) -> None:
        cls.primer_bed_file = primer_bed_file

    def _post_execution(self) -> None:
        pass


class VCFintersectWrapper(CMDwrapperBase):

    exec_name = "vcfintersect"

    def __init__(
        self,
        input_files: _FileArtifacts,
        vcf_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        union: bool,
        *args: str
    ) -> None:
        (vcf_file_1, ) = input_files.path_from_cwd()
        (vcf_file_2, ) = vcf_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        self.output_files = vcf_files.coupled_files(
            ".vcf",
            exec_name="union" if union else "intersect"
        )
        super().__init__(
            input_files,
            *args,
            "-r", os.path.abspath(refSeq_file),
            "-u" if union else "-i", vcf_file_1,
            vcf_file_2
        )

    def _post_execution(self) -> None:
        (vcf_file, ) = self.output_files.file_path
        with open(vcf_file, "w") as f:
            f.writelines(self.stdout)


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
