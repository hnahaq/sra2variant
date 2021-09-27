import os
from copy import deepcopy
from collections import defaultdict

import numpy as np
import pysam
import vcf

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


"""
gatk CreateSequenceDictionary -R NC_045512.2.fasta
gatk AddOrReplaceReadGroups I=SRR15432222.sort.bam O=SRR15432222_with_RG.sort.bam RGID=4 RGLB=lib1 RGPL=illumina RGSM=1 RGPU=unit1 CREATE_INDEX=True
gatk HaplotypeCaller -R NC_045512.2.fasta -I SRR15432222.sort.bam -O SRR15432222.vcf -ERC GVCF -ploidy 1
"""


class LoFreqViterbiWraper(CMDwrapperBase):

    exec_name: str = "lofreq viterbi"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts
    ) -> None:
        (sort_bam_file, ) = input_files.path_from_cwd()
        self.output_files = input_files.coupled_files(("_lofreq_viterbi", ))
        (viterbi_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--defqual", "2",
            *self.exec_args,
            "--ref", os.path.abspath(refSeq_file),
            "--out", viterbi_file,
            sort_bam_file,
        )

    def post_execution(self) -> None:
        (viterbi_file, ) = self.output_files.file_path
        self.output_files = self.output_files.coupled_files(
            ("_lofreq_viterbi.sort.bam", )
        )
        (sort_bam_file, ) = self.output_files.file_path
        pysam.sort("-@", self.threads, "-o", sort_bam_file, viterbi_file)
        pysam.index(sort_bam_file)


class LoFreqIndelQualWrapper(CMDwrapperBase):

    exec_name: str = "lofreq indelqual"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts
    ) -> None:
        self.output_files = input_files.coupled_files(("_lofreq_indelqual", ))
        (sort_bam_file, ) = input_files.path_from_cwd()
        (indelqual_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            *self.exec_args,
            "--dindel",
            "--verbose",
            "--ref", os.path.abspath(refSeq_file),
            "-o", indelqual_file,
            sort_bam_file,
        )

    def post_execution(self) -> None:
        (indelqual_file, ) = self.output_files.file_path
        pysam.index(indelqual_file)


class LoFreqCallWrapper(CMDwrapperBase):

    exec_name: str = "lofreq call"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts
    ) -> None:
        self.output_files = input_files.coupled_files(("_lofreq.vcf", ))
        (indelqual_file, ) = input_files.path_from_cwd()
        (vcf_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--min-cov", "5",
            "--max-depth", "1000000",
            "--min-bq", "30",
            "--min-alt-bq", "30",
            "--min-mq", "20",
            "--max-mq", "255",
            "--min-jq", "0",
            "--min-alt-jq", "0",
            "--def-alt-jq", "0",
            "--sig", "0.0005",
            "--bonf", "dynamic",
            "--no-default-filter",
            "--no-default-filter",
            *self.exec_args,
            "--call-indels",
            "--verbose",
            "--ref", os.path.abspath(refSeq_file),
            "-o", vcf_file,
            indelqual_file
        )

    def post_execution(self) -> None:
        (vcf_file, ) = self.output_files.file_path
        res = []
        with open(vcf_file) as f:
            prev_pos = -1
            curr_pos = 0
            temp_records = defaultdict(list)
            for record in vcf.Reader(f):
                curr_pos = record.POS
                if "DP4" in record.INFO:
                    total_count = sum(record.INFO["DP4"])
                else:
                    total_count = sum(record.INFO["I16"][:4])
                    
                for allele in record.ALT:
                    reference = record.REF
                    allele_seq = allele.sequence
                    if len(reference) > len(allele_seq):
                        variant_type = "Deletion"
                        reference = reference[len(allele_seq):]
                        variant_len = len(reference)
                        allele_seq = "-"
                        record.POS += 1
                    elif len(record.REF) < len(allele.sequence):
                        variant_type = "Insertion"
                        allele_seq = allele_seq[len(reference):]
                        variant_len = len(allele_seq)
                        reference = "-"
                        record.POS += 1
                    else:
                        variant_len = len(allele_seq)
                        if variant_len == 1:
                            variant_type = "SNV"
                        else:
                            variant_type = "MNV"
                    allele_frequency = record.INFO["AF"]
                    curr_record = {
                        "Reference Position": record.POS,
                        "Type": variant_type,
                        "Length": variant_len,
                        "Reference": reference,
                        "Allele": allele_seq,
                        "Linkage": "",
                        "Count": int(allele_frequency * total_count),
                        "Coverage": int(total_count),
                        "Frequency": allele_frequency * 100,
                        "Forward/reverse balance": record.INFO["SB"],
                        "Average quality": record.QUAL,
                        "Overlapping annotations": "N/A",
                        "Coding region change": "N/A",
                        "Amino acid change": "N/A"
                    }
                    if (curr_pos - prev_pos <= 1 and
                        curr_record["Type"] == "SNV" and
                        prev_record["Type"] == "SNV"):
                        temp_records[curr_record["Reference Position"]].append(curr_record)
                    else:
                        if prev_pos > 0:
                            res.extend(resolve_MNV(temp_records))
                        temp_records = defaultdict(list)
                        temp_records[curr_record["Reference Position"]].append(curr_record)
                prev_record = curr_record
                prev_pos = curr_pos
            if len(temp_records):
                res.extend(resolve_MNV(temp_records))
        self.output_files._write_result(res)


class BCFtoolsMpileupWrapper(CMDwrapperBase):

    exec_name: str = "bcftools mpileup"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts
    ) -> None:
        self.output_files = input_files.coupled_files(("_bcftools_mpileup", ))
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
            *self.exec_args,
            "--annotate", "INFO/AD,INFO/ADF,INFO/ADR",
            "--output", mpileup_file,
            sort_bam_file
        )

    def post_execution(self) -> None:
        pass


# class BCFtoolsCallWrapper(CMDwrapperBase):

#     exec_name: str = "bcftools call"

#     def __init__(
#         self,
#         input_files: _FileArtifacts,
#     ) -> None:
#         self.output_files = input_files.coupled_files(("_bcftools_call", ))
#         (mpileup_file, ) = input_files.path_from_cwd()
#         (call_bcf_file, ) = self.output_files.path_from_cwd()
#         super().__init__(
#             input_files,
#             "--threads", self.threads,
#             "--output-type", "b",
#             "--prior", "0",
#             "--ploidy", "1",
#             *self.exec_args,
#             "--keep-alts",
#             "--variants-only",
#             "--multiallelic-caller",
#             "--output", call_bcf_file,
#             mpileup_file
#         )

#     def post_execution(self) -> None:
#         pass


class BCFtoolsNormWrapper(CMDwrapperBase):

    exec_name: str = "bcftools norm"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts
    ) -> None:
        self.output_files = input_files.coupled_files(("_bcftools_norm", ))
        (call_bcf_file, ) = input_files.path_from_cwd()
        (norm_vcf_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--threads", self.threads,
            "--multiallelics", "-any",
            "--output-type", "v",
            "--fasta-ref", os.path.abspath(refSeq_file),
            *self.exec_args,
            "--output", norm_vcf_file,
            call_bcf_file
        )

    def post_execution(self) -> None:
        (vcf_file, ) = self.output_files.file_path
        res = []
        with open(vcf_file) as f:
            for record in vcf.Reader(f):
                if "DP4" in record.INFO:
                    total_count = sum(record.INFO["DP4"])
                else:
                    total_count = sum(record.INFO["I16"][:4])
                for n, allele in enumerate(record.ALT):
                    if allele.type == "*":
                        continue
                    reference = record.REF
                    allele_seq = allele.sequence
                    if len(reference) > len(allele_seq):
                        variant_type = "Deletion"
                        reference = reference[len(allele_seq):]
                        variant_len = len(reference)
                        allele_seq = "-"
                    elif len(record.REF) < len(allele.sequence):
                        variant_type = "Insertion"
                        allele_seq = allele_seq[len(reference):]
                        variant_len = len(allele_seq)
                        reference = "-"
                    else:
                        variant_len = len(allele_seq)
                        if variant_len == 1:
                            variant_type = "SNV"
                        else:
                            variant_type = "MNV"
                    allele_count = record.INFO["AD"][n + 1]
                    res.append({
                        "Reference Position": record.POS,
                        "Type": variant_type,
                        "Length": variant_len,
                        "Reference": reference,
                        "Allele": allele_seq,
                        "Linkage": "",
                        "Count": allele_count,
                        "Coverage": int(total_count),
                        "Frequency": allele_count / total_count,
                        "Forward/reverse balance": "N/A",
                        "Average quality": record.INFO["QS"][n + 1],
                        "Overlapping annotations": "N/A",
                        "Coding region change": "N/A",
                        "Amino acid change": "N/A"
                    })
        self.output_files._write_result(res)


def resolve_MNV(temp_records: dict) -> list:
    
    if len(temp_records) == 1:
        (records, ) = temp_records.values()
        return records
    
    temp_records = list(temp_records.values())
    sort_i = list(range(len(temp_records)))
    for records in temp_records:
        for _ in records:
            swap_flag0 = False
            i = 0
            while i < len(records) - 1:
                if records[i]["Frequency"] < records[i + 1]["Frequency"]:
                    records[i], records[i + 1] = records[i + 1], records[i]
                    swap_flag0 = True
                i += 1
            if not swap_flag0:
                break
        j = 0
        while j < len(temp_records) - 1:
            if len(temp_records[j]) < len(temp_records[j + 1]):
                temp_records[j], temp_records[j + 1] = temp_records[j + 1], temp_records[j]
                sort_i[j], sort_i[j + 1] = sort_i[j + 1], sort_i[j]
            j += 1
    
    for _ in temp_records:
        swap_flag1 = False
        k = 0
        while k < len(temp_records) - 1:
            if (temp_records[k][0]["Frequency"] < temp_records[k + 1][0]["Frequency"] and
                len(temp_records[k]) == len(temp_records[k + 1])):
                temp_records[k], temp_records[k + 1] = temp_records[k + 1], temp_records[k]
                sort_i[k], sort_i[k + 1] = sort_i[k + 1], sort_i[k]
                swap_flag1 = True
            k += 1
        if not swap_flag1:
            break
            
    ref_records = temp_records[0]
    res = deepcopy(temp_records[0])
    for records in temp_records[1:]:
        init_compare = 0
        for rec in records:
            min_diff = float('inf')
            compare_start, compare_end = init_compare, init_compare + 1
            ref_freq = 0
            for n, ref_rec in enumerate(ref_records, start=init_compare):
                ref_freq = ref_rec["Frequency"]
                diff = abs(ref_freq - rec["Frequency"])
                if diff < min_diff:
                    min_diff = diff
                    compare_start, compare_end = n, n + 1
            ref_acc_freq = 0
            for n, ref_rec in enumerate(ref_records, start=init_compare):
                ref_acc_freq += ref_rec["Frequency"]
                diff = abs(ref_acc_freq - rec["Frequency"])
                if diff < min_diff:
                    min_diff = diff
                    compare_start, compare_end = init_compare, n + 1
                
            init_compare = compare_end
                    
            for ref_rec_index in range(compare_start, compare_end):
                if type(res[ref_rec_index]["Reference Position"]) is int:
                    res[ref_rec_index]["Reference Position"] = [
                        res[ref_rec_index]["Reference Position"],
                        rec["Reference Position"]
                    ]
                else:
                    res[ref_rec_index]["Reference Position"].append(rec["Reference Position"])

                res[ref_rec_index]["Type"] = "MNV"
                res[ref_rec_index]["Length"] += 1
                res[ref_rec_index]["Reference"] += rec["Reference"]
                res[ref_rec_index]["Allele"] += rec["Allele"]

#                 res[ref_rec_index]["Coverage"] += rec["Coverage"]
#                 res[ref_rec_index]["Frequency"] += rec["Frequency"]
                res[ref_rec_index]["Average quality"] += rec["Average quality"]
    
    for ref_rec in res:
        if (type(ref_rec["Reference Position"]) is int):
            continue
        else:
            sort_i2 = np.argsort(ref_rec["Reference Position"])
            all(i == j for i, j in zip(sort_i, sort_i2))
            ref_rec["Reference Position"] = min(ref_rec["Reference Position"])
            ref_rec["Reference"] = "".join(ref_rec["Reference"][i] for i in sort_i)
            ref_rec["Allele"] = "".join(ref_rec["Allele"][i] for i in sort_i)
    #         ref_rec["Coverage"] /= len(sort_i)
    #         ref_rec["Coverage"] = int(ref_rec["Coverage"])
            ref_rec["Average quality"] /= len(sort_i)
    return res
