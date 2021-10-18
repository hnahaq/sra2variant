import os
import sys
from multiprocessing import Pool

from sra2variant.pipeline.cmd_wrapper import _FileArtifacts, CMDwrapperBase
from sra2variant.pipeline.sra2fastq import FastqDumpWrapper
from sra2variant.pipeline.fastq_trim import FastpWrapper
from sra2variant.pipeline.reads_mapping import (
    BWAindexWrapper,
    BWAmemWrapper
)
from sra2variant.pipeline.sam_trim import (
    SAMtoolsSortWrapper,
    SAMtoolsIndexWrapper,
    SAMtoolsViewWrapper,
    PicardMarkDuplicatedWrapper
)
from sra2variant.pipeline.variant_calling import (
    LoFreqFaidxWrapper,
    LoFreqViterbiWraper,
    LoFreqIndelQualWrapper,
    LoFreqCallWrapper,
    LoFreqFilterWrapper
)
from sra2variant.vcfparser.vcf2csv import vcf2csv
from sra2variant.sra2variant_cmd import (
    common_flags,
    common_check,
    init_sra_file,
    ErrorTolerance,
)


def run_workflow(
    refSeq_file: str,
    sra_files: tuple,
    output_dir: str,
    csv_dir: str,
    cores: int,
    threads: int,
    error_dir: str,
    max_errors: int = 0,
):
    fasta_file = _FileArtifacts(
        refSeq_file,
        cwd=os.path.dirname(refSeq_file),
        working_id=os.path.basename(refSeq_file)
    )
    if not fasta_file.exist():
        raise FileNotFoundError(f"Reference genome {refSeq_file} not found")
    bwa_index = BWAindexWrapper(fasta_file)
    LoFreqFaidxWrapper(fasta_file).execute_cmd()

    CMDwrapperBase.set_threads(threads)
    BWAmemWrapper.set_base_index(bwa_index)

    with Pool(cores) as p:
        p.starmap(
            func=__workflow,
            iterable=((fn, output_dir, csv_dir, error_dir, max_errors)
                      for fn in sra_files)
        )


def __workflow(
    sra_file: str,
    output_dir: str,
    csv_dir: str,
    error_dir: str,
    max_errors: int,
) -> None:
    sra_file: _FileArtifacts = init_sra_file(
        sra_file,
        output_dir,
        csv_dir
    )
    if sra_file.result_exists():
        print(f"Result for {sra_file} exists in {sra_file.result_file()}")
        return None
    sra_file.create_cwd()
    task_log_file = sra_file.log_file()
    try:
        fastq_files = FastqDumpWrapper(sra_file).execute_cmd()
        if len(fastq_files) != 2:
            print(f"{sra_file} doesn't contain paired end reads")
            return None
        fastq_files = FastpWrapper(
            fastq_files,
            "--detect_adapter_for_pe"
        ).execute_cmd()
        sam_files = BWAmemWrapper(fastq_files).execute_cmd()
        refSeq_file = BWAmemWrapper.bwa_index.output_files

        bam_files = SAMtoolsSortWrapper(sam_files).execute_cmd()
        bam_files = SAMtoolsIndexWrapper(bam_files).execute_cmd()
        bam_files = SAMtoolsViewWrapper(
            bam_files,
            "-q", "20",
            "-f", "3"
        ).execute_cmd()

        bam_files = PicardMarkDuplicatedWrapper(
            bam_files,
            "REMOVE_DUPLICATES=true",
            "ASSUME_SORTED=true",
            "DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES",
            "OPTICAL_DUPLICATE_PIXEL_DISTANCE=100",
            "VALIDATION_STRINGENCY=LENIENT",
        ).execute_cmd()
        bam_files = LoFreqViterbiWraper(
            bam_files,
            refSeq_file,
            "--defqual", "2"
        ).execute_cmd()
        bam_files = SAMtoolsSortWrapper(bam_files, "--no-PG").execute_cmd()
        bam_files = LoFreqIndelQualWrapper(
            bam_files,
            refSeq_file,
            "--dindel"
        ).execute_cmd()
        bam_files = SAMtoolsIndexWrapper(bam_files).execute_cmd()
        vcf_files = LoFreqCallWrapper(
            bam_files,
            refSeq_file,
            "--call-indels",
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
        ).execute_cmd()
        vcf_files = LoFreqFilterWrapper(
            vcf_files,
            "--no-defaults",
            "--print-all",
            "-a", "0.0",
            "-A", "0.0",
            "-b", "fdr",
            "-c", "0.001"
        ).execute_cmd()
        vcf2csv(vcf_files)
    except Exception as e:
        error_tolerance = ErrorTolerance(error_dir, max_errors)
        error_tolerance.handle(e, task_log_file)


def main(sysargs=sys.argv[1:]):
    parser = common_flags("WGS PE pipeline")
    params, _ = common_check(sysargs, parser)

    run_workflow(**params)
