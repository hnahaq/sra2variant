import pathlib
import sys
import logging
from multiprocessing import Pool

from sra2variant.pipeline.cmd_wrapper import CMDwrapperBase, ErrorTolerance
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
from sra2variant.artifacts.base_file import _FileArtifacts
from sra2variant.artifacts.fasta_file import init_fasta_file
from sra2variant.artifacts.sra_file import init_sra_file
from sra2variant.artifacts.fastq_file import init_fastq_pair
from sra2variant.vcfparser.vcf2csv import vcf2csv
from sra2variant.cmdutils.parser_maker import (
    common_flags,
    sra_flags,
    fastq_flags,
)
from sra2variant.cmdutils.parser_checker import (
    check_common_flags,
    check_sra_flags,
    check_fastq_flags,
)


def prep_workflow(
    refSeq_file: str,
    threads: int,
    max_errors: int = 0,
):
    fasta_file = init_fasta_file(refSeq_file)
    bwa_index = BWAindexWrapper(fasta_file)
    LoFreqFaidxWrapper(fasta_file).execute_cmd()

    CMDwrapperBase.set_threads(threads)
    ErrorTolerance.set_max_errors(max_errors)
    BWAmemWrapper.set_base_index(bwa_index)


def workflow_from_sra(
    sra_file: str,
    output_dir: str,
    csv_dir: str,
    error_dir: str,
) -> None:
    sra_file: _FileArtifacts = init_sra_file(
        sra_file=sra_file,
        output_dir=output_dir,
        csv_dir=csv_dir
    )
    task_log_file = sra_file.log_file()
    try:
        if not sra_file.create_cwd():
            return None
        fastq_files = FastqDumpWrapper(sra_file).execute_cmd()
        # line 74-76 repaired by hanna,2023/7/28
        # if len(fastq_files) != 2:
            # logging.warning(f"{sra_file} doesn't contain paired end reads")
            # return None
        __workflow(fastq_files)
    except Exception as e:
        error_tolerance = ErrorTolerance(error_dir, task_log_file)
        error_tolerance.handle(e)


def workflow_from_fastq(
    fastq_forward: pathlib.Path,
    forward_ext: str,
    reverse_ext: str,
    output_dir: str,
    csv_dir: str,
    error_dir: str,
):
    fastq_files: _FileArtifacts = init_fastq_pair(
        fastq_forward=fastq_forward,
        forward_ext=forward_ext,
        reverse_ext=reverse_ext,
        output_dir=output_dir,
        csv_dir=csv_dir
    )
    task_log_file = fastq_files.log_file()
    try:
        if not fastq_files.create_cwd():
            return None
        __workflow(fastq_files)
    except Exception as e:
        error_tolerance = ErrorTolerance(error_dir, task_log_file)
        error_tolerance.handle(e)


def __workflow(fastq_files: _FileArtifacts) -> None:
    # line 110-122 repaired by hanna,2023/7/28
    # fastq_files = FastpWrapper(
        # fastq_files,
        # "--detect_adapter_for_pe"
    # ).execute_cmd()
    if len(fastq_files) == 1:
        fastq_files = FastpWrapper(
            fastq_files,
        ).execute_cmd()
    elif len(fastq_files) == 2:
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


def main(sysargs=sys.argv[1:]):
    parser = common_flags("WGS PE pipeline")
    parser = sra_flags(parser)

    params, args = check_common_flags(sysargs, parser)
    params = check_sra_flags(args, params)

    prep_workflow(
        refSeq_file=params["refSeq_file"],
        threads=params["threads"],
        max_errors=params["max_errors"]
    )

    output_dir = params["output_dir"]
    csv_dir = params["csv_dir"]
    error_dir = params["error_dir"]

    with Pool(params["cores"]) as p:
        p.starmap(func=workflow_from_sra,
                  iterable=((fn, output_dir, csv_dir, error_dir)
                            for fn in params["sra_files"]))


def main_fastq(sysargs=sys.argv[1:]):
    parser = common_flags("WGS PE pipeline from fastq pairs")
    parser = fastq_flags(parser)

    params, args = check_common_flags(sysargs, parser)
    params = check_fastq_flags(args, params)

    prep_workflow(
        refSeq_file=params["refSeq_file"],
        threads=params["threads"],
        max_errors=params["max_errors"]
    )

    func_extra_args = (
        params["forward_ext"],
        params["reverse_ext"],
        params["output_dir"],
        params["csv_dir"],
        params["error_dir"]
    )

    with Pool(params["cores"]) as p:
        p.starmap(func=workflow_from_fastq,
                  iterable=((p, *func_extra_args)
                            for p in params["fastq_forward"]))
