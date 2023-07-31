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
)
from sra2variant.pipeline.variant_calling import (
    LoFreqFaidxWrapper,
    LoFreqViterbiWraper,
    LoFreqIndelQualWrapper,
    LoFreqCallWrapper,
    LoFreqFilterWrapper,
    IvarTrimWrapper,
    IvarGetMaskedWrapper,
    IvarRemoveReadsWrapper,
    VCFintersectWrapper,
    sanitize_bed,
    prepare_amplicon_info,
    write_amplicon_info_file,
    completemask,
)
from sra2variant.artifacts.base_file import _FileArtifacts
from sra2variant.artifacts.fasta_file import init_fasta_file
from sra2variant.artifacts.bed_file import init_bed_file
from sra2variant.artifacts.tsv_file import init_tsv_file
from sra2variant.artifacts.sra_file import init_sra_file
from sra2variant.artifacts.fastq_file import init_fastq_pair
from sra2variant.vcfparser.vcf2csv import vcf2csv
from sra2variant.cmdutils.parser_checker import (
    check_common_flags,
    check_artic_flags,
    check_sra_flags,
    check_fastq_flags,
)
from sra2variant.cmdutils.parser_maker import (
    common_flags,
    artic_flags,
    sra_flags,
    fastq_flags,
)


def prep_workflow(**params) -> None:

    refSeq_file: str = params["refSeq_file"]
    primer_file: str = params["primer_file"]
    amplicon_info_file: str = params["amplicon_info_file"]
    threads: int = params["threads"]
    max_errors: int = params["max_errors"]

    fasta_file = init_fasta_file(refSeq_file)
    bed_file = init_bed_file(primer_file)
    tsv_file = init_tsv_file(amplicon_info_file)

    bed_file = sanitize_bed(bed_file)
    tsv_file_2 = prepare_amplicon_info(bed_file, tsv_file)
    tsv_file_3 = write_amplicon_info_file(bed_file)

    IvarTrimWrapper.set_primer_file(bed_file)
    IvarTrimWrapper.set_amplicon_info(tsv_file_2)

    IvarGetMaskedWrapper.set_primer_file(bed_file)
    IvarGetMaskedWrapper.set_amplicon_info(tsv_file_3)

    IvarRemoveReadsWrapper.set_primer_file(bed_file)

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
        # line 103-105 repaired by hanna,2023/7/28
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
    # line 139-151 repaired by hanna,2023/7/28
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
        "-f", "1",
        "-F", "268"
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
    bam_files = IvarTrimWrapper(
        bam_files,
        "-x", "0",
        "-e",
        "-m", "1",
        "-q", "0",
        "-s", "4",
    ).execute_cmd()
    bam_files = SAMtoolsSortWrapper(bam_files).execute_cmd()
    bam_files = SAMtoolsIndexWrapper(bam_files).execute_cmd()
    vcf_files_1 = LoFreqCallWrapper(
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
    ).execute_cmd()
    vcf_files = LoFreqFilterWrapper(
        vcf_files_1,
        "--no-defaults",
        "--verbose",
        "-v", "5",
        "-V", "0",
        "-a", "0.1",
        "-A", "0.9",
    ).execute_cmd()
    ivar_getmasked = IvarGetMaskedWrapper(vcf_files)
    ivar_getmasked.execute_cmd()
    masked_primers_file = completemask(ivar_getmasked)
    bam_files = IvarRemoveReadsWrapper(
        bam_files,
        masked_primers_file,
    ).execute_cmd()
    bam_files = SAMtoolsIndexWrapper(bam_files).execute_cmd()
    vcf_files_2 = LoFreqCallWrapper(
        bam_files,
        refSeq_file,
        "--call-indels",
        "--min-cov", "20",
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
    ).execute_cmd()
    vcf_files_3 = VCFintersectWrapper(
        vcf_files_2,
        vcf_files_1,
        refSeq_file,
        False,
        "-v",
        "-w", "0",
    ).execute_cmd()
    vcf_files = VCFintersectWrapper(
        vcf_files_3,
        vcf_files_2,
        refSeq_file,
        True,
        "-w", "0"
    ).execute_cmd()
    vcf_files = LoFreqFilterWrapper(
        vcf_files,
        "--no-defaults",
        "--verbose",
        "--print-all",
        "-v", "0",
        "-V", "0",
        "-a", "0.0",
        "-A", "0.0",
        "-b", "fdr",
        "-c", "0.001"
    ).execute_cmd()
    vcf2csv(vcf_files)


def main(sysargs=sys.argv[1:]):
    parser = common_flags("ARTIC PE pipeline")
    parser = artic_flags(parser)
    parser = sra_flags(parser)

    params, args = check_common_flags(sysargs, parser)
    params = check_artic_flags(args, params)
    params = check_sra_flags(args, params)

    prep_workflow(**params)

    func_extra_args = (
        params["output_dir"],
        params["csv_dir"],
        params["error_dir"]
    )

    with Pool(params["cores"]) as p:
        p.starmap(func=workflow_from_sra,
                  iterable=((fn, *func_extra_args)
                            for fn in params["sra_files"]))


def main_fastq(sysargs=sys.argv[1:]):
    parser = common_flags("ARTIC PE pipeline from fastq pairs")
    parser = artic_flags(parser)
    parser = fastq_flags(parser)

    params, args = check_common_flags(sysargs, parser)
    params = check_artic_flags(args, params)
    params = check_fastq_flags(args, params)

    prep_workflow(**params)

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
