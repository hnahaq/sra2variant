import os
import sys
import glob
import argparse
from functools import partial
from multiprocessing import Pool, cpu_count

from sra2variant.pipeline.cmd_wrapper import _FileArtifacts, CMDwrapperBase
from sra2variant.pipeline.sra2fastq import FastqDumpWrapper
from sra2variant.pipeline.fastq_trim import FastpWrapper
from sra2variant.pipeline.reads_mapping import (
    BWAindexWrapper,
    BWAmemWrapper
)
from sra2variant.pipeline.sam_trim import PicardMarkDuplicatedWrapper
from sra2variant.pipeline.variant_calling import (
    LoFreqViterbiWraper,
    LoFreqIndelQualWrapper,
    LoFreqCallWrapper
)
from sra2variant import __version__


def workflow(
    sra_file: str,
    output_dir: str,
    csv_dir: str,
) -> None:
    sra_id = os.path.basename(sra_file)
    if sra_file.endswith(".sra"):
        (sra_id, _) = os.path.splitext(sra_id)
    sra_file: _FileArtifacts = _FileArtifacts(
        file_path=(sra_file, ),
        cwd=os.path.join(output_dir, sra_id),
        working_id=sra_id,
        res_dir=csv_dir
    )
    if sra_file.result_exists():
        print(f"Result of {sra_file} already exists in {sra_file.result_file()}")
        return None
    sra_file.create_cwd()
    fastq_files = FastqDumpWrapper(sra_file).execute_cmd()
    if len(fastq_files) != 2:
        return None
    fastq_files = FastpWrapper(fastq_files).execute_cmd()
    bam_files = BWAmemWrapper(fastq_files).execute_cmd()

    refSeq_file = BWAmemWrapper.refSeq_file

    bam_files = PicardMarkDuplicatedWrapper(bam_files).execute_cmd()
    bam_files = LoFreqViterbiWraper(bam_files, refSeq_file).execute_cmd()
    bam_files = LoFreqIndelQualWrapper(bam_files, refSeq_file).execute_cmd()
    vcf_files = LoFreqCallWrapper(bam_files, refSeq_file).execute_cmd()
    del vcf_files


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="SRA to CSV pipeline")
    parser.add_argument("-v", "--version", action="version",
                        version=f"sra2variant {__version__}")
    parser.add_argument("-r", "--refSeq", required=True,
                        help="File path of reference genome")
    parser.add_argument("-i", "--input", required=True,
                        help="Input SRA directory")
    parser.add_argument("-o", "--output", default="csv",
                        help="Output CSV directory")
    parser.add_argument("-e", "--temp", default="output",
                        help="Temporary output directory")
    parser.add_argument("-s", "--sra", default=None,
                        help="List of SRA ids or path to a file")
    parser.add_argument("-c", "--cores", type=int, default=cpu_count(),
                        help="Number of cores for spawning process")
    parser.add_argument("-t", "--threads", type=int, default=2,
                        help="Number of threads for each process")
    parser.add_argument("--SRAext", default="sra",
                        help="File extension for searching sra files in INPUT directory")
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    args = parser.parse_args()

    refSeq_file = args.refSeq
    sra_dir = args.input
    csv_dir = args.output
    output_dir = args.temp
    sra_ids = args.sra
    cores = args.cores
    threads = args.threads
    sra_ext = args.SRAext

    if sra_dir is None:
        raise Exception("SRA directory not provided.")
    if not os.path.exists(sra_dir):
        raise Exception(f"{sra_dir} does not exist.")

    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if sra_ids is not None:
        if os.path.exists(sra_ids):
            # Assume the given argument is a file path
            with open(sra_ids) as f:
                sra_files = [os.path.join(sra_dir, f"{row.strip()}.{sra_ext}")
                             for row in f]
        else:
            # Assume the given argument is a list of SRA ids
            sra_files = [os.path.join(sra_dir, f"{i}.{sra_ext}")
                         for i in sra_ids.split(",")]
    else:
        # Use all SRA file when no ids provided
        sra_files = glob.glob(os.path.join(sra_dir, f"*.{sra_ext}"))
    print("Number of sra file to process", len(sra_files))

    CMDwrapperBase.set_threads(threads)
    bwaindex = BWAindexWrapper(refSeq_file)
    BWAmemWrapper.set_base_index(bwaindex)
    workflow2 = partial(
        workflow,
        output_dir=output_dir,
        csv_dir=csv_dir,
    )
    with Pool(cores) as p:
        p.map(workflow2, sra_files)
