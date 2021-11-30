import os
import shutil
import sys
import glob
import pathlib
import argparse
import logging

from sra2variant.pipeline.cmd_wrapper import _FileArtifacts
from sra2variant import __version__

logging.basicConfig(
    format="[%(asctime)s]: %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO
)


def common_flags(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version",
                        version=f"sra2variant {__version__}")
    parser.add_argument("-r", "--refSeq", required=True,
                        help="File path of reference genome")
    parser.add_argument("-i", "--input", required=True,
                        help="Input SRA directory")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Output directory")
    parser.add_argument("-c", "--cores", type=int, default=os.cpu_count(),
                        help="Number of paralle tasks")
    parser.add_argument("-t", "--threads", type=int, default=2,
                        help="Number of threads for each task")
    parser.add_argument("--maxErrors", type=int, default=10,
                        help="Max number of errors before stop the program")
    return parser


def artic_flags(parser: argparse.ArgumentParser):
    parser.add_argument("-p", "--primers", required=True,
                        help="ARTIC primer bed files")
    parser.add_argument("-a", "--amplicon", required=True,
                        help="ARTIC primers to amplicon assignments")
    return parser


def sra_flags(parser: argparse.ArgumentParser):
    parser.add_argument("-s", "--sra", default=None,
                        help="List of SRA ids or path to a file")
    parser.add_argument("--SRAext", default="sra",
                        help="File extension for searching sra files in INPUT directory")
    return parser


def fastq_flags(parser: argparse.ArgumentParser):
    parser.add_argument("--forward", default="_1.fastq",
                        help="The trailing name for forward fastq files in INPUT directory")
    parser.add_argument("--reverse", default="_2.fastq",
                        help="The trailing name for reverse fastq files in INPUT directory")
    return parser


def check_common_flags(sysargs, parser: argparse.ArgumentParser):
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    refSeq_file = args.refSeq
    out_dir = args.outdir
    cores = args.cores
    threads = args.threads
    max_errors = args.maxErrors

    if not os.path.exists(refSeq_file):
        raise FileNotFoundError(f"{refSeq_file} does not exist.")

    csv_dir = os.path.join(out_dir, "csv")
    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    output_dir = os.path.join(out_dir, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    error_dir = os.path.join(out_dir, "errors")
    # Need to clear the error_dir for error counting
    if os.path.exists(error_dir):
        shutil.rmtree(error_dir)
    os.makedirs(error_dir)
    return {
        "refSeq_file": refSeq_file,
        "output_dir": output_dir,
        "csv_dir": csv_dir,
        "cores": cores,
        "threads": threads,
        "error_dir": error_dir,
        "max_errors": max_errors,
    }, args


def check_artic_flags(args, params: dict):
    primer_file = args.primers
    amplicon_info_file = args.amplicon

    if not os.path.exists(primer_file):
        raise FileNotFoundError(f"{primer_file} does not exist.")

    if not os.path.exists(amplicon_info_file):
        raise FileNotFoundError(f"{amplicon_info_file} does not exist.")

    return {
        **params,
        "primer_file": primer_file,
        "amplicon_info_file": amplicon_info_file
    }


def check_sra_flags(args, params: dict):
    sra_dir = args.input
    sra_ids = args.sra
    sra_ext = args.SRAext

    if not os.path.exists(sra_dir):
        raise FileNotFoundError(f"{sra_dir} does not exist.")

    if sra_ext:
        sra_ext = "." + sra_ext

    if sra_ids is not None:
        if os.path.exists(sra_ids):
            # Assume the given argument is a file path
            with open(sra_ids) as f:
                sra_files = [os.path.join(sra_dir, f"{row.strip()}{sra_ext}")
                             for row in f]
        else:
            # Assume the given argument is a list of SRA ids
            sra_files = [os.path.join(sra_dir, f"{i}{sra_ext}")
                         for i in sra_ids.split(",")]
    else:
        # Use all SRA files when no id list provided
        sra_files = glob.glob(os.path.join(sra_dir, f"*{sra_ext}"))
    # logging.info(f"Number of sra files to process: {len(sra_files)}")
    return {**params, "sra_files": sra_files}


def check_fastq_flags(args, params: dict):
    fastq_dir = args.input
    forward_ext = args.forward
    reverse_ext = args.reverse

    if not os.path.exists(fastq_dir):
        raise FileNotFoundError(f"{fastq_dir} does not exist.")

    if forward_ext == reverse_ext:
        raise ValueError(
            f"Forward ({forward_ext}) and reverse ({reverse_ext}) trailing name can't be the same for pair end fastq")

    # forward_ext_length = len(forward_ext)
    # fastq_pairs = []
    # for forward_f in pathlib.Path(fastq_dir).glob(f"**/*{forward_ext}"):
    #     forward_f = str(forward_f)
    #     fastq_id = forward_f[:-forward_ext_length]
    #     reverse_f = fastq_id + reverse_ext
    #     fastq_id = os.path.basename(fastq_id)
    #     fastq_pairs.append((forward_f, reverse_f, fastq_id))
    # logging.info(f"Number of fastq pairs to process: {len(fastq_pairs)}")
    # return {**params, "fastq_pairs": fastq_pairs}
    return {
        **params,
        "fastq_forward": pathlib.Path(fastq_dir).glob(f"**/*{forward_ext}"),
        "forward_ext": forward_ext,
        "reverse_ext": reverse_ext
    }


def init_sra_file(
    sra_file: str,
    output_dir: str,
    csv_dir: str,
) -> _FileArtifacts:
    sra_id = os.path.basename(sra_file)
    if sra_file.endswith(".sra"):
        (sra_id, _) = os.path.splitext(sra_id)
    return _FileArtifacts(
        sra_file,
        working_id=sra_id,
        cwd=os.path.join(output_dir, sra_id),
        res_dir=csv_dir
    )


def init_fastq_pair(
    fastq_forward: pathlib.Path,
    forward_ext: str,
    reverse_ext: str,
    output_dir: str,
    csv_dir: str
) -> _FileArtifacts:
    fastq_forward = str(fastq_forward)
    fastq_id = fastq_forward[:-len(forward_ext)]
    fastq_reverse = fastq_id + reverse_ext
    fastq_id = os.path.basename(fastq_id)
    return _FileArtifacts(
        fastq_forward,
        fastq_reverse,
        working_id=fastq_id,
        cwd=os.path.join(output_dir, fastq_id),
        res_dir=csv_dir
    )


class ErrorTolerance:

    max_errors: int = 0

    __slots__ = ["error_dir", "task_log_file"]

    def __init__(self, error_dir: str, task_log_file: str) -> None:
        self.error_dir: str = error_dir
        self.task_log_file: str = task_log_file

    @classmethod
    def set_max_errors(cls, max_errors: int) -> None:
        cls.max_errors = max_errors

    def handle(self, e: Exception) -> None:
        with open(self.task_log_file, "a") as f:
            f.write(f"Raised error: {e}")
        copied_log_file = os.path.join(
            self.error_dir,
            os.path.basename(self.task_log_file)
        )
        shutil.copyfile(self.task_log_file, copied_log_file)
        n_errors = len(os.listdir(self.error_dir))
        if n_errors > self.max_errors:
            raise RuntimeError(f"{n_errors} errors occurred exceeding maximum")
