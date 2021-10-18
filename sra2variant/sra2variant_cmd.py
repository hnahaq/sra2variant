import os
import shutil
import sys
import glob
import argparse

from sra2variant.pipeline.cmd_wrapper import _FileArtifacts
from sra2variant import __version__


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
    parser.add_argument("-s", "--sra", default=None,
                        help="List of SRA ids or path to a file")
    parser.add_argument("-c", "--cores", type=int, default=os.cpu_count(),
                        help="Number of paralle tasks")
    parser.add_argument("-t", "--threads", type=int, default=2,
                        help="Number of threads for each task")
    parser.add_argument("--SRAext", default="sra",
                        help="File extension for searching sra files in INPUT directory")
    parser.add_argument("--maxErrors", type=int, default=10,
                        help="Max number of errors before stop the program")
    return parser


def common_check(sysargs, parser: argparse.ArgumentParser):
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    refSeq_file = args.refSeq
    sra_dir = args.input
    out_dir = args.outdir
    sra_ids = args.sra
    cores = args.cores
    threads = args.threads
    sra_ext = args.SRAext
    max_errors = args.maxErrors

    if not os.path.exists(refSeq_file):
        raise FileNotFoundError(f"{refSeq_file} does not exist.")

    if not os.path.exists(sra_dir):
        raise FileNotFoundError(f"{sra_dir} does not exist.")

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
    print("Number of sra file to process:", len(sra_files))
    return {
        "refSeq_file": refSeq_file,
        "sra_files": sra_files,
        "output_dir": output_dir,
        "csv_dir": csv_dir,
        "cores": cores,
        "threads": threads,
        "error_dir": error_dir,
        "max_errors": max_errors,
    }, args


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


class ErrorTolerance:

    __slots__ = ["error_dir", "max_errors"]

    def __init__(self, error_dir: str, max_errors: int) -> None:
        self.error_dir: str = error_dir
        self.max_errors: int = max_errors

    def handle(self, e: Exception, task_log_file: str) -> None:
        with open(task_log_file, "a") as f:
            f.write(f"Raised error: {e}")
        copied_log_file = os.path.join(
            self.error_dir,
            os.path.basename(task_log_file)
        )
        shutil.copyfile(task_log_file, copied_log_file)
        if len(os.listdir(self.error_dir)) > self.max_errors:
            raise RuntimeError(f"Over {self.max_errors} errors occurred")
