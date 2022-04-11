import os
import sys
import glob
import shutil
import pathlib
import argparse


def check_base_parser(sysargs, parser: argparse.ArgumentParser):
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    input_dir = args.input
    out_dir = args.outdir
    max_errors = args.maxErrors
    input_ext = args.inputExt

    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"{input_dir} does not exist.")

    csv_dir = os.path.join(out_dir, "csv")
    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    error_dir = os.path.join(out_dir, "errors")
    # Need to clear the error_dir for error counting
    if os.path.exists(error_dir):
        shutil.rmtree(error_dir)
    os.makedirs(error_dir)

    if input_ext:
        input_ext = "." + input_ext

    return {
        "input_dir": input_dir,
        "out_dir": out_dir,
        "csv_dir": csv_dir,
        "error_dir": error_dir,
        "max_errors": max_errors,
        "input_ext": input_ext
    }, args


def check_common_flags(sysargs, parser: argparse.ArgumentParser):
    params, args = check_base_parser(sysargs, parser)

    refSeq_file = args.refSeq
    cores = args.cores
    threads = args.threads
    out_dir = params.pop("out_dir")

    if not os.path.exists(refSeq_file):
        raise FileNotFoundError(f"{refSeq_file} does not exist.")

    output_dir = os.path.join(out_dir, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return {
        "refSeq_file": refSeq_file,
        "output_dir": output_dir,
        "cores": cores,
        "threads": threads,
        **params,
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
    sra_ids = args.sra
    sra_dir = params.pop("input_dir")

    sra_ext = args.SRAext  # TODO deprecate --SRAext with --inputExt for general purpose
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
    forward_ext = args.forward
    reverse_ext = args.reverse
    fastq_dir = params.pop("input_dir")

    if forward_ext == reverse_ext:
        raise ValueError(
            f"Forward ({forward_ext}) and reverse ({reverse_ext}) trailing name can't be the same")

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
