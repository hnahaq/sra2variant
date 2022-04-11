import os
import argparse

from sra2variant import __version__


def base_parser(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version",
                        version=f"sra2variant {__version__}")
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Output directory")
    parser.add_argument("--maxErrors", type=int, default=10,
                        help="Max number of errors before stop the program")
    parser.add_argument("--inputExt", default=None,
                        help="File extension for searching target files in INPUT directory")
    return parser


def common_flags(description):
    parser = base_parser(description)
    parser.add_argument("-r", "--refSeq", required=True,
                        help="File path of reference genome")
    parser.add_argument("-c", "--cores", type=int, default=os.cpu_count(),
                        help="Number of paralle tasks")
    parser.add_argument("-t", "--threads", type=int, default=2,
                        help="Number of threads for each task")
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
    parser.add_argument("--SRAext", default="sra",  # TODO deprecate --SRAext with --inputExt for general purpose
                        help="File extension for searching sra files in INPUT directory")
    return parser


def fastq_flags(parser: argparse.ArgumentParser):
    parser.add_argument("--forward", default="_1.fastq",
                        help="The trailing name for forward fastq files in INPUT directory")
    parser.add_argument("--reverse", default="_2.fastq",
                        help="The trailing name for reverse fastq files in INPUT directory")
    return parser
