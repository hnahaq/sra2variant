import sys
import pathlib

from sra2variant.pipeline.cmd_wrapper import ErrorTolerance
from sra2variant.artifacts.base_file import _FileArtifacts
from sra2variant.artifacts.vcf_file import init_vcf_file
from sra2variant.vcfparser.vcf2csv import vcf2csv
from sra2variant.cmdutils.parser_maker import base_parser
from sra2variant.cmdutils.parser_checker import check_base_parser


def main(sysargs=sys.argv[1:]):
    parser = base_parser("Convert VCF to CSV")

    params, _ = check_base_parser(sysargs, parser)

    vcf_dir = params["input_dir"]
    csv_dir = params["csv_dir"]
    error_dir = params["error_dir"]
    max_errors = params["max_errors"]

    ErrorTolerance.set_max_errors(max_errors)

    for vcf_file in pathlib.Path(vcf_dir).glob(f"**/*.vcf"):
        vcf_file: _FileArtifacts = init_vcf_file(
            vcf_file=vcf_file,
            csv_dir=csv_dir
        )
        task_log_file = vcf_file.log_file()
        try:
            vcf2csv(vcf_file)
        except Exception as e:
            error_tolerance = ErrorTolerance(error_dir, task_log_file)
            error_tolerance.handle(e)
