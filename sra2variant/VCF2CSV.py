# import sys
# import pathlib

# from sra2variant.pipeline.cmd_wrapper import _FileArtifacts
# from sra2variant.vcfparser.vcf2csv import vcf2csv
# from sra2variant.sra2variant_cmd import (
#     base_parser,
#     vcf_flags,
#     check_base_parser,
#     check_vcf_flags,
# )


# def main(sysargs=sys.argv[1:]):
#     parser = base_parser("Convert VCF to CSV")
#     parser = vcf_flags(parser)

#     params, args = check_base_parser(sysargs, parser)
#     params = check_vcf_flags(args, params)

#     csv_dir = params["csv_dir"],
#     error_dir = params["error_dir"]
#     max_errors = params["max_errors"]

#     for vcf_file in pathlib.Path(vcf_dir).glob(f"**/*.vcf"):
#         print(vcf_file)
