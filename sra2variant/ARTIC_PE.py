import os
import sys
import glob
import re
import argparse
from multiprocessing import Pool, cpu_count

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
)
from sra2variant.vcfparser.vcf2csv import vcf2csv
from sra2variant import __version__


def run_workflow(
    refSeq_file: str,
    primer_file: str,
    amplicon_info_file: str,
    sra_files: tuple,
    output_dir: str,
    csv_dir: str,
    cores: int,
    threads: int
):
    fasta_file = _FileArtifacts(
        refSeq_file,
        cwd=os.path.dirname(refSeq_file),
        working_id=os.path.basename(refSeq_file)
    )
    if not fasta_file.exist():
        raise FileNotFoundError(f"Reference genome {refSeq_file} not found")
    bed_file = _FileArtifacts(
        primer_file,
        cwd=os.path.dirname(primer_file),
    )
    if not bed_file.exist():
        raise FileNotFoundError(f"Primer file {primer_file} not found")

    tsv_file = _FileArtifacts(
        amplicon_info_file,
        cwd=os.path.dirname(amplicon_info_file)
    )
    if not tsv_file.exist():
        raise FileNotFoundError(f"Primer file {amplicon_info_file} not found")
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
    BWAmemWrapper.set_base_index(bwa_index)

    with Pool(cores) as p:
        p.starmap(__workflow, ((fn, output_dir, csv_dir) for fn in sra_files))


def __workflow(
    sra_file: str,
    output_dir: str,
    csv_dir: str,
) -> None:
    sra_id = os.path.basename(sra_file)
    if sra_file.endswith(".sra"):
        (sra_id, _) = os.path.splitext(sra_id)
    sra_file: _FileArtifacts = _FileArtifacts(
        sra_file,
        working_id=sra_id,
        cwd=os.path.join(output_dir, sra_id),
        res_dir=csv_dir
    )
    if sra_file.result_exists():
        print(f"Result for {sra_file} exists in {sra_file.result_file()}")
        return None
    sra_file.create_cwd()
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
        vcf_files_1,
        vcf_files_2,
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


def prepare_amplicon_info(
    primer_bed_file: _FileArtifacts,
    amplicon_info_file: _FileArtifacts,
) -> _FileArtifacts:
    # This function is a direct copy from:
    # https://github.com/galaxyproject/tools-iuc/tree/master/tools/ivar/
    (bed_file, ) = primer_bed_file.file_path
    primer_starts = {}
    with open(bed_file) as i:
        for line in i:
            f = line.strip().split('\t')
            try:
                if f[5] == '+':
                    primer_starts[f[3]] = int(f[1])
                elif f[5] == '-':
                    primer_starts[f[3]] = int(f[2]) - 1
                else:
                    raise ValueError()
            except (IndexError, ValueError):
                sys.exit(
                    'Primer BED file needs to be TAB-separated with the '
                    'following columns: '
                    'chrom, chromStart, chromEnd, name, score, strand, '
                    'where "chromStart", "chromEnd" need to be integer values '
                    'and "strand" needs to be either "+" or "-".'
                )
    (info_file, ) = amplicon_info_file.file_path
    with open(info_file) as i:
        ret_lines = []
        for line in i:
            first = last = None
            for pname in line.strip().split('\t'):
                try:
                    primer_start = primer_starts[pname]
                except KeyError:
                    sys.exit(
                        'Amplicon info with primer name not found in '
                        f'primer BED file: "{pname}"'
                    )
                if first is None or primer_start < primer_starts[first]:
                    first = pname
                if last is None or primer_start > primer_starts[last]:
                    last = pname
            if first == last:
                sys.exit(
                    line
                    + 'is not a proper amplicon info line.'
                )
            ret_lines.append(f'{first}\t{last}\n')
    res = amplicon_info_file.coupled_files(
        ".tsv",
        exec_name="prepare_amplicon_info"
    )
    (res_file, ) = res.file_path
    with open(res_file, "w") as o:
        o.writelines(ret_lines)
    return res


def write_amplicon_info_file(primer_bed_file: _FileArtifacts) -> _FileArtifacts:
    # This function is a direct copy from:
    # https://github.com/galaxyproject/tools-iuc/tree/master/tools/ivar/
    AMPLICON_PAT = re.compile(r'.*_(?P<num>\d+).*_(?P<name>L(?:EFT)?|R(?:IGHT)?)')
    (primer_file, ) = primer_bed_file.file_path
    with open(primer_file) as bed_file:
        amplicon_sets = {}
        for line in bed_file:
            fields = line.strip().split('\t')
            start = int(fields[1])
            name = fields[3]
            re_match = AMPLICON_PAT.match(name)
            if re_match is None:
                raise ValueError(
                    '{} does not match expected amplicon name format'.format(name)
                )
            amplicon_id = int(re_match.group('num'))
            amplicon_set = amplicon_sets.get(amplicon_id, [])
            amplicon_set.append((name, start))
            amplicon_sets[amplicon_id] = amplicon_set

    # write amplicons sorted by number with primers sorted by start position
    res = primer_bed_file.coupled_files(
        ".tsv",
        exec_name="write_amplicon_info_file"
    )
    (tsv_file, ) = res.file_path
    with open(tsv_file, "w") as amplicon_info_file:
        for id in sorted(amplicon_sets):
            amplicon_info = '\t'.join(
                [name for name, start in sorted(
                    amplicon_sets[id], key=lambda x: x[1]
                )]
            ) + '\n'
            amplicon_info_file.write(amplicon_info)
    return res


def completemask(ivar_getmasked: IvarGetMaskedWrapper) -> _FileArtifacts:
    # This function is a direct copy from:
    # https://github.com/galaxyproject/tools-iuc/tree/master/tools/ivar/
    (masked_primers_file, ) = ivar_getmasked.output_files.file_path
    with open(masked_primers_file) as i:
        getmasked_output = i.readline().strip()

    (tsv_file, ) = ivar_getmasked.amplicon_info_file.file_path
    if not getmasked_output:
        pass
        # print()
        # print('No affected primer binding sites found!')
    else:
        masked_primers = getmasked_output.split('\t')
        with open(tsv_file) as i:
            amplicon_data = [line.strip().split('\t') for line in i]

        masked_complete = []
        for primer in masked_primers:
            for amplicon in amplicon_data:
                if primer in amplicon:
                    masked_complete += amplicon
        result = '\t'.join(sorted(set(masked_complete)))
        # print()
        # print('Removing reads primed with any of:')
        # print(result)
        with open(masked_primers_file, 'w') as o:
            o.write(result + '\n')
    return ivar_getmasked.output_files


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="SRA to CSV pipeline")
    parser.add_argument("-v", "--version", action="version",
                        version=f"sra2variant {__version__}")
    parser.add_argument("-r", "--refSeq", required=True,
                        help="File path of reference genome")
    parser.add_argument("-p", "--primers", required=True,
                        help="ARTIC primer bed files")
    parser.add_argument("-a", "--amplicon", required=True,
                        help="ARTIC primers to amplicon assignments")
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
    primer_file = args.primers
    amplicon_info_file = args.amplicon
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

    run_workflow(
        refSeq_file=refSeq_file,
        primer_file=primer_file,
        amplicon_info_file=amplicon_info_file,
        sra_files=sra_files,
        output_dir=output_dir,
        csv_dir=csv_dir,
        cores=cores,
        threads=threads
    )
