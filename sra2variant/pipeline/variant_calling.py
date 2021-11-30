import os
import re
import sys

from sra2variant.pipeline.cmd_wrapper import (
    _FileArtifacts,
    CMDwrapperBase,
)


class LoFreqFaidxWrapper(CMDwrapperBase):

    exec_name: str = "lofreq faidx"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        self.output_files = None
        (refSeq_file, ) = input_files.file_path
        super().__init__(
            input_files,
            *args,
            refSeq_file
        )

    def _post_execution(self) -> None:
        self.output_files = self.input_files


class LoFreqViterbiWraper(CMDwrapperBase):

    exec_name: str = "lofreq viterbi"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        (sort_bam_file, ) = input_files.path_from_cwd()
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name="viterbi"
        )
        (viterbi_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            *args,
            "--ref", os.path.abspath(refSeq_file),
            "--out", viterbi_file,
            sort_bam_file,
        )

    def _post_execution(self) -> None:
        pass


class LoFreqIndelQualWrapper(CMDwrapperBase):

    exec_name: str = "lofreq indelqual"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name="indelqual"
        )
        (sort_bam_file, ) = input_files.path_from_cwd()
        (indelqual_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            *args,
            "--ref", os.path.abspath(refSeq_file),
            "--out", indelqual_file,
            sort_bam_file,
        )

    def _post_execution(self) -> None:
        pass


class LoFreqCallWrapper(CMDwrapperBase):

    exec_name: str = "lofreq call-parallel"

    def __init__(
        self,
        input_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.output_files = input_files.coupled_files(
            ".vcf",
            exec_name="lofreq call"
        )
        (indelqual_file, ) = input_files.path_from_cwd()
        (vcf_file, ) = self.output_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        super().__init__(
            input_files,
            "--pp-threads", self.threads,
            *args,
            "--verbose",
            "--ref", os.path.abspath(refSeq_file),
            "--out", vcf_file,
            indelqual_file
        )

    def _post_execution(self) -> None:
        pass


class LoFreqFilterWrapper(CMDwrapperBase):

    exec_name: str = "lofreq filter"

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        self.output_files = input_files.coupled_files(
            ".vcf",
            exec_name=self.exec_name
        )
        (vcf_file, ) = input_files.path_from_cwd()
        (filter_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-i", vcf_file,
            "-o", filter_file,
        )

    def _post_execution(self) -> None:
        pass


class IvarTrimWrapper(CMDwrapperBase):

    exec_name = "ivar trim"
    primer_bed_file: _FileArtifacts = None
    amplicon_info_file: _FileArtifacts = None

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (sort_bam_file, ) = input_files.path_from_cwd()
        bed_file_flag = ()
        if self.primer_bed_file:
            (bed_file, ) = self.primer_bed_file.file_path
            bed_file_flag = ("-b", os.path.abspath(bed_file))
        amplicon_info_flag = ()
        if self.amplicon_info_file:
            (tsv_file, ) = self.amplicon_info_file.file_path
            amplicon_info_flag = ("-f", os.path.abspath(tsv_file))
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (trimmed_bam_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            *bed_file_flag,
            *amplicon_info_flag,
            "-i", sort_bam_file,
            "-p", trimmed_bam_file[:-4]
        )

    @classmethod
    def set_primer_file(cls, primer_bed_file: _FileArtifacts) -> None:
        cls.primer_bed_file = primer_bed_file

    @classmethod
    def set_amplicon_info(cls, amplicon_info_file: _FileArtifacts) -> None:
        cls.amplicon_info_file = amplicon_info_file

    def _post_execution(self) -> None:
        pass


class IvarGetMaskedWrapper(CMDwrapperBase):

    exec_name = "ivar getmasked"
    primer_bed_file: _FileArtifacts = None
    amplicon_info_file: _FileArtifacts = None

    def __init__(self, input_files: _FileArtifacts, *args: str) -> None:
        (vcf_file, ) = input_files.path_from_cwd()
        (bed_file, ) = self.primer_bed_file.file_path
        (tsv_file, ) = self.amplicon_info_file.file_path
        self.output_files = input_files.coupled_files(
            ".txt",
            exec_name=self.exec_name
        )
        (masked_primers, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-b", os.path.abspath(bed_file),
            "-f", os.path.abspath(tsv_file),
            "-i", vcf_file,
            "-p", masked_primers[:-4]
        )

    @classmethod
    def set_primer_file(cls, primer_bed_file: _FileArtifacts) -> None:
        cls.primer_bed_file = primer_bed_file

    @classmethod
    def set_amplicon_info(cls, amplicon_info_file: _FileArtifacts) -> None:
        cls.amplicon_info_file = amplicon_info_file

    def _post_execution(self) -> None:
        pass


class IvarRemoveReadsWrapper(CMDwrapperBase):

    exec_name = "ivar removereads"
    primer_bed_file: _FileArtifacts = None

    def __init__(
        self,
        input_files: _FileArtifacts,
        masked_primers_file: _FileArtifacts,
        *args: str
    ) -> None:
        (sort_bam_file, ) = input_files.file_path
        (masked_primers, ) = masked_primers_file.file_path
        (bed_file, ) = self.primer_bed_file.file_path
        self.output_files = input_files.coupled_files(
            ".bam",
            exec_name=self.exec_name
        )
        (trimmed_bam_file, ) = self.output_files.path_from_cwd()
        super().__init__(
            input_files,
            *args,
            "-t", os.path.abspath(masked_primers),
            "-b", os.path.abspath(bed_file),
            "-i", sort_bam_file,
            "-p", trimmed_bam_file[:-4]
        )

    @classmethod
    def set_primer_file(cls, primer_bed_file: _FileArtifacts) -> None:
        cls.primer_bed_file = primer_bed_file

    def _post_execution(self) -> None:
        pass


class VCFintersectWrapper(CMDwrapperBase):

    exec_name = "vcfintersect"

    def __init__(
        self,
        input_files: _FileArtifacts,
        vcf_files: _FileArtifacts,
        refSeq_files: _FileArtifacts,
        union: bool,
        *args: str
    ) -> None:
        (vcf_file_1, ) = input_files.path_from_cwd()
        (vcf_file_2, ) = vcf_files.path_from_cwd()
        (refSeq_file, ) = refSeq_files.file_path
        self.output_files = vcf_files.coupled_files(
            ".vcf",
            exec_name="union" if union else "intersect"
        )
        super().__init__(
            input_files,
            *args,
            "-r", os.path.abspath(refSeq_file),
            "-u" if union else "-i", vcf_file_1,
            vcf_file_2
        )

    def _post_execution(self) -> None:
        (vcf_file, ) = self.output_files.file_path
        with open(vcf_file, "w") as f:
            f.writelines(self.stdout)


def sanitize_bed(primer_bed_file: _FileArtifacts) -> None:
    # This function is a direct copy from:
    # https://github.com/galaxyproject/tools-iuc/tree/master/tools/ivar/
    (bed_file, ) = primer_bed_file.file_path
    with open(bed_file) as i:
        bed_data = i.readlines()

    res = primer_bed_file.coupled_files(
        ".bed",
        exec_name="sanitize_bed"
    )
    (sanitize_bed_file, ) = res.file_path
    sanitized_data = []
    try:
        for record in bed_data:
            fields = record.split('\t')
            sanitized_data.append(
                '\t'.join(fields[:4] + ['60'] + fields[5:])
            )
    except IndexError:
        pass  # leave column number issue to getmasked
    else:
        with open(sanitize_bed_file, 'w') as o:
            o.writelines(sanitized_data)
    return res


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
    AMPLICON_PAT = re.compile(
        r'.*_(?P<num>\d+).*_(?P<name>L(?:EFT)?|R(?:IGHT)?)'
    )
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
                    f'{name} does not match expected amplicon name format'
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
