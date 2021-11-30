import os
import glob
import shutil
import filecmp
import unittest

from sra2variant import WGS_PE, ARTIC_PE
from sra2variant.pipeline.sra2fastq import FastqDumpWrapper, PrefetchWrapper

THIS_DIR = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "sra2variant"
)
DATA_DIR = os.path.join(THIS_DIR, "data")
TEMP_DIR = os.path.join(THIS_DIR, "temp")


class TestWorkFlow(unittest.TestCase):
    def setUp(self):
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)
        self.refSeq_file = os.path.join(DATA_DIR, "NC_045512.2.fasta")
        self.primer_file = os.path.join(DATA_DIR, "ARTIC_nCoV-2019_v3.bed")
        self.amplicon_info_file = os.path.join(
            DATA_DIR, "ARTIC_amplicon_info_v3.tsv")
        # self.output_dir = os.path.join(TEMP_DIR, "output")
        # self.csv_dir = os.path.join(TEMP_DIR, "csv")
        # self.error_dir = os.path.join(TEMP_DIR, "errors")
        # if not os.path.exists(self.output_dir):
        #     os.makedirs(self.output_dir)
        # if not os.path.exists(self.csv_dir):
        #     os.makedirs(self.csv_dir)
        # if not os.path.exists(self.error_dir):
        #     os.makedirs(self.error_dir)

    def test_ARTIC_PE(self):
        sra_id = "SRR14388832"
        reads_dir = os.path.join(TEMP_DIR, "artic_reads")
        if not os.path.exists(reads_dir):
            os.makedirs(reads_dir)
        sra_file = PrefetchWrapper(sra_id, reads_dir).execute_cmd()
        FastqDumpWrapper(sra_file).execute_cmd()

        ARTIC_PE.main([
            "-r", self.refSeq_file,
            "-p", self.primer_file,
            "-a", self.amplicon_info_file,
            "-i", reads_dir,
            "-o", TEMP_DIR,
            "-c", "2",
            "-t", str(os.cpu_count()),
            "--maxErrors", "0"
        ])
        csv_file1 = os.path.join(TEMP_DIR, "csv", sra_id + ".csv")
        csv_file2 = os.path.join(TEMP_DIR, "csv", sra_id + "_2.csv")
        os.rename(csv_file1, csv_file2)
        ARTIC_PE.main_fastq([
            "-r", self.refSeq_file,
            "-p", self.primer_file,
            "-a", self.amplicon_info_file,
            "-i", reads_dir,
            "-o", TEMP_DIR,
            "-c", "2",
            "-t", str(os.cpu_count()),
            "--maxErrors", "0",
            "--forward", "_1.fastq",
            "--reverse", "_2.fastq"
        ])
        self.assertTrue(filecmp.cmp(csv_file1, csv_file2))
        shutil.rmtree(reads_dir)

    def test_WGS_PE(self):
        sra_id = "ERR4989943"
        reads_dir = os.path.join(TEMP_DIR, "wgs_reads")
        if not os.path.exists(reads_dir):
            os.makedirs(reads_dir)
        sra_file = PrefetchWrapper(sra_id, reads_dir).execute_cmd()
        FastqDumpWrapper(sra_file).execute_cmd()

        WGS_PE.main([
            "-r", self.refSeq_file,
            "-i", reads_dir,
            "-o", TEMP_DIR,
            "-c", "2",
            "-t", str(os.cpu_count()),
            "--maxErrors", "0"
        ])
        csv_file1 = os.path.join(TEMP_DIR, "csv", sra_id + ".csv")
        csv_file2 = os.path.join(TEMP_DIR, "csv", sra_id + "_2.csv")
        os.rename(csv_file1, csv_file2)
        WGS_PE.main_fastq([
            "-r", self.refSeq_file,
            "-i", reads_dir,
            "-o", TEMP_DIR,
            "-c", "2",
            "-t", str(os.cpu_count()),
            "--maxErrors", "0",
            "--forward", "_1.fastq",
            "--reverse", "_2.fastq"
        ])
        self.assertTrue(filecmp.cmp(csv_file1, csv_file2))
        shutil.rmtree(reads_dir)

    def test_error_tolerance(self):
        sra_id = "SRR16026846"
        reads_dir = os.path.join(TEMP_DIR, "artic_reads")
        if not os.path.exists(reads_dir):
            os.makedirs(reads_dir)
        PrefetchWrapper(sra_id, reads_dir).execute_cmd()

        ARTIC_PE.main([
            "-r", self.refSeq_file,
            "-p", self.primer_file,
            "-a", self.amplicon_info_file,
            "-i", reads_dir,
            "-o", TEMP_DIR,
            "-c", "2",
            "-t", str(os.cpu_count()),
            "--maxErrors", "1"
        ])
        for fn in glob.glob(os.path.join(DATA_DIR, "_*")):
            os.remove(fn)
        shutil.rmtree(reads_dir)

    def tearDown(self):
        for fn in glob.glob(os.path.join(DATA_DIR, "NC_045512.2.fasta.*")):
            os.remove(fn)
        for fn in glob.glob(os.path.join(DATA_DIR, "*log.txt")):
            os.remove(fn)
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    unittest.main()
