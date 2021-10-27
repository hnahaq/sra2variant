import os
import glob
import shutil
import unittest

from sra2variant import WGS_PE, ARTIC_PE
from sra2variant.pipeline.sra2fastq import PrefetchWrapper

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
        self.output_dir = os.path.join(TEMP_DIR, "output")
        self.csv_dir = os.path.join(TEMP_DIR, "csv")
        self.error_dir = os.path.join(TEMP_DIR, "errors")
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.csv_dir):
            os.makedirs(self.csv_dir)
        if not os.path.exists(self.error_dir):
            os.makedirs(self.error_dir)

    def test_ARTIC_PE(self):
        sra_id = "SRR14388832"
        sra_file = os.path.join(DATA_DIR, f"{sra_id}.sra")
        if not os.path.exists(sra_file):
            prefetch = PrefetchWrapper(sra_id, DATA_DIR)
            prefetch.execute_cmd()

        ARTIC_PE.run_workflow(
            refSeq_file=self.refSeq_file,
            primer_file=self.primer_file,
            amplicon_info_file=self.amplicon_info_file,
            sra_files=[sra_file],
            output_dir=self.output_dir,
            csv_dir=self.csv_dir,
            cores=2,
            threads=os.cpu_count(),
            error_dir=self.error_dir,
        )
        for fn in glob.glob(os.path.join(DATA_DIR, "_*")):
            os.remove(fn)
        os.remove(sra_file)

    def test_WGS_PE(self):
        sra_id = "ERR4989943"
        sra_file = os.path.join(DATA_DIR, f"{sra_id}.sra")
        if not os.path.exists(sra_file):
            prefetch = PrefetchWrapper(sra_id, DATA_DIR)
            prefetch.execute_cmd()

        WGS_PE.run_workflow(
            self.refSeq_file,
            [sra_file],
            self.output_dir,
            self.csv_dir,
            cores=2,
            threads=os.cpu_count(),
            error_dir=self.error_dir,
        )
        os.remove(sra_file)

    def test_error_tolerance(self):
        sra_id = "SRR16026846"
        sra_file = os.path.join(DATA_DIR, f"{sra_id}.sra")
        if not os.path.exists(sra_file):
            prefetch = PrefetchWrapper(sra_id, DATA_DIR)
            prefetch.execute_cmd()

        ARTIC_PE.run_workflow(
            refSeq_file=self.refSeq_file,
            primer_file=self.primer_file,
            amplicon_info_file=self.amplicon_info_file,
            sra_files=[sra_file],
            output_dir=self.output_dir,
            csv_dir=self.csv_dir,
            cores=2,
            threads=os.cpu_count(),
            error_dir=self.error_dir,
            max_errors=1
        )
        for fn in glob.glob(os.path.join(DATA_DIR, "_*")):
            os.remove(fn)
        os.remove(sra_file)

    def tearDown(self):
        for fn in glob.glob(os.path.join(DATA_DIR, "NC_045512.2.fasta.*")):
            os.remove(fn)
        for fn in glob.glob(os.path.join(DATA_DIR, "*log.txt")):
            os.remove(fn)
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    unittest.main()
