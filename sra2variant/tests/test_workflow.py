import os
import glob
import shutil
import unittest

from sra2variant import WGS_PE
from sra2variant.pipeline.sra2fastq import PrefetchWrapper

THIS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(THIS_DIR, os.pardir, "data")
TEMP_DIR = os.path.join(THIS_DIR, os.pardir, "temp")


class TestWorkFlow(unittest.TestCase):
    def setUp(self):
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)
        self.refSeq_file = os.path.join(DATA_DIR, "NC_045512.2.fasta")
        self.output_dir = os.path.join(TEMP_DIR, "output")
        self.csv_dir = os.path.join(TEMP_DIR, "csv")
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.csv_dir):
            os.makedirs(self.csv_dir)

    # def test_SE_workflow(self):
    #     sra_id = "ERR4238190"

    def test_PE_workflow(self):
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
            threads=os.cpu_count()
        )
        
        os.remove(sra_file)

    def tearDown(self):
        for fn in glob.glob(os.path.join(DATA_DIR, "NC_045512.2.fasta.*")):
            os.remove(fn)
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    unittest.main()
