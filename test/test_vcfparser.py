import os
import glob
import shutil
import unittest

from sra2variant.pipeline.cmd_wrapper import _FileArtifacts
from sra2variant.vcfparser.vcf2csv import vcf2csv, AmbiguousSNPcombination

THIS_DIR = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "sra2variant"
)
DATA_DIR = os.path.join(THIS_DIR, "data")
TEMP_DIR = os.path.join(THIS_DIR, "temp")


class TestVCFparser(unittest.TestCase):

    def setUp(self):
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)
        os.mkdir(TEMP_DIR)

    def test_vcf2csv(self):
        for fp in glob.glob(os.path.join(DATA_DIR, "*.vcf")):
            fn = os.path.basename(fp)
            cwd = os.path.dirname(fp)
            (working_id, _) = os.path.splitext(fn)
            vcf_file = _FileArtifacts(
                fp,
                working_id=working_id,
                cwd=cwd,
                res_dir=TEMP_DIR
            )
            self.assertWarns(AmbiguousSNPcombination, vcf2csv, vcf_file)

    def tearDown(self):
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    unittest.main()
