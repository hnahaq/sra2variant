import os
import glob
import shutil
import unittest
from functools import partial
from multiprocessing import Pool, cpu_count

from sra2variant import WGS_PE
from sra2variant.pipeline.cmd_wrapper import CMDwrapperBase
from sra2variant.pipeline.sra2fastq import PrefetchWrapper
from sra2variant.pipeline.reads_mapping import (
    BWAindexWrapper,
    BWAmemWrapper,
)

THIS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(THIS_DIR, os.pardir, "data")
TEMP_DIR = os.path.join(THIS_DIR, os.pardir, "temp")


class TestWorkFlow(unittest.TestCase):
    def setUp(self):
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)
        self.refSeq_file = os.path.join(DATA_DIR, "SARS-CoV-2_refSeq.fasta")
        self.output_dir = os.path.join(TEMP_DIR, "output")
        self.csv_dir = os.path.join(TEMP_DIR, "csv")
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.csv_dir):
            os.makedirs(self.csv_dir)

    # def test_SE_workflow(self):
    #     sra_id = "ERR4238190"
    #     sra_file = _FileArtifacts(
    #         file_path=(os.path.join(DATA_DIR, f"{sra_id}.sra"), ),
    #         cwd=os.path.join(self.output_dir, sra_id),
    #         working_id=sra_id,
    #         res_dir=self.csv_dir
    #     )
    #     if sra_file.exist():
    #         fq_dump = FastqDumpWrapper(sra_file)
    #         fastq_files = fq_dump.execute_cmd()
    #         self.assertTrue(fastq_files.exist(),
    #                         "SE fastq-dump failed")

    #         TrimmomaticWrapper.set_exec_args(("SLIDINGWINDOW:4:0", ))
    #         trimmomatic = TrimmomaticWrapper(fastq_files)
    #         fastq_files = trimmomatic.execute_cmd()
    #         self.assertTrue(fastq_files.exist(),
    #                         "SE trimmomatic failed")

    #         refSeq_file = _FileArtifacts(
    #             (os.path.join(DATA_DIR, "SARS-CoV-2_refSeq.fasta"),),
    #             cwd=DATA_DIR,
    #             working_id="SARS-CoV-2_refSeq"
    #         )

    #         ht2_build = Hisat2BuildWrapper(refSeq_file)
    #         ht2_build.execute_cmd()
    #         Hisat2Wrapper.set_exec_args(("--no-spliced-alignment", "--no-unal"))
    #         Hisat2Wrapper.set_base_index(ht2_build)
    #         hisat2 = Hisat2Wrapper(fastq_files)
    #         bam_files = hisat2.execute_cmd()

    #         mpileup = BCFtoolsMpileupWrapper(bam_files, refSeq_file)
    #         print(mpileup)
    #         mpileup_file = mpileup.execute_cmd()
    #         self.assertIsNotNone(mpileup_file, "PE BCFtools mpileup failed")
    #         bcf_norm = BCFtoolsNormWrapper(mpileup_file, refSeq_file)
    #         print(bcf_norm)
    #         vcf_file = bcf_norm.execute_cmd()
    #         self.assertIsNotNone(vcf_file, "PE BCFtools norm failed")

    #         bt2_build = Bowtie2BuildWrapper(refSeq_file)
    #         bt2_build.execute_cmd()
    #         Bowtie2Wrapper.set_base_index(bt2_build)
    #         bowtie2 = Bowtie2Wrapper(fastq_files)
    #         bam_files = bowtie2.execute_cmd()


    #         lofreq_IQ = LoFreqIndelQualWrapper(bam_files, refSeq_file)
    #         print(lofreq_IQ)
    #         bam_files = lofreq_IQ.execute_cmd()
    #         self.assertIsNotNone(bam_files, "PE LoFreq indelqual failed")
    #         lofreq_call = LoFreqCallWrapper(bam_files, refSeq_file)
    #         print(lofreq_call)
    #         vcf_file = lofreq_call.execute_cmd()
    #         self.assertIsNotNone(vcf_file, "PE LoFreq call failed")

    def test_PE_workflow(self):
        CMDwrapperBase.set_threads(cpu_count())
        sra_id = "ERR4989943"
        sra_file = os.path.join(DATA_DIR, f"{sra_id}.sra")
        refSeq_file = os.path.join(DATA_DIR, "SARS-CoV-2_refSeq.fasta")
        if not os.path.exists(sra_file):
            prefetch = PrefetchWrapper(sra_id, DATA_DIR)
            print(prefetch)
            prefetch.execute_cmd()
        bwaindex = BWAindexWrapper(refSeq_file)
        BWAmemWrapper.set_base_index(bwaindex)
        workflow2 = partial(
            WGS_PE.workflow,
            output_dir=self.output_dir,
            csv_dir=self.csv_dir,
        )
        with Pool(2) as p:
            p.map(workflow2, [sra_file])
        os.remove(sra_file)


    def tearDown(self):
        for fn in glob.glob(os.path.join(DATA_DIR, "SARS-CoV-2_refSeq.fasta.*")):
            os.remove(fn)
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    unittest.main()
