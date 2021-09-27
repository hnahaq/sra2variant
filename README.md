# A pipline to extract variant in SRA of SARS-CoV-2


## Installation

Use the following command to install

```bash
conda env update --file environment.yml
```

## Quickstart

Download the reference genome in `fasta` format and `sra` files.

To use pipeline for WGS paired end reads
```bash
conda activate sra2variant

sra2variant-WGS-PE -i /path/to/sra_dir/ -r ./refSeq/SARS-CoV-2_refSeq.fasta
```

Other pipelines are under development
