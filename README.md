# A pipline to extract variant in SRA of SARS-CoV-2


## Installation

Download the repository from GitHub, `unzip` the file and `cd` into the directory. Use the following command to install

```bash
conda env update --file environment.yml
pip install .
```

## Quickstart

Download the reference genome in `fasta` format and `sra` files.

To use pipeline for WGS paired end reads
```bash
sra2variant-WGS-PE -i /path/to/sra_dir/ -r ./refSeq/SARS-CoV-2_refSeq.fasta
```

Other pipelines are under development
