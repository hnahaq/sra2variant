# A pipline to detect variant in sequencing data of SARS-CoV-2


## Installation

[Bioconda](https://bioconda.github.io/user/install.html#install-conda) is required to install the tool

Use the following command to install and activate the environment

```bash
conda env update --file https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/environment.yml
conda activate sra2variant
```

## Quickstart

1. First a reference genome in `fasta` format is needed.

```bash
mkdir ./reference
wget https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/sra2variant/data/NC_045512.2.fasta
mv NC_045512.2.fasta ./reference
```

2. Download reads files in `sra` format and store them in a separate directory.

```bash
mkdir ./sample_reads
prefetch -o ./sample_reads/ERR4989943.sra ERR4989943
```

3. Use the pipeline for WGS paired end reads
```bash
sra2variant-WGS-PE -i ./sample_reads -r ./reference/NC_045512.2.fasta
```

Other pipelines are under development
