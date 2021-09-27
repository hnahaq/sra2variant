# A pipline to extract variant in SRA of SARS-CoV-2


## Installation

Use the following command to install

```bash
conda env update --file https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/environment.yml
```

## Quickstart

1. Download the reference genome in `fasta` format and store as `/path/to/SARS-CoV-2_refSeq.fasta`.
2. Download reads files in `sra` format and store them in `/path/to/sra_dir/`.

To use pipeline for WGS paired end reads
```bash
conda activate sra2variant

sra2variant-WGS-PE -i /path/to/sra_dir/ -r /path/to/SARS-CoV-2_refSeq.fasta
```

Other pipelines are under development
