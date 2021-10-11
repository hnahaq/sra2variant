# A pipline to detect variant in sequencing data of SARS-CoV-2


## Installation

[Bioconda](https://bioconda.github.io/user/install.html#install-conda) is required to install the tool

Use the following command to install and activate the environment

```bash
conda env update --file https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/environment.yml
conda activate sra2variant
```

## Quickstart

### WGS PE
1. First a reference genome in `fasta` format is needed. The following command downloads and store the genome as `reference/NC_045512.2.fasta`.

```bash
mkdir ./reference
wget https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/sra2variant/data/NC_045512.2.fasta
mv NC_045512.2.fasta ./reference
```

2. Download reads files in `sra` format and store them in a separate directory. Here two `sra` files are stored in `./wgs_reads` directory.

```bash
mkdir ./wgs_reads
prefetch -o ./wgs_reads/SRR14119630.sra SRR14119630
prefetch -o ./wgs_reads/SRR14119629.sra SRR14119629
```

3. Use the pipeline for WGS paired end reads. In this example, We use the reference genome `./reference/NC_045512.2.fasta` to analyze all `sra` files in `./wgs_reads` directory.

```bash
sra2variant-WGS-PE -i ./wgs_reads -r ./reference/NC_045512.2.fasta
```

### ARTIC PE

1. First a reference genome in `fasta` format, artic primer in `bed` format and a amplicon assignment in `tsv` format are needed. The following command downloads and store the files in `reference`.

```bash
mkdir ./reference

wget https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/sra2variant/data/NC_045512.2.fasta
mv NC_045512.2.fasta ./reference

wget https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/sra2variant/data/ARTIC_nCoV-2019_v3.bed
mv ARTIC_nCoV-2019_v3.bed ./reference

wget https://raw.githubusercontent.com/wuaipinglab/sra2variant/main/sra2variant/data/ARTIC_amplicon_info_v3.tsv
mv ARTIC_amplicon_info_v3.tsv ./reference
```

2. Download reads files in `sra` format and store them in a separate directory. Here two `sra` files are stored in `./artic_reads` directory.

```bash
mkdir ./artic_reads
prefetch -o ./artic_reads/SRR14388832.sra SRR14388832
prefetch -o ./artic_reads/SRR14398873.sra SRR14398873
```

3. Use the pipeline for WGS paired end reads. In this example, We use the reference genome `./reference/NC_045512.2.fasta` to analyze all `sra` files in `./artic_reads` directory.

```bash
sra2variant-ARTIC-PE -i ./artic_reads/ \
                     -r ./reference/NC_045512.2.fasta \
                     -p reference/ARTIC_nCoV-2019_v3.bed \
                     -a reference/ARTIC_amplicon_info_v3.tsv
```

Other pipelines are under development
