# Ribo remover

This is a small utility to remove mammalian ribosomal RNA content from FASTQ files, particularly for RNA-seq.
Ribo-remover produces filtered FASTQ files without the rRNA content, using NCBI's `blastn` utility to identify rRNA.
It supports paired end reads.

## Installation

Using python 3.9 or greater, and preferably using a virtual environment, run the following:

``` shell
pip install git+https://github.com/tgbrooks/ribo_remover
```

## Usage

``` shell
# For single-end reads
ribo-remover --input-fastqs R1.fastq.gz --output-fastqs filtered.R1.fastq.gz

# For paired-end reads
ribo-remover --input-fastqs R1.fastq.gz R2.fastq.gz --output-fastqs filtered.R1.fastq.gz filtered.R2.fastq.gz

# To output a file containing ribo content stats as a CSV file
ribo-remover --input-fastqs R1.fastq.gz R2.fastq.gz --output-fastqs filtered.R1.fastq.gz filtered.R2.fastq.gz --stats-file ribostats.csv
```

Both input and output can be plain `fastq` files or gzipped files ending in `.gz`, which will be automatically decompressed (for input files) and compressed (for output files) based off the provided file names.
