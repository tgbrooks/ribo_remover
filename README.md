# Ribo remover

This is a small utility to remove mammalian ribosomal RNA content from FASTQ files, particularly for RNA-seq.
Ribo-remover produces filtered FASTQ files without the rRNA content, using NCBI's `blastn` utility to identify rRNA, similar to [PORT](https://github.com/itmat/Normalization).
It supports single end or paired end reads.

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

# To also output a file containing ribo content stats as a CSV file
ribo-remover --input-fastqs R1.fastq.gz R2.fastq.gz --output-fastqs filtered.R1.fastq.gz filtered.R2.fastq.gz --stats-file ribostats.csv
```

Both input and output can be plain `fastq` files or gzipped files ending in `.gz`, which will be automatically decompressed (for input files) and compressed (for output files) based off the provided file names.

## Multi-threading

You can request multi-threading with the `--threads` option, but note that this is the number of threads used _per fastq file_.
If you have paired end reads, then two instances of `blastn` will be run, each using the provided number of threads.
In addition, the input/output files may need computational resources to be gzipped.
This means that you generally want to have at least `2n+1` cores available when requesting `--threads n` with paired ends and `n+1` cores if using single ends.
Performance benefits also drop off after `--threads 2`.

## UMI tags in separate files

Sometimes, unique molecule identifier (UMI) tags are provided as separate files, for example `R1.fastq`, `R2.fastq` along side UMI tag files `I1.fastq` and `I2.fastq`.
In this case, any ribosomal reads will want to be removed from the UMI tag files as well.
This should be possible by passing all fastq files (including UMIs) together to `--input-fastqs` and giving corresponding outputs for `--output-fastqs`.
This will BLAST each UMI tag against the ribosomal database as well, but typical UMI tags are too short to give a significant hit and so should have no real effect on the ribosomal removal.
However, this functionality has not been extensively tested.
