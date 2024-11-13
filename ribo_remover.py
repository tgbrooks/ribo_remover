import argparse

parser = argparse.ArgumentParser(
    description = "Ribo Remover\nUsing NCBI blastn to compare fastq files to a ribo database and remove any that match"
)
parser.add_argument("--input_fastq", help="filename of fastq to filter. If ends in .gz, then assumed to be gzipped.")
parser.add_argument("--output_fastq", help="filename to write filtered fastq to. If ends in .gz, will output as gzipped")
parser.add_argument("--ribo_db", help="blast db of ribo to filter against", default="blast_db/ribodb")
parser.add_argument("--num_threads", help="number of threads to use for blast", default=1, type=int)

args = parser.parse_args()

import gzip
import subprocess
import time

RIBO_DB = args.ribo_db

# TODO: package the blast executable
BLAST_DIR = "/project/itmatlab/SOFTWARE/PORT/PORT-0.8.5e-beta/norm_scripts/ncbi-blast-2.2.30+"

#TODO: there is actually a more recent one of these
RIBO_SEQUENCES = "/project/itmatlab/SOFTWARE/PORT/PORT-0.8.5e-beta/norm_scripts/rRNA_mm9.fa"
#TODO: regenerate the ribo db
# CODE TO MAKE THE RIBO DB:
#f"{BLAST_DIR}/bin/makeblastdb -dbtype nucl -in {RIBO_SEQUENCES} -out {RIBO_DB}"
E_VALUE_THRESHOLD = 1e-7

def maybe_gzip_open(filename, *args):
    if filename.endswith(".gz"):
        return gzip.open(filename, *args, compresslevel=6)
    else:
        return open(filename, *args)

CAT = "zcat" if args.input_fastq.endswith(".gz") else "cat"
cmd = f"{CAT} {args.input_fastq} | " \
    "sed -n '1~4s/^@/>/p;2~4p' | " \
    f"{BLAST_DIR}/bin/blastn -task blastn -db {RIBO_DB} -query - -outfmt '10 qseqid' -evalue {E_VALUE_THRESHOLD} -num_threads {args.num_threads}"
print(cmd)

blast = subprocess.Popen(
    cmd,
    shell=True,
    stdout = subprocess.PIPE,
    stderr = subprocess.PIPE,
)

start_time = time.time()
num_filtered = 0
num_unfiltered = 0
filtered_id = None
DONE = "!!!DONE!!!"
with maybe_gzip_open(args.input_fastq, "rb") as in_fastq, maybe_gzip_open(args.output_fastq, "wb") as out_fastq:
    while True:
        header_line: str = in_fastq.readline()
        if not header_line:
            # We're done reading the file
            break
        seq_line = in_fastq.readline()
        separator_line = in_fastq.readline()
        quality_line = in_fastq.readline()

        id = header_line.removeprefix(b"@").split(b" ")[0]

        filtered = filtered_id == id
        while filtered_id != DONE:
            blastout = blast.stdout.readline()
            if not blastout:
                # We have finished with all blast output, everything else is NOT filtered
                filtered_id = DONE
                break
            filtered_id = blastout.strip()
            if filtered_id == id:
                filtered = True
            else:
                # Any ID other than our current one indicates BLAST is ahead of our position
                # and so we need to stop looking to see if this is filtered
                break
        if filtered:
            num_filtered += 1
            continue

        num_unfiltered += 1
        out_fastq.write(header_line)
        out_fastq.write(seq_line)
        out_fastq.write(separator_line)
        out_fastq.write(quality_line)

end_time = time.time()
print(f"Out of {num_filtered + num_unfiltered} total reads, {num_filtered} ribo found and {num_unfiltered} remain")
print(f"Done in {end_time - start_time:0.1f} seconds")
