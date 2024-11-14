import argparse

parser = argparse.ArgumentParser(
    description = "Ribo Remover\nUses NCBI blastn to compare fastq files to a ribo database of mammalian ribo and remove any that match"
)
parser.add_argument("--input_fastqs", help="filenames of fastq files to filter. If ends in .gz, then assumed to be gzipped. If more than one fastq, all are assumed to have the same reads (such as for paired ends) in the same order and are filtered out if any read matches the database", nargs="+")
parser.add_argument("--output_fastqs", help="filenames to write filtered fastq to. If ends in .gz, will output as gzipped. One filename per input fastq", nargs="+")
parser.add_argument("--ribo_db", help="blast db of ribo to filter against", default="blast_db/ribodb")
parser.add_argument("--num_threads", help="number of threads to use per blastn instance (one for each input fastq). Additional threads are also used for gzip and other aspects.", default=1, type=int)

args = parser.parse_args()

assert len(args.input_fastqs) == len(args.output_fastqs), "Must have same number of input fastqs as output fastqs"


import gzip
import subprocess
import time
import contextlib
import sys

RIBO_DB = args.ribo_db

# TODO: package the blast executable
#BLAST_DIR = "/project/itmatlab/SOFTWARE/PORT/PORT-0.8.5e-beta/norm_scripts/ncbi-blast-2.2.30+"
BLASTN_EXE = "./blastn"

# CODE TO MAKE THE RIBO DB:
#makeblastdb -dbtype nucl -in {RIBO_SEQUENCES} -out {RIBO_DB}
E_VALUE_THRESHOLD = 1e-7

def maybe_gzip_open(filename, *args):
    if filename.endswith(".gz"):
        return gzip.open(filename, *args, compresslevel=6)
    else:
        return open(filename, *args)

# Start blastn on each input fastq
blast_procs = []
for input_fastq in args.input_fastqs:
    CAT = "zcat" if input_fastq.endswith(".gz") else "cat"
    cmd = f"{CAT} {input_fastq} | " \
        "sed -n '1~4s/^@/>/p;2~4p' | " \
        f"{BLASTN_EXE} -task blastn -db {RIBO_DB} -query - -outfmt '10 qseqid' -evalue {E_VALUE_THRESHOLD} -num_threads {args.num_threads}"
    print(cmd, file=sys.stderr)

    blast = subprocess.Popen(
        cmd,
        shell=True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
    )
    blast_procs.append(blast)

start_time = time.time()
num_filtered = 0
num_unfiltered = 0
filtered_ids: list[None | str] = [None for _ in blast_procs]
DONE = "!!!DONE!!!"
with contextlib.ExitStack() as stack:
    # Open input fastqs
    in_fastqs = [stack.enter_context(maybe_gzip_open(input_fastq, "rb")) for input_fastq in args.input_fastqs]
    out_fastqs = [stack.enter_context(maybe_gzip_open(output_fastq, "wb")) for output_fastq in args.output_fastqs]

    while True:
        # Read one read from each input fastq
        header_lines = [in_fastq.readline() for in_fastq in in_fastqs]
        if all(not header for header in header_lines):
            # We're done reading the files
            break
        seq_lines = [in_fastq.readline() for in_fastq in in_fastqs]
        separator_lines = [in_fastq.readline() for in_fastq in in_fastqs]
        quality_lines = [in_fastq.readline() for in_fastq in in_fastqs]

        # Extract read id
        ids = [header_line.removeprefix(b"@").split(b" ", maxsplit=1)[0]
                for header_line in header_lines]
        assert len(set(ids)) == 1, f"FASTQ files did not have matching read ids at: {', '.join(h.decode() for h in header_lines)}"
        id = ids[0]

        # Check if previous line from blastn outputs match this read
        filtered = any(filtered_id == id for filtered_id in filtered_ids)

        # Advance to next output of each blastn process
        for idx, blast in enumerate(blast_procs):
            while filtered_ids[idx] != DONE:
                blastout = blast.stdout.readline()

                if not blastout:
                    # We have finished with all blast output, everything else is NOT filtered
                    filtered_ids[idx] = DONE
                    # Join on the blast process, which should now be done, to check for blast errors
                    _, err = blast.communicate()
                    if blast.returncode != 0:
                        raise Exception(f"BLAST Failed with return code {blast.returncode}:\n{err.decode()}")
                    break

                filtered_ids[idx] = blastout.strip()
                if filtered_ids[idx] == id:
                    filtered = True
                else:
                    # Any ID other than our current one indicates BLAST is ahead of our position
                    # and so we need to stop looking to see if this is filtered
                    break

        if filtered:
            num_filtered += 1
        else:
            num_unfiltered += 1
            # Output to each fastq
            for (header, seq, sep, qual, out) in zip(header_lines, seq_lines, separator_lines, quality_lines, out_fastqs):
                out.write(header)
                out.write(seq)
                out.write(sep)
                out.write(qual)

end_time = time.time()
print(f"Out of {num_filtered + num_unfiltered} total reads, {num_filtered} ribo found and {num_unfiltered} remain", file=sys.stderr)
print(f"Done in {end_time - start_time:0.1f} seconds", file=sys.stderr)
