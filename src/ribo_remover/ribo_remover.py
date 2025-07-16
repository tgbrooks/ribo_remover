import gzip
import subprocess
import time
import contextlib
import sys
from typing import Union
from importlib import resources

# CODE TO MAKE THE RIBO DB:
# makeblastdb -dbtype nucl -in {RIBO_SEQUENCES} -out {RIBO_DB}

E_VALUE_THRESHOLD = 1e-7


def maybe_gzip_open(filename: str, *args):
    if filename.endswith(".gz"):
        return gzip.open(filename, *args, compresslevel=6)
    else:
        return open(filename, *args)


def ribo_remover(
    input_fastqs: [str],
    output_fastqs: [str],
    stats_file: Union[str, None] = None,
    num_threads: int = 1,
):
    assert len(input_fastqs) == len(
        output_fastqs
    ), "Must have same number of input fastqs as output fastqs"

    # We open the stats_file now so that the script fails immediately if it's not openable
    stats_file = open(stats_file, "wt") if stats_file is not None else None

    # Start blastn on each input fastq
    blast_procs = []
    for input_fastq in input_fastqs:
        CAT = "zcat" if input_fastq.endswith(".gz") else "cat"
        with resources.as_file(
            resources.files("ribo_remover.data").joinpath("blastn")
        ) as blastn_exe, resources.as_file(
            resources.files("ribo_remover.data").joinpath("blast_db")
        ) as ribo_db:
            cmd = (
                f"{CAT} {input_fastq} | "
                "fastq-converter | "
                f"{blastn_exe} -task blastn -db {ribo_db}/ribodb -query - -outfmt '10 qseqid' -evalue {E_VALUE_THRESHOLD} -num_threads {num_threads} -num_alignments 1"
            )

            blast = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            blast_procs.append(blast)

    start_time = time.time()
    num_filtered = 0
    num_unfiltered = 0
    read_num = 0
    filtered_read_nums: list[None | str] = [None for _ in blast_procs]
    DONE = "!!!DONE!!!"
    with contextlib.ExitStack() as stack:
        # Open input fastqs
        in_fastqs = [
            stack.enter_context(maybe_gzip_open(input_fastq, "rb"))
            for input_fastq in input_fastqs
        ]
        out_fastqs = [
            stack.enter_context(maybe_gzip_open(output_fastq, "wb"))
            for output_fastq in output_fastqs
        ]

        while True:
            # Read one read from each input fastq
            header_lines = [in_fastq.readline() for in_fastq in in_fastqs]
            if all(not header for header in header_lines):
                # We're done reading the files
                break
            seq_lines = [in_fastq.readline() for in_fastq in in_fastqs]
            separator_lines = [in_fastq.readline() for in_fastq in in_fastqs]
            quality_lines = [in_fastq.readline() for in_fastq in in_fastqs]
            read_num += 1

            # Extract read id
            ids = [
                header_line.removeprefix(b"@").split(b" ", maxsplit=1)[0].strip()
                for header_line in header_lines
            ]
            assert (
                len(set(ids)) == 1
            ), f"FASTQ files did not have matching read ids at: {', '.join(h.decode() for h in header_lines)}"

            # Check if previous line from blastn outputs match this read
            filtered = any(
                filtered_read_num == read_num
                for filtered_read_num in filtered_read_nums
            )

            # Advance to next output of each blastn process
            for idx, blast in enumerate(blast_procs):
                while (
                    filtered_read_nums[idx] is None
                    or filtered_read_nums[idx] == read_num
                ):
                    blastout = blast.stdout.readline()

                    if not blastout:
                        # We have finished with all blast output, everything else is NOT filtered
                        filtered_read_nums[idx] = DONE
                        # Join on the blast process, which should now be done, to check for blast errors
                        _, err = blast.communicate()
                        if blast.returncode != 0:
                            raise Exception(
                                f"BLAST Failed with return code {blast.returncode}:\n{err.decode()}"
                            )
                        break

                    filtered_read_nums[idx] = int(blastout.strip())
                    if filtered_read_nums[idx] == read_num:
                        filtered = True

            if filtered:
                num_filtered += 1
            else:
                num_unfiltered += 1
                # Output to each fastq
                for header, seq, sep, qual, out in zip(
                    header_lines, seq_lines, separator_lines, quality_lines, out_fastqs
                ):
                    out.write(header)
                    out.write(seq)
                    out.write(sep)
                    out.write(qual)

    end_time = time.time()
    total = num_filtered + num_unfiltered
    pct_filtered = num_filtered / total
    pct_unfiltered = num_unfiltered / total
    print(
        f"Out of {total} total reads, {num_filtered} ({pct_filtered:0.1%}) ribo found and {num_unfiltered} ({pct_unfiltered:0.1%}) remain",
        file=sys.stderr,
    )
    print(f"Done in {end_time - start_time:0.1f} seconds", file=sys.stderr)
    if stats_file is not None:
        stats_file.write("class,number,percent,files\n")
        stats_file.write(
            f"filtered,{num_filtered},{pct_filtered:0.1%},{' '.join(input_fastqs)}\n"
        )
        stats_file.write(
            f"unfiltered,{num_unfiltered},{pct_unfiltered:0.1%},{' '.join(input_fastqs)}\n"
        )
        stats_file.close()
