def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Ribo Remover\nUses NCBI blastn to compare fastq files to a ribo database of mammalian ribo and remove any that match"
    )
    parser.add_argument(
        "--input-fastqs",
        help="filenames of fastq files to filter. If ends in .gz, then assumed to be gzipped. If more than one fastq, all are assumed to have the same reads (such as for paired ends) in the same order and are filtered out if any read matches the database",
        nargs="+",
    )
    parser.add_argument(
        "--output-fastqs",
        help="filenames to write filtered fastq to. If ends in .gz, will output as gzipped. One filename per input fastq",
        nargs="+",
    )
    # TODO: support alternative ribo dbs
    # parser.add_argument("--ribo_db", help="blast db of ribo to filter against. Construct with `makeblastdb -dbtype nucl -in my_sequences.fa -out my_db_name`", default="blast_db/ribodb")
    parser.add_argument(
        "--num-threads",
        help="number of threads to use per blastn instance (one for each input fastq). Additional threads are also used for gzip and other aspects.",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--stats-file",
        help="filename where to output statistics as a CSV file",
        default=None,
    )

    args = parser.parse_args()

    from ribo_remover.ribo_remover import ribo_remover

    ribo_remover(
        input_fastqs=args.input_fastqs,
        output_fastqs=args.output_fastqs,
        num_threads=args.num_threads,
        stats_file=args.stats_file,
    )
