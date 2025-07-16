"""
Converts a fastq file read in from stdin to a fasta file on stdout
but also replaces read IDs with just 1,2,3...
This removes difficulties with handling reads (particularly non-unique read ids)
"""


def main():
    import sys

    next_id = 1
    while True:
        id = sys.stdin.readline()
        seq = sys.stdin.readline()
        _spacer = sys.stdin.readline()  # skip for fasta
        _qual = sys.stdin.readline()  # skip for fasta
        if not id:
            # nothing more to read
            break
        assert id.startswith("@")
        sys.stdout.write(f">{next_id}\n")
        sys.stdout.write(seq)
        next_id += 1
