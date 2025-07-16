import pathlib

from ribo_remover.ribo_remover import ribo_remover


def test_one_fastq(tmp_path):
    out_R1 = tmp_path / "out.R1.fastq"
    stats_file = tmp_path / "stats.txt"
    ribo_remover(
        ["tests/resources/R1.fastq"],
        [str(out_R1)],
        stats_file=stats_file,
    )
    out_contents = out_R1.read_text()
    ref_contents = pathlib.Path("tests/resources/one_fastq.out.R1.fastq").read_text()
    assert out_contents == ref_contents

    out_stat = stats_file.read_text()
    ref_stat = pathlib.Path("tests/resources/one_fastq.stats.txt").read_text()
    assert ref_stat == out_stat


def test_two_fastqs(tmp_path):
    out_R1 = tmp_path / "out.R1.fastq"
    out_R2 = tmp_path / "out.R2.fastq"
    stats_file = tmp_path / "stats.txt"
    ribo_remover(
        ["tests/resources/R1.fastq", "tests/resources/R2.fastq"],
        [str(out_R1), str(out_R2)],
        stats_file=stats_file,
    )
    out_contents = out_R1.read_text()
    ref_contents = pathlib.Path("tests/resources/two_fastqs.out.R1.fastq").read_text()
    assert out_contents == ref_contents

    out_contents = out_R2.read_text()
    ref_contents = pathlib.Path("tests/resources/two_fastqs.out.R2.fastq").read_text()
    assert out_contents == ref_contents

    out_stat = stats_file.read_text()
    ref_stat = pathlib.Path("tests/resources/two_fastqs.stats.txt").read_text()
    assert ref_stat == out_stat
