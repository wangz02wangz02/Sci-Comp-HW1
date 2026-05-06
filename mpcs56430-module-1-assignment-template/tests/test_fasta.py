import inspect
import textwrap
from pathlib import Path

import pytest

from module_1_assignment.fasta import (
    FastaFormatError,
    format_record,
    parse_fasta,
    summarize_fasta,
)


def write(tmp_path: Path, text: str) -> str:
    path = tmp_path / "input.fasta"
    path.write_text(textwrap.dedent(text).lstrip())
    return str(path)


def test_parse_fasta_is_a_generator() -> None:
    assert inspect.isgeneratorfunction(parse_fasta)


def test_parse_fasta_streams_multiple_records(tmp_path: Path) -> None:
    path = write(
        tmp_path,
        """
        >dna_one
        ACGT
        >dna_two
        AAAA
        CCCC
        """,
    )
    records = list(parse_fasta(path))
    assert records == [("dna_one", "ACGT"), ("dna_two", "AAAACCCC")]


def test_parse_fasta_uppercases_input(tmp_path: Path) -> None:
    path = write(tmp_path, ">low\nacgt\n")
    assert list(parse_fasta(path)) == [("low", "ACGT")]


def test_parse_fasta_rejects_invalid_character(tmp_path: Path) -> None:
    path = write(tmp_path, ">bad\nACGT123\n")
    with pytest.raises(FastaFormatError, match="invalid character"):
        list(parse_fasta(path))


def test_parse_fasta_rejects_blank_line_inside_record(tmp_path: Path) -> None:
    path = write(tmp_path, ">x\nACGT\n\nACGT\n")
    with pytest.raises(FastaFormatError, match="blank line"):
        list(parse_fasta(path))


def test_parse_fasta_rejects_data_before_header(tmp_path: Path) -> None:
    path = write(tmp_path, "ACGT\n>x\nACGT\n")
    with pytest.raises(FastaFormatError, match="before any '>' header"):
        list(parse_fasta(path))


def test_parse_fasta_rejects_empty_header(tmp_path: Path) -> None:
    path = write(tmp_path, ">\nACGT\n")
    with pytest.raises(FastaFormatError, match="empty FASTA header"):
        list(parse_fasta(path))


def test_parse_fasta_accepts_protein(tmp_path: Path) -> None:
    path = write(tmp_path, ">p\nMKVLWAALLVTFLAGCQA*\n")
    [(_, sequence)] = list(parse_fasta(path))
    assert sequence.startswith("MKVL")


def test_format_record_uses_basename() -> None:
    assert format_record("/tmp/example.fa", "id", "ACGT") == "- example.fa | id | 4 | ACGT"


def test_summarize_fasta_streams_without_buffering(tmp_path: Path) -> None:
    path = write(
        tmp_path,
        """
        >a
        ACGT
        >b
        AAAACCCC
        >c
        MKLW*
        """,
    )
    summary = summarize_fasta(path)
    assert summary.record_count == 3
    assert summary.total_length == 4 + 8 + 5
    assert summary.alphabets == {"DNA": 2, "protein": 1}
    assert summary.longest == ("b", 8)
