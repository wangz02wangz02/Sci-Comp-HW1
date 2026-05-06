"""FASTA parsing with IUPAC validation and generator-based streaming.

The parser yields one ``(description, sequence)`` record at a time so that
files containing many sequences do not need to be held in memory all at once.
Validation runs incrementally as each line is read: the set of plausible
alphabets (DNA, RNA, protein) is narrowed character-by-character, so an
illegal residue raises ``FastaFormatError`` on the *line* it appears on rather
than after the entire record has been buffered.

For workloads that only need aggregate statistics (e.g. record counts or
alphabet breakdowns) the ``summarize_fasta`` helper streams the file without
ever retaining a full sequence in memory.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Iterator, Tuple

# IUPAC-compliant alphabets. Gap ('-') is allowed everywhere.
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN-")
RNA_ALPHABET = frozenset("ACGURYSWKMBDHVN-")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBZXU*-")

ALPHABETS: dict[str, frozenset[str]] = {
    "DNA": DNA_ALPHABET,
    "RNA": RNA_ALPHABET,
    "protein": PROTEIN_ALPHABET,
}

FastaRecord = Tuple[str, str]


class FastaFormatError(ValueError):
    """Raised when a FASTA file is malformed or contains invalid characters."""


def _classify_final(possible: set[str]) -> str:
    """Pick the narrowest matching alphabet (DNA > RNA > protein)."""
    for name in ("DNA", "RNA", "protein"):
        if name in possible:
            return name
    return "protein"  # unreachable when called after a successful parse


def _validate_line(
    line_chars: str,
    possible: set[str],
    file_path: str,
    line_number: int,
) -> None:
    """Narrow ``possible`` to alphabets that still admit every char in ``line_chars``.

    Raises ``FastaFormatError`` (with the offending character and line number)
    the moment no alphabet remains.
    """
    for char in line_chars:
        still_valid = {name for name in possible if char in ALPHABETS[name]}
        if not still_valid:
            raise FastaFormatError(
                f"{file_path}:{line_number}: invalid character {char!r} "
                f"(no IUPAC alphabet matches the sequence so far)"
            )
        possible.intersection_update(still_valid)


def parse_fasta(file_path: str) -> Iterator[FastaRecord]:
    """Yield ``(description, sequence)`` tuples from a FASTA file.

    The file is opened once and consumed line-by-line. Only one record is
    materialized at a time (the in-flight description plus the chunks of its
    sequence); previously yielded records are released to the consumer and not
    retained by the parser. Each character is validated as it is read.
    """
    description: str | None = None
    sequence_parts: list[str] = []
    possible: set[str] = set(ALPHABETS)
    saw_header = False

    with open(file_path, "r") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                if description is not None and sequence_parts:
                    raise FastaFormatError(
                        f"{file_path}:{line_number}: blank line inside a record"
                    )
                continue

            if line.startswith(">"):
                saw_header = True
                if description is not None:
                    yield _finalize(description, sequence_parts, possible, file_path)
                description = line[1:].strip()
                if not description:
                    raise FastaFormatError(
                        f"{file_path}:{line_number}: empty FASTA header"
                    )
                sequence_parts = []
                possible = set(ALPHABETS)
            else:
                if description is None:
                    raise FastaFormatError(
                        f"{file_path}:{line_number}: sequence data before any '>' header"
                    )
                upper = line.upper()
                _validate_line(upper, possible, file_path, line_number)
                sequence_parts.append(upper)

    if not saw_header:
        raise FastaFormatError(f"{file_path}: no FASTA header found")

    if description is not None:
        yield _finalize(description, sequence_parts, possible, file_path)


def _finalize(
    description: str,
    sequence_parts: list[str],
    possible: set[str],
    file_path: str,
) -> FastaRecord:
    if not sequence_parts:
        raise FastaFormatError(
            f"{file_path}: record '{description}' has no sequence data"
        )
    sequence = "".join(sequence_parts)
    # ``possible`` is non-empty because _validate_line raises before it could empty.
    _ = _classify_final(possible)
    return description, sequence


def format_record(file_path: str, description: str, sequence: str) -> str:
    """Format one record in the assignment's required output format."""
    return f"- {os.path.basename(file_path)} | {description} | {len(sequence)} | {sequence}"


@dataclass
class FastaSummary:
    """Aggregate statistics computed without holding any full sequence in memory."""

    record_count: int = 0
    total_length: int = 0
    alphabets: dict[str, int] = field(default_factory=dict)
    longest: tuple[str, int] | None = None  # (description, length)


def summarize_fasta(file_path: str) -> FastaSummary:
    """Stream the file and return aggregate stats; never retains a full sequence.

    Useful when a file is too large to fit in memory but you still want a
    record count, total residue count, alphabet breakdown, and the longest
    record's identifier.
    """
    summary = FastaSummary()
    description: str | None = None
    running_length = 0
    possible: set[str] = set(ALPHABETS)
    saw_header = False

    def commit() -> None:
        nonlocal running_length
        if description is None:
            return
        if running_length == 0:
            raise FastaFormatError(
                f"{file_path}: record '{description}' has no sequence data"
            )
        alphabet = _classify_final(possible)
        summary.record_count += 1
        summary.total_length += running_length
        summary.alphabets[alphabet] = summary.alphabets.get(alphabet, 0) + 1
        if summary.longest is None or running_length > summary.longest[1]:
            summary.longest = (description, running_length)
        running_length = 0

    with open(file_path, "r") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                saw_header = True
                commit()
                description = line[1:].strip()
                if not description:
                    raise FastaFormatError(
                        f"{file_path}:{line_number}: empty FASTA header"
                    )
                possible = set(ALPHABETS)
            else:
                if description is None:
                    raise FastaFormatError(
                        f"{file_path}:{line_number}: sequence data before any '>' header"
                    )
                upper = line.upper()
                _validate_line(upper, possible, file_path, line_number)
                running_length += len(upper)

    if not saw_header:
        raise FastaFormatError(f"{file_path}: no FASTA header found")
    commit()
    return summary
