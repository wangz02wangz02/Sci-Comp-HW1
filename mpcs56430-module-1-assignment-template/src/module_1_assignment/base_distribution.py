"""Base composition analysis and plotting for nucleotide sequences."""

from __future__ import annotations

from collections import Counter
from typing import Iterable, Sequence


DEFAULT_BASES: tuple[str, ...] = ("A", "G", "C", "T")
DEFAULT_COLORS: tuple[str, ...] = ("#1f77b4", "#2ca02c", "#d62728", "#ff7f0e")


def read_sequence(file_path: str) -> str:
    """Read a single FASTA file and return the concatenated, uppercased sequence.

    Header lines (starting with '>') and blank lines are ignored. Suitable for
    files containing one record (e.g. a whole-genome FASTA).
    """
    with open(file_path, "r") as handle:
        return "".join(
            line.strip().upper()
            for line in handle
            if line.strip() and not line.startswith(">")
        )


def count_bases(sequence: str, bases: Sequence[str] = DEFAULT_BASES) -> list[int]:
    """Return the count of each base in ``bases`` within ``sequence``."""
    counts = Counter(sequence)
    return [counts.get(base, 0) for base in bases]


def plot_base_distribution(
    sequence: str,
    bases: Sequence[str] = DEFAULT_BASES,
    colors: Iterable[str] = DEFAULT_COLORS,
    title: str = "Base Distribution",
) -> None:
    """Render a labeled bar chart of base counts using the inline backend."""
    import matplotlib.pyplot as plt

    values = count_bases(sequence, bases)
    total = sum(values)
    percents = [v / total * 100 if total else 0.0 for v in values]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(list(bases), values, color=list(colors))
    ax.set_title(title)
    ax.set_xlabel("Base")
    ax.set_ylabel("Count")

    headroom = max(values) * 0.01 if values else 0
    for bar, pct in zip(bars, percents):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + headroom,
            f"{height:,} ({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    fig.tight_layout()
    plt.show()
