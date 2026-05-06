"""Per-position sensitivity analysis of the standard genetic code."""

from __future__ import annotations

from dataclasses import dataclass


CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

NUCLEOTIDES: tuple[str, ...] = ("A", "C", "G", "T")


@dataclass
class SensitivityResult:
    changes: list[int]
    totals: list[int]

    @property
    def fractions(self) -> list[float]:
        return [c / t if t else 0.0 for c, t in zip(self.changes, self.totals)]

    @property
    def most_sensitive_position(self) -> int:
        fractions = self.fractions
        return fractions.index(max(fractions)) + 1

    @property
    def least_sensitive_position(self) -> int:
        fractions = self.fractions
        return fractions.index(min(fractions)) + 1

    @property
    def sensitivity_ratio(self) -> float:
        fractions = self.fractions
        return max(fractions) / min(fractions)


def analyze_codon_sensitivity() -> SensitivityResult:
    """For each codon position, count how many single-nucleotide subs change the AA.

    Stop codons are excluded as starting points; mutations that introduce a stop
    codon are counted as amino-acid changes.
    """
    sense_codons = [codon for codon, aa in CODON_TABLE.items() if aa != "*"]
    changes = [0, 0, 0]
    totals = [0, 0, 0]

    for codon in sense_codons:
        original_aa = CODON_TABLE[codon]
        for pos in range(3):
            for nt in NUCLEOTIDES:
                if nt == codon[pos]:
                    continue
                mutated = codon[:pos] + nt + codon[pos + 1:]
                totals[pos] += 1
                if CODON_TABLE[mutated] != original_aa:
                    changes[pos] += 1

    return SensitivityResult(changes=changes, totals=totals)


def format_sensitivity_report(result: SensitivityResult) -> str:
    """Render a fixed-width report of the sensitivity analysis."""
    lines = [
        "Codon Position Sensitivity Analysis",
        "=" * 55,
        f"{'Position':<12} {'Changed':<10} {'Total':<10} {'Fraction':<10}",
        "-" * 55,
    ]
    for pos in range(3):
        frac = result.fractions[pos]
        lines.append(
            f"  {pos + 1:<10} {result.changes[pos]:<10} {result.totals[pos]:<10} {frac:<10.4f}"
        )
    lines.append("")
    lines.append(
        f"Most sensitive position:  {result.most_sensitive_position} "
        f"({max(result.fractions):.4f})"
    )
    lines.append(
        f"Least sensitive position: {result.least_sensitive_position} "
        f"({min(result.fractions):.4f})"
    )
    lines.append(f"Ratio (most / least):     {result.sensitivity_ratio:.2f}x")
    return "\n".join(lines)
