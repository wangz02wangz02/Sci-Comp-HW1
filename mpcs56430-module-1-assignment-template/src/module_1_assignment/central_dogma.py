"""Central dogma simulation: transcription, reverse complement, translation, and ORF finding."""

from __future__ import annotations

from typing import List, Tuple


RNA_CODON_TABLE: dict[str, str] = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})

OrfRecord = Tuple[int, int, int, int, str]


class CentralDogma:
    """Model the DNA → RNA → protein pipeline for a single sequence record."""

    def __init__(self, description: str, sequence: str) -> None:
        self.description = description
        self.sequence = sequence.upper()

    def transcribe(self, dna_sequence: str | None = None) -> str:
        dna = (dna_sequence or self.sequence).upper()
        return dna.replace("T", "U")

    def reverse_complement(self, dna_sequence: str | None = None) -> str:
        dna = (dna_sequence or self.sequence).upper()
        complement = str.maketrans("ACGT", "TGCA")
        return dna.translate(complement)[::-1]

    def translate(self, rna_sequence: str) -> str:
        protein: list[str] = []
        for index in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[index:index + 3]
            amino_acid = RNA_CODON_TABLE.get(codon)
            if amino_acid is None or amino_acid == "*":
                break
            protein.append(amino_acid)
        return "".join(protein)

    def find_orfs(self) -> List[OrfRecord]:
        """Return all ORFs across the six reading frames as (frame, start, end, length, protein)."""
        orfs: list[OrfRecord] = []
        sequence_length = len(self.sequence)

        def scan_strand(dna_sequence: str, frame_label: int, reverse: bool = False) -> None:
            offset = abs(frame_label) - 1
            for start_index in range(offset, len(dna_sequence) - 2, 3):
                if dna_sequence[start_index:start_index + 3] != "ATG":
                    continue

                for stop_index in range(start_index + 3, len(dna_sequence) - 2, 3):
                    stop_codon = dna_sequence[stop_index:stop_index + 3]
                    if stop_codon not in STOP_CODONS:
                        continue

                    coding_dna = dna_sequence[start_index:stop_index]
                    protein = self.translate(self.transcribe(coding_dna))
                    if protein:
                        if reverse:
                            start_base = sequence_length - start_index
                            end_base = sequence_length - (stop_index + 2)
                        else:
                            start_base = start_index + 1
                            end_base = stop_index + 3

                        orfs.append((frame_label, start_base, end_base, len(protein), protein))
                    break

        for frame_label in (1, 2, 3):
            scan_strand(self.sequence, frame_label, reverse=False)

        reverse_sequence = self.reverse_complement()
        for frame_label in (-1, -2, -3):
            scan_strand(reverse_sequence, frame_label, reverse=True)

        return orfs


def format_orfs(description: str, orfs: List[OrfRecord]) -> str:
    """Format ORFs in the assignment's required output layout."""
    lines = [f">{description}", "frame | start | stop | length | sequence"]
    for frame, start, stop, length, protein in orfs:
        lines.append(f"* {frame} | {start} | {stop} | {length} | {protein}")
    return "\n".join(lines)
