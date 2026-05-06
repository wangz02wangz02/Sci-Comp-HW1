"""Top-level package for the scientific computing project."""

from .base_distribution import (
    DEFAULT_BASES,
    count_bases,
    plot_base_distribution,
    read_sequence,
)
from .central_dogma import CentralDogma, format_orfs
from .codon_sensitivity import (
    SensitivityResult,
    analyze_codon_sensitivity,
    format_sensitivity_report,
)
from .fasta import FastaFormatError, FastaSummary, format_record, parse_fasta, summarize_fasta
from .monte_carlo_pi import plot_pi_simulation, simulate_pi
from .physics import compute_kinetic_energy

__all__ = [
    "CentralDogma",
    "DEFAULT_BASES",
    "FastaFormatError",
    "FastaSummary",
    "SensitivityResult",
    "analyze_codon_sensitivity",
    "compute_kinetic_energy",
    "count_bases",
    "format_orfs",
    "format_record",
    "format_sensitivity_report",
    "parse_fasta",
    "plot_base_distribution",
    "plot_pi_simulation",
    "read_sequence",
    "simulate_pi",
    "summarize_fasta",
]
