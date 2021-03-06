"""
phipkit - analysis package for Phage Immunoprecipitation Sequencing (phip-seq)
"""

__version__ = "0.0.1"
from .antigen_analysis import AntigenAnalysis
from .call_antigens import call_antigens
from .call_hits import call_hits
from .score import compute_scores

__all__ = [
    "AntigenAnalysis",
    "call_antigens",
    "compute_scores",
]
