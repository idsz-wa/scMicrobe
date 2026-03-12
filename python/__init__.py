"""
scMicrobe Python package.
"""

__version__ = "0.1.0"

from .bam_parser import BAMParser, ReadRecord
from .host_filter import HostFilter
from .microbe_aligner import MicrobeAligner, TaxonomyAssignment
from .quantifier import Quantifier
from .contamination_filter import ContaminationFilter
from .output_writer import OutputWriter

__all__ = [
    'BAMParser',
    'ReadRecord',
    'HostFilter',
    'MicrobeAligner',
    'TaxonomyAssignment',
    'Quantifier',
    'ContaminationFilter',
    'OutputWriter',
]
