"""
Quantifier module for scMicrobe.

Handles UMI deduplication and quantification of microbial reads per cell.
"""

import logging
from typing import List, Dict, Tuple, Set
from collections import defaultdict
from dataclasses import dataclass
import numpy as np
import scipy.sparse as sp

from bam_parser import ReadRecord
from microbe_aligner import TaxonomyAssignment


@dataclass
class UMIKey:
    """Unique key for UMI deduplication."""
    barcode: str
    species: str
    umi: str
    
    def __hash__(self):
        return hash((self.barcode, self.species, self.umi))
    
    def __eq__(self, other):
        return (self.barcode == other.barcode and 
                self.species == other.species and 
                self.umi == other.umi)


class Quantifier:
    """
    Quantify microbial abundance per cell.
    
    Supports two methods:
    1. Read counting: Simple read count per cell-species
    2. UMI deduplicated counting: Count unique UMIs per cell-species
    """
    
    def __init__(
        self,
        method: str = 'umi',
        umi_dedup_method: str = 'exact'
    ):
        """
        Initialize quantifier.
        
        Args:
            method: Quantification method ('reads' or 'umi')
            umi_dedup_method: UMI deduplication method ('exact' or 'hamming')
        """
        self.method = method
        self.umi_dedup_method = umi_dedup_method
        self.logger = logging.getLogger('Quantifier')
        
        # Statistics
        self.stats = {
            'total_reads': 0,
            'valid_reads': 0,
            'reads_without_barcode': 0,
            'reads_without_umi': 0,
            'unique_cell_species_pairs': 0,
            'total_umis': 0,
            'deduplicated_umis': 0
        }
    
    def quantify(
        self,
        aligned_reads: List[Tuple[ReadRecord, TaxonomyAssignment]]
    ) -> sp.csr_matrix:
        """
        Quantify microbial abundance per cell.
        
        Args:
            aligned_reads: List of (ReadRecord, TaxonomyAssignment) tuples
            
        Returns:
            Sparse matrix (microbes x cells)
        """
        self.logger.info(f"Quantifying {len(aligned_reads)} aligned reads")
        
        if not aligned_reads:
            self.logger.warning("No aligned reads to quantify")
            return sp.csr_matrix((0, 0))
        
        # Filter and collect valid reads
        valid_records = self._filter_valid_reads(aligned_reads)
        
        if self.method == 'umi':
            count_dict = self._quantify_by_umi(valid_records)
        else:
            count_dict = self._quantify_by_reads(valid_records)
        
        # Build sparse matrix
        matrix = self._build_matrix(count_dict)
        
        self.logger.info(f"Quantification complete: {matrix.shape}")
        
        return matrix
    
    def _filter_valid_reads(
        self,
        aligned_reads: List[Tuple[ReadRecord, TaxonomyAssignment]]
    ) -> List[Tuple[ReadRecord, TaxonomyAssignment]]:
        """Filter reads with valid barcodes and UMIs."""
        valid = []
        
        for read, assignment in aligned_reads:
            self.stats['total_reads'] += 1
            
            # Skip unclassified reads
            if not assignment.classified:
                continue
            
            # Check for barcode
            if not read.barcode:
                self.stats['reads_without_barcode'] += 1
                continue
            
            # Check for UMI if using UMI method
            if self.method == 'umi' and not read.umi:
                self.stats['reads_without_umi'] += 1
                continue
            
            valid.append((read, assignment))
            self.stats['valid_reads'] += 1
        
        self.logger.info(f"Valid reads: {self.stats['valid_reads']}/{self.stats['total_reads']}")
        
        return valid
    
    def _quantify_by_umi(
        self,
        valid_records: List[Tuple[ReadRecord, TaxonomyAssignment]]
    ) -> Dict[Tuple[str, str], int]:
        """
        Quantify by unique UMIs per cell-species.
        
        Unique key: (cell barcode, species, UMI)
        """
        umi_set: Set[UMIKey] = set()
        
        for read, assignment in valid_records:
            key = UMIKey(
                barcode=read.barcode,
                species=assignment.species,
                umi=read.umi
            )
            umi_set.add(key)
        
        self.stats['total_umis'] = len(valid_records)
        self.stats['deduplicated_umis'] = len(umi_set)
        
        # Count UMIs per cell-species
        count_dict = defaultdict(int)
        for key in umi_set:
            count_dict[(key.species, key.barcode)] += 1
        
        self.stats['unique_cell_species_pairs'] = len(count_dict)
        
        self.logger.info(f"UMI deduplication: {self.stats['total_umis']} -> {self.stats['deduplicated_umis']}")
        
        return dict(count_dict)
    
    def _quantify_by_reads(
        self,
        valid_records: List[Tuple[ReadRecord, TaxonomyAssignment]]
    ) -> Dict[Tuple[str, str], int]:
        """
        Quantify by read counts per cell-species.
        """
        count_dict = defaultdict(int)
        
        for read, assignment in valid_records:
            key = (assignment.species, read.barcode)
            count_dict[key] += 1
        
        self.stats['unique_cell_species_pairs'] = len(count_dict)
        
        return dict(count_dict)
    
    def _build_matrix(
        self,
        count_dict: Dict[Tuple[str, str], int]
    ) -> sp.csr_matrix:
        """
        Build sparse matrix from count dictionary.
        
        Returns:
            Sparse matrix with shape (n_microbes, n_cells)
        """
        # Get unique microbes and cells
        microbes = sorted(set(k[0] for k in count_dict.keys()))
        cells = sorted(set(k[1] for k in count_dict.keys()))
        
        # Create mappings
        microbe_to_idx = {m: i for i, m in enumerate(microbes)}
        cell_to_idx = {c: i for i, c in enumerate(cells)}
        
        # Build matrix data
        rows = []
        cols = []
        data = []
        
        for (species, barcode), count in count_dict.items():
            rows.append(microbe_to_idx[species])
            cols.append(cell_to_idx[barcode])
            data.append(count)
        
        # Create sparse matrix
        matrix = sp.csr_matrix(
            (data, (rows, cols)),
            shape=(len(microbes), len(cells)),
            dtype=np.int32
        )
        
        # Store mappings for later use
        self.microbes = microbes
        self.cells = cells
        
        return matrix
    
    def get_feature_names(self) -> List[str]:
        """Get microbe feature names (species)."""
        return getattr(self, 'microbes', [])
    
    def get_cell_names(self) -> List[str]:
        """Get cell barcode names."""
        return getattr(self, 'cells', [])
    
    def get_stats(self) -> Dict:
        """Get quantification statistics."""
        return self.stats.copy()


class UMIDeduplicator:
    """
    Advanced UMI deduplication with error correction.
    """
    
    def __init__(self, method: str = 'exact', hamming_threshold: int = 1):
        """
        Initialize UMI deduplicator.
        
        Args:
            method: Deduplication method ('exact' or 'hamming')
            hamming_threshold: Hamming distance threshold for clustering
        """
        self.method = method
        self.hamming_threshold = hamming_threshold
        self.logger = logging.getLogger('UMIDeduplicator')
    
    def deduplicate(self, umis: List[str]) -> List[str]:
        """
        Deduplicate UMIs.
        
        Args:
            umis: List of UMI sequences
            
        Returns:
            List of unique UMIs
        """
        if self.method == 'exact':
            return list(set(umis))
        elif self.method == 'hamming':
            return self._deduplicate_hamming(umis)
        else:
            raise ValueError(f"Unknown deduplication method: {self.method}")
    
    def _deduplicate_hamming(self, umis: List[str]) -> List[str]:
        """
        Deduplicate UMIs using Hamming distance clustering.
        
        UMIs within threshold distance are considered duplicates.
        """
        if not umis:
            return []
        
        # Start with unique UMIs
        unique_umis = []
        
        for umi in umis:
            is_duplicate = False
            
            for unique_umi in unique_umis:
                if self._hamming_distance(umi, unique_umi) <= self.hamming_threshold:
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                unique_umis.append(umi)
        
        return unique_umis
    
    def _hamming_distance(self, s1: str, s2: str) -> int:
        """Calculate Hamming distance between two strings."""
        if len(s1) != len(s2):
            # Pad shorter string
            max_len = max(len(s1), len(s2))
            s1 = s1.ljust(max_len, 'N')
            s2 = s2.ljust(max_len, 'N')
        
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def main():
    """Test the quantifier."""
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    # Create test data
    test_reads = [
        (ReadRecord('r1', 'ATCG', 'IIII', 'cell1', 'UMI1', 0, True, False, False), 
         TaxonomyAssignment('r1', 'tx1', 'E.coli', 'Escherichia', 0.9, 100)),
        (ReadRecord('r2', 'ATCG', 'IIII', 'cell1', 'UMI1', 0, True, False, False),  # Duplicate UMI
         TaxonomyAssignment('r2', 'tx1', 'E.coli', 'Escherichia', 0.9, 100)),
        (ReadRecord('r3', 'GGGG', 'IIII', 'cell1', 'UMI2', 0, True, False, False),
         TaxonomyAssignment('r3', 'tx1', 'E.coli', 'Escherichia', 0.9, 100)),
        (ReadRecord('r4', 'TTTT', 'IIII', 'cell2', 'UMI3', 0, True, False, False),
         TaxonomyAssignment('r4', 'tx2', 'S.aureus', 'Staphylococcus', 0.9, 100)),
        (ReadRecord('r5', 'CCCC', 'IIII', 'cell2', 'UMI4', 0, True, False, False),
         TaxonomyAssignment('r5', 'tx2', 'S.aureus', 'Staphylococcus', 0.9, 100)),
    ]
    
    # Test UMI quantification
    print("\nTesting UMI quantification:")
    quantifier_umi = Quantifier(method='umi')
    matrix_umi = quantifier_umi.quantify(test_reads)
    print(f"  Matrix shape: {matrix_umi.shape}")
    print(f"  Features: {quantifier_umi.get_feature_names()}")
    print(f"  Cells: {quantifier_umi.get_cell_names()}")
    print(f"  Stats: {quantifier_umi.get_stats()}")
    print(f"  Matrix:\n{matrix_umi.toarray()}")
    
    # Test read quantification
    print("\nTesting read quantification:")
    quantifier_reads = Quantifier(method='reads')
    matrix_reads = quantifier_reads.quantify(test_reads)
    print(f"  Matrix shape: {matrix_reads.shape}")
    print(f"  Stats: {quantifier_reads.get_stats()}")
    print(f"  Matrix:\n{matrix_reads.toarray()}")


if __name__ == '__main__':
    main()
