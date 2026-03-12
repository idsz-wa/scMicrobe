"""
Contamination Filter module for scMicrobe.

Implements cell-level contamination filtering:
1. Ambient background estimation
2. Cell specificity scoring
3. Statistical filtering
"""

import logging
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import numpy as np
import scipy.sparse as sp
from scipy import stats


class ContaminationFilter:
    """
    Filter environmental contamination from microbial signals.
    
    Contamination can arise from:
    - Ambient RNA in the sequencing library
    - Reagent contamination
    - Cross-cell contamination
    """
    
    def __init__(
        self,
        ambient_cutoff: float = 0.1,
        min_cells: int = 3,
        min_reads_per_cell: int = 1,
        specificity_threshold: float = 2.0,
        method: str = 'ambient'
    ):
        """
        Initialize contamination filter.
        
        Args:
            ambient_cutoff: Quantile threshold for ambient RNA (0-1)
            min_cells: Minimum number of cells for microbe detection
            min_reads_per_cell: Minimum reads per cell for detection
            specificity_threshold: Minimum specificity score for retention
            method: Filtering method ('ambient', 'binomial', or 'both')
        """
        self.ambient_cutoff = ambient_cutoff
        self.min_cells = min_cells
        self.min_reads_per_cell = min_reads_per_cell
        self.specificity_threshold = specificity_threshold
        self.method = method
        self.logger = logging.getLogger('ContaminationFilter')
        
        # Statistics
        self.stats = {
            'input_microbes': 0,
            'filtered_ambient': 0,
            'filtered_low_prevalence': 0,
            'filtered_low_specificity': 0,
            'retained_microbes': 0
        }
    
    def filter(self, matrix: sp.csr_matrix) -> sp.csr_matrix:
        """
        Apply contamination filtering to microbe-cell matrix.
        
        Args:
            matrix: Sparse matrix (microbes x cells)
            
        Returns:
            Filtered sparse matrix
        """
        self.logger.info(f"Contamination filtering on matrix: {matrix.shape}")
        
        if matrix.shape[0] == 0:
            return matrix
        
        self.stats['input_microbes'] = matrix.shape[0]
        
        # Identify microbes to keep
        keep_mask = np.ones(matrix.shape[0], dtype=bool)
        
        if self.method in ['ambient', 'both']:
            ambient_mask = self._filter_ambient(matrix)
            keep_mask &= ambient_mask
            self.stats['filtered_ambient'] = np.sum(~ambient_mask)
        
        if self.method in ['prevalence', 'both']:
            prevalence_mask = self._filter_by_prevalence(matrix)
            keep_mask &= prevalence_mask
            self.stats['filtered_low_prevalence'] = np.sum(~prevalence_mask)
        
        # Apply specificity filter
        specificity_mask = self._filter_by_specificity(matrix)
        keep_mask &= specificity_mask
        self.stats['filtered_low_specificity'] = np.sum(~specificity_mask)
        
        self.stats['retained_microbes'] = np.sum(keep_mask)
        
        self.logger.info(f"Retained {self.stats['retained_microbes']}/{self.stats['input_microbes']} microbes")
        
        # Return filtered matrix
        return matrix[keep_mask, :]
    
    def _filter_ambient(self, matrix: sp.csr_matrix) -> np.ndarray:
        """
        Filter ambient RNA contamination.
        
        Microbes with high ambient signal (present in many cells at low levels)
        are considered contaminants.
        """
        n_microbes = matrix.shape[0]
        keep = np.ones(n_microbes, dtype=bool)
        
        # Convert to dense for easier manipulation (for small matrices)
        # For large matrices, use sparse operations
        if n_microbes > 10000:
            return self._filter_ambient_sparse(matrix)
        
        dense_matrix = matrix.toarray()
        
        for i in range(n_microbes):
            row = dense_matrix[i, :]
            
            # Skip if no signal
            if np.sum(row) == 0:
                keep[i] = False
                continue
            
            # Calculate ambient score
            # High fraction of cells with low counts = ambient
            nonzero_cells = np.sum(row > 0)
            total_counts = np.sum(row)
            
            if total_counts == 0:
                keep[i] = False
                continue
            
            # Mean count per positive cell
            mean_count = total_counts / nonzero_cells if nonzero_cells > 0 else 0
            
            # Fraction of cells with this microbe
            frac_cells = nonzero_cells / matrix.shape[1]
            
            # Ambient signature: high prevalence, low mean count
            # Compare to cutoff
            if frac_cells > self.ambient_cutoff and mean_count < self.min_reads_per_cell * 2:
                # This looks like ambient contamination
                # But keep if it has some high-count cells (specific signal)
                max_count = np.max(row)
                if max_count < self.min_reads_per_cell * 5:
                    keep[i] = False
        
        return keep
    
    def _filter_ambient_sparse(self, matrix: sp.csr_matrix) -> np.ndarray:
        """Sparse version of ambient filtering for large matrices."""
        n_microbes = matrix.shape[0]
        keep = np.ones(n_microbes, dtype=bool)
        
        for i in range(n_microbes):
            row = matrix[i, :]
            
            # Get non-zero entries
            nonzero = row.nonzero()[1]
            nonzero_cells = len(nonzero)
            total_counts = row.sum()
            
            if total_counts == 0:
                keep[i] = False
                continue
            
            mean_count = total_counts / nonzero_cells if nonzero_cells > 0 else 0
            frac_cells = nonzero_cells / matrix.shape[1]
            
            if frac_cells > self.ambient_cutoff and mean_count < self.min_reads_per_cell * 2:
                max_count = row.max()
                if max_count < self.min_reads_per_cell * 5:
                    keep[i] = False
        
        return keep
    
    def _filter_by_prevalence(self, matrix: sp.csr_matrix) -> np.ndarray:
        """
        Filter by cell prevalence.
        
        Microbes detected in fewer than min_cells are removed.
        """
        # Count cells per microbe
        cells_per_microbe = np.diff(matrix.indptr)
        
        return cells_per_microbe >= self.min_cells
    
    def _filter_by_specificity(self, matrix: sp.csr_matrix) -> np.ndarray:
        """
        Filter by cell specificity.
        
        Microbes should be enriched in specific cells rather than
        uniformly distributed across all cells.
        """
        n_microbes = matrix.shape[0]
        keep = np.ones(n_microbes, dtype=bool)
        
        # Convert to dense for calculation
        dense_matrix = matrix.toarray()
        
        for i in range(n_microbes):
            row = dense_matrix[i, :]
            
            if np.sum(row) == 0:
                keep[i] = False
                continue
            
            # Calculate specificity score
            # Using coefficient of variation (CV)
            mean_count = np.mean(row)
            std_count = np.std(row)
            
            if mean_count == 0:
                keep[i] = False
                continue
            
            cv = std_count / mean_count
            
            # Also calculate Gini coefficient for inequality
            gini = self._gini_coefficient(row)
            
            # Combined specificity score
            specificity = cv * (1 + gini)
            
            if specificity < self.specificity_threshold:
                keep[i] = False
        
        return keep
    
    def _gini_coefficient(self, x: np.ndarray) -> float:
        """
        Calculate Gini coefficient for measuring inequality.
        
        Values range from 0 (perfect equality) to 1 (perfect inequality).
        """
        # Sort values
        sorted_x = np.sort(x)
        n = len(x)
        
        if n == 0 or np.sum(x) == 0:
            return 0
        
        # Calculate Gini
        cumsum = np.cumsum(sorted_x)
        return (n + 1 - 2 * np.sum(cumsum) / cumsum[-1]) / n
    
    def estimate_ambient_profile(
        self,
        matrix: sp.csr_matrix,
        empty_droplets: Optional[List[int]] = None
    ) -> Dict[int, float]:
        """
        Estimate ambient RNA profile from empty droplets or low-count cells.
        
        Args:
            matrix: Microbe-cell matrix
            empty_droplets: Indices of empty droplets (if known)
            
        Returns:
            Dictionary mapping microbe index to ambient fraction
        """
        if empty_droplets is None:
            # Identify low-count cells as proxy for empty droplets
            total_counts_per_cell = np.array(matrix.sum(axis=0)).flatten()
            threshold = np.percentile(total_counts_per_cell, 10)
            empty_droplets = np.where(total_counts_per_cell <= threshold)[0].tolist()
        
        self.logger.info(f"Estimating ambient from {len(empty_droplets)} cells")
        
        # Calculate ambient profile
        ambient_counts = np.array(matrix[:, empty_droplets].sum(axis=1)).flatten()
        ambient_total = np.sum(ambient_counts)
        
        if ambient_total == 0:
            return {}
        
        ambient_profile = {}
        for i, count in enumerate(ambient_counts):
            ambient_profile[i] = count / ambient_total
        
        return ambient_profile
    
    def correct_ambient(
        self,
        matrix: sp.csr_matrix,
        ambient_profile: Dict[int, float]
    ) -> sp.csr_matrix:
        """
        Correct for ambient RNA contamination.
        
        Subtract expected ambient counts from observed counts.
        """
        # This is a simplified correction
        # Full implementation would use more sophisticated methods
        # like those in CellBender or SoupX
        
        corrected = matrix.copy().astype(np.float32)
        
        for cell_idx in range(matrix.shape[1]):
            cell_total = matrix[:, cell_idx].sum()
            
            for microbe_idx, ambient_frac in ambient_profile.items():
                expected_ambient = cell_total * ambient_frac
                observed = matrix[microbe_idx, cell_idx]
                
                # Subtract expected ambient
                corrected_val = max(0, observed - expected_ambient)
                corrected[microbe_idx, cell_idx] = corrected_val
        
        return corrected
    
    def get_stats(self) -> Dict:
        """Get filtering statistics."""
        return self.stats.copy()


class DoubletDetector:
    """
    Detect cells with multiple microbial signatures (potential doublets).
    """
    
    def __init__(self, max_microbes_per_cell: int = 3):
        """
        Initialize doublet detector.
        
        Args:
            max_microbes_per_cell: Maximum expected microbes per cell
        """
        self.max_microbes_per_cell = max_microbes_per_cell
        self.logger = logging.getLogger('DoubletDetector')
    
    def detect(self, matrix: sp.csr_matrix) -> List[int]:
        """
        Detect potential doublets.
        
        Returns:
            List of cell indices that are potential doublets
        """
        doublets = []
        
        # Count microbes per cell
        microbes_per_cell = np.array((matrix > 0).sum(axis=0)).flatten()
        
        # Cells with too many microbes are potential doublets
        for cell_idx, n_microbes in enumerate(microbes_per_cell):
            if n_microbes > self.max_microbes_per_cell:
                doublets.append(cell_idx)
        
        self.logger.info(f"Detected {len(doublets)} potential doublets")
        
        return doublets


def main():
    """Test the contamination filter."""
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    # Create test matrix
    # 5 microbes x 10 cells
    data = np.array([
        # cell1-10
        [10, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Specific to cell1
        [0, 8, 0, 0, 0, 0, 0, 0, 0, 0],   # Specific to cell2
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],   # Ambient (all cells, low count)
        [0, 0, 0, 0, 0, 0, 0, 0, 5, 0],   # Specific to cell9
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # No signal
    ])
    
    matrix = sp.csr_matrix(data)
    
    print("Test matrix:")
    print(data)
    print(f"Shape: {matrix.shape}")
    
    # Test filtering
    filter_obj = ContaminationFilter(
        ambient_cutoff=0.5,
        min_cells=2,
        specificity_threshold=1.0
    )
    
    filtered = filter_obj.filter(matrix)
    
    print(f"\nFiltered matrix shape: {filtered.shape}")
    print(f"Filtered matrix:\n{filtered.toarray()}")
    print(f"Stats: {filter_obj.get_stats()}")


if __name__ == '__main__':
    main()
