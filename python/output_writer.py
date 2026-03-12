"""
Output Writer module for scMicrobe.

Writes output in formats compatible with Seurat and Scanpy:
- matrix.mtx (Matrix Market format)
- barcodes.tsv (cell barcodes)
- features.tsv (microbe features)
- microbe_cell_matrix.tsv (dense TSV format)
"""

import os
import logging
from typing import List, Dict, Optional
from pathlib import Path
import json

import numpy as np
import scipy.sparse as sp
import scipy.io
import pandas as pd
import anndata


class OutputWriter:
    """
    Write scMicrobe output files.
    """
    
    def __init__(
        self,
        output_dir: str,
        prefix: str = 'scmicro'
    ):
        """
        Initialize output writer.
        
        Args:
            output_dir: Output directory path
            prefix: Prefix for output files
        """
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.logger = logging.getLogger('OutputWriter')
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Store feature and cell names
        self.features: List[str] = []
        self.cells: List[str] = []
    
    def write_matrix(
        self,
        matrix: sp.csr_matrix,
        features: Optional[List[str]] = None,
        cells: Optional[List[str]] = None,
        feature_key: str = 'species'
    ):
        """
        Write microbe-cell matrix in multiple formats.
        
        Args:
            matrix: Sparse matrix (microbes x cells)
            features: List of feature names (microbes)
            cells: List of cell barcodes
            feature_key: Type of feature (species, genus, etc.)
        """
        n_microbes, n_cells = matrix.shape
        
        # Generate default names if not provided
        if features is None:
            features = [f"{feature_key}_{i}" for i in range(n_microbes)]
        if cells is None:
            cells = [f"cell_{i}" for i in range(n_cells)]
        
        self.features = features
        self.cells = cells
        
        self.logger.info(f"Writing matrix: {n_microbes} microbes x {n_cells} cells")
        
        # Write MTX format (for Seurat/Scanpy)
        self._write_mtx(matrix)
        
        # Write TSV format (dense, human-readable)
        self._write_tsv(matrix)
        
        # Write H5AD format (for Scanpy)
        self._write_h5ad(matrix, features, cells)
    
    def _write_mtx(self, matrix: sp.csr_matrix):
        """
        Write matrix in Matrix Market format.
        
        Files written:
        - matrix.mtx
        - barcodes.tsv
        - features.tsv
        """
        mtx_dir = self.output_dir / 'mtx'
        mtx_dir.mkdir(exist_ok=True)
        
        # Write matrix
        mtx_path = mtx_dir / 'matrix.mtx'
        self.logger.info(f"Writing MTX: {mtx_path}")
        
        # Try to use scipy.io.mmwrite, fallback to custom implementation
        try:
            scipy.io.mmwrite(str(mtx_path), matrix.T)
        except (AttributeError, ImportError):
            # Fallback: write MTX format manually
            self._write_mtx_manual(matrix, mtx_path)
        
        # Write barcodes
        barcodes_path = mtx_dir / 'barcodes.tsv'
        with open(barcodes_path, 'w') as f:
            for cell in self.cells:
                f.write(f"{cell}\n")
        
        # Write features
        features_path = mtx_dir / 'features.tsv'
        with open(features_path, 'w') as f:
            for feature in self.features:
                f.write(f"{feature}\n")
        
        self.logger.info(f"MTX files written to {mtx_dir}")
    
    def _write_mtx_manual(self, matrix: sp.csr_matrix, mtx_path: Path):
        """
        Manually write Matrix Market format (coordinate format).
        """
        # Transpose to cells x features
        matrix_t = matrix.T.tocoo()
        
        n_rows, n_cols = matrix_t.shape
        n_nnz = matrix_t.nnz
        
        with open(mtx_path, 'w') as f:
            # Write header
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write(f"{n_rows} {n_cols} {n_nnz}\n")
            
            # Write data (1-indexed for MTX format)
            for i, j, v in zip(matrix_t.row, matrix_t.col, matrix_t.data):
                f.write(f"{i+1} {j+1} {int(v)}\n")
    
    def _write_tsv(self, matrix: sp.csr_matrix):
        """
        Write dense matrix in TSV format.
        
        Format:
        - First column: microbe name
        - Remaining columns: cell barcodes
        """
        tsv_path = self.output_dir / f'{self.prefix}_matrix.tsv'
        self.logger.info(f"Writing TSV: {tsv_path}")
        
        # Convert to dense for TSV output
        dense_matrix = matrix.toarray()
        
        # Create DataFrame
        df = pd.DataFrame(
            dense_matrix,
            index=self.features,
            columns=self.cells
        )
        
        # Write TSV
        df.to_csv(tsv_path, sep='\t')
        
        self.logger.info(f"TSV file written: {tsv_path}")
    
    def _write_h5ad(
        self,
        matrix: sp.csr_matrix,
        features: List[str],
        cells: List[str]
    ):
        """
        Write H5AD format for Scanpy.
        
        Args:
            matrix: Sparse matrix (microbes x cells)
            features: List of feature names
            cells: List of cell barcodes
        """
        h5ad_path = self.output_dir / f'{self.prefix}.h5ad'
        self.logger.info(f"Writing H5AD: {h5ad_path}")
        
        # Create AnnData object
        # AnnData expects cells x features, so we transpose
        adata = anndata.AnnData(
            X=matrix.T,
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(index=features)
        )
        
        # Add metadata
        adata.obs['cell_barcode'] = cells
        adata.var['species'] = features
        
        # Write file
        adata.write_h5ad(h5ad_path)
        
        self.logger.info(f"H5AD file written: {h5ad_path}")
    
    def write_stats(self, stats: Dict):
        """
        Write pipeline statistics to JSON.
        
        Args:
            stats: Dictionary of statistics
        """
        stats_path = self.output_dir / f'{self.prefix}_stats.json'
        self.logger.info(f"Writing stats: {stats_path}")
        
        # Convert non-serializable values
        serializable_stats = {}
        for key, value in stats.items():
            if isinstance(value, (np.integer, np.floating)):
                serializable_stats[key] = float(value)
            elif isinstance(value, tuple):
                serializable_stats[key] = list(value)
            else:
                serializable_stats[key] = value
        
        with open(stats_path, 'w') as f:
            json.dump(serializable_stats, f, indent=2)
        
        self.logger.info(f"Stats written: {stats_path}")
    
    def write_qc_report(
        self,
        matrix: sp.csr_matrix,
        output_path: Optional[str] = None
    ):
        """
        Write QC report with summary statistics.
        
        Args:
            matrix: Microbe-cell matrix
            output_path: Path for QC report
        """
        if output_path is None:
            output_path = self.output_dir / f'{self.prefix}_qc.txt'
        
        self.logger.info(f"Writing QC report: {output_path}")
        
        # Calculate QC metrics
        n_microbes, n_cells = matrix.shape
        
        # Microbes per cell
        microbes_per_cell = np.array((matrix > 0).sum(axis=0)).flatten()
        
        # Reads per cell
        reads_per_cell = np.array(matrix.sum(axis=0)).flatten()
        
        # Cells per microbe
        cells_per_microbe = np.array((matrix > 0).sum(axis=1)).flatten()
        
        # Reads per microbe
        reads_per_microbe = np.array(matrix.sum(axis=1)).flatten()
        
        # Write report
        with open(output_path, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("scMicrobe QC Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write("Matrix Summary:\n")
            f.write(f"  Microbes: {n_microbes:,}\n")
            f.write(f"  Cells: {n_cells:,}\n")
            f.write(f"  Total counts: {matrix.sum():,}\n")
            f.write(f"  Sparsity: {1 - matrix.nnz / (n_microbes * n_cells):.4f}\n")
            f.write("\n")
            
            f.write("Microbes per Cell:\n")
            f.write(f"  Mean: {np.mean(microbes_per_cell):.2f}\n")
            f.write(f"  Median: {np.median(microbes_per_cell):.2f}\n")
            f.write(f"  Max: {np.max(microbes_per_cell)}\n")
            f.write(f"  Min: {np.min(microbes_per_cell)}\n")
            f.write("\n")
            
            f.write("Reads per Cell:\n")
            f.write(f"  Mean: {np.mean(reads_per_cell):.2f}\n")
            f.write(f"  Median: {np.median(reads_per_cell):.2f}\n")
            f.write(f"  Max: {np.max(reads_per_cell)}\n")
            f.write(f"  Min: {np.min(reads_per_cell)}\n")
            f.write("\n")
            
            f.write("Cells per Microbe:\n")
            f.write(f"  Mean: {np.mean(cells_per_microbe):.2f}\n")
            f.write(f"  Median: {np.median(cells_per_microbe):.2f}\n")
            f.write(f"  Max: {np.max(cells_per_microbe)}\n")
            f.write(f"  Min: {np.min(cells_per_microbe)}\n")
            f.write("\n")
            
            f.write("Reads per Microbe:\n")
            f.write(f"  Mean: {np.mean(reads_per_microbe):.2f}\n")
            f.write(f"  Median: {np.median(reads_per_microbe):.2f}\n")
            f.write(f"  Max: {np.max(reads_per_microbe)}\n")
            f.write(f"  Min: {np.min(reads_per_microbe)}\n")
            f.write("\n")
        
        self.logger.info(f"QC report written: {output_path}")
    
    def write_cell_metadata(
        self,
        metadata: pd.DataFrame,
        output_path: Optional[str] = None
    ):
        """
        Write cell metadata.
        
        Args:
            metadata: DataFrame with cell metadata
            output_path: Output path
        """
        if output_path is None:
            output_path = self.output_dir / f'{self.prefix}_cell_metadata.tsv'
        
        self.logger.info(f"Writing cell metadata: {output_path}")
        metadata.to_csv(output_path, sep='\t')


def main():
    """Test the output writer."""
    import sys
    import tempfile
    
    logging.basicConfig(level=logging.INFO)
    
    # Create test matrix
    data = np.array([
        [10, 0, 5],
        [0, 8, 0],
        [1, 1, 2],
    ])
    matrix = sp.csr_matrix(data)
    
    features = ['E.coli', 'S.aureus', 'B.subtilis']
    cells = ['cell_A', 'cell_B', 'cell_C']
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"\nTest output directory: {tmpdir}")
        
        writer = OutputWriter(tmpdir, prefix='test')
        
        # Write matrix
        writer.write_matrix(matrix, features, cells)
        
        # Write stats
        writer.write_stats({
            'total_reads': 1000,
            'filtered_reads': 100,
            'final_matrix_shape': matrix.shape
        })
        
        # Write QC report
        writer.write_qc_report(matrix)
        
        # List output files
        print("\nOutput files:")
        for f in Path(tmpdir).iterdir():
            print(f"  {f.name}")
        
        mtx_dir = Path(tmpdir) / 'mtx'
        if mtx_dir.exists():
            print("\nMTX files:")
            for f in mtx_dir.iterdir():
                print(f"  {f.name}")


if __name__ == '__main__':
    main()
