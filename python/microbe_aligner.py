"""
Microbe Aligner module for scMicrobe.

Supports two alignment strategies:
1. Kraken2: Fast taxonomic classification using k-mer matching
2. minimap2: Accurate alignment allowing sensitive microbial detection
"""

import os
import logging
import tempfile
import subprocess
import json
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from collections import defaultdict

from bam_parser import ReadRecord


@dataclass
class TaxonomyAssignment:
    """Taxonomic assignment result for a read."""
    read_id: str
    taxid: str
    species: str
    genus: str
    confidence: float
    alignment_length: int
    num_kmers: Optional[int] = None
    classified: bool = True


class MicrobeAligner:
    """
    Align reads to microbial reference and assign taxonomy.
    """
    
    def __init__(
        self,
        aligner: str = 'minimap2',
        microbe_ref: Optional[str] = None,
        kraken_db: Optional[str] = None,
        minimap2_index: Optional[str] = None,
        min_mapq: int = 1,
        threads: int = 4,
        confidence_threshold: float = 0.0
    ):
        """
        Initialize microbe aligner.
        
        Args:
            aligner: Alignment method ('minimap2' or 'kraken2')
            microbe_ref: Path to microbial reference FASTA
            kraken_db: Path to Kraken2 database
            minimap2_index: Path to pre-built minimap2 index
            min_mapq: Minimum MAPQ for valid alignment
            threads: Number of threads
            confidence_threshold: Minimum confidence for taxonomy assignment
        """
        self.aligner = aligner
        self.microbe_ref = microbe_ref
        self.kraken_db = kraken_db
        self.minimap2_index = minimap2_index or microbe_ref
        self.min_mapq = min_mapq
        self.threads = threads
        self.confidence_threshold = confidence_threshold
        self.logger = logging.getLogger('MicrobeAligner')
        
        # Taxonomy mapping (would be loaded from database)
        self.taxonomy_map = {}
        
        # Statistics
        self.stats = {
            'total_reads': 0,
            'classified_reads': 0,
            'unclassified_reads': 0,
            'unique_species': set()
        }
        
        # Validate inputs
        self._validate_inputs()
    
    def _validate_inputs(self):
        """Validate input parameters."""
        if self.aligner == 'kraken2':
            if not self.kraken_db:
                raise ValueError("Kraken2 alignment requires --kraken-db")
            if not os.path.exists(self.kraken_db):
                raise FileNotFoundError(f"Kraken2 database not found: {self.kraken_db}")
        elif self.aligner == 'minimap2':
            if not self.minimap2_index:
                raise ValueError("Minimap2 alignment requires --microbe-ref or --minimap2-index")
    
    def align(self, reads: List[ReadRecord]) -> List[Tuple[ReadRecord, TaxonomyAssignment]]:
        """
        Align reads to microbial reference and assign taxonomy.
        
        Args:
            reads: List of ReadRecord objects
            
        Returns:
            List of (ReadRecord, TaxonomyAssignment) tuples
        """
        if not reads:
            return []
        
        self.logger.info(f"Aligning {len(reads)} reads using {self.aligner}")
        
        if self.aligner == 'kraken2':
            return self._align_with_kraken2(reads)
        elif self.aligner == 'minimap2':
            return self._align_with_minimap2(reads)
        else:
            raise ValueError(f"Unknown aligner: {self.aligner}")
    
    def _align_with_kraken2(
        self,
        reads: List[ReadRecord]
    ) -> List[Tuple[ReadRecord, TaxonomyAssignment]]:
        """
        Classify reads using Kraken2.
        """
        # Write reads to temporary FASTQ
        with tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.fastq',
            delete=False
        ) as tmp_fastq:
            tmp_fastq_path = tmp_fastq.name
            for read in reads:
                tmp_fastq.write(read.to_fastq())
        
        results = []
        
        try:
            # Run Kraken2
            output_path = tmp_fastq_path + '.kraken2'
            report_path = tmp_fastq_path + '.kreport2'
            
            cmd = [
                'kraken2',
                '--db', self.kraken_db,
                '--threads', str(self.threads),
                '--output', output_path,
                '--report', report_path,
                '--use-names',
                tmp_fastq_path
            ]
            
            self.logger.info(f"Running Kraken2: {' '.join(cmd)}")
            
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )
            
            if process.returncode != 0:
                self.logger.error(f"Kraken2 failed: {process.stderr}")
                raise RuntimeError(f"Kraken2 alignment failed: {process.stderr}")
            
            # Parse Kraken2 output
            assignments = self._parse_kraken2_output(output_path)
            
            # Match assignments with reads
            read_dict = {r.read_id: r for r in reads}
            
            for read_id, assignment in assignments.items():
                if read_id in read_dict:
                    read = read_dict[read_id]
                    results.append((read, assignment))
                    
                    # Update stats
                    self.stats['total_reads'] += 1
                    if assignment.classified:
                        self.stats['classified_reads'] += 1
                        self.stats['unique_species'].add(assignment.species)
                    else:
                        self.stats['unclassified_reads'] += 1
            
            self.logger.info(f"Kraken2 classified {self.stats['classified_reads']} reads")
            
        finally:
            # Cleanup
            for path in [tmp_fastq_path, output_path, report_path]:
                if os.path.exists(path):
                    os.unlink(path)
        
        return results
    
    def _parse_kraken2_output(self, output_path: str) -> Dict[str, TaxonomyAssignment]:
        """Parse Kraken2 output file."""
        assignments = {}
        
        with open(output_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                
                classified = parts[0] == 'C'
                read_id = parts[1]
                taxid = parts[2]
                seq_len = int(parts[3])
                taxonomy = parts[4] if len(parts) > 4 else 'unclassified'
                
                # Parse taxonomy string (e.g., "root;cellular organisms;Bacteria;...")
                species = 'unclassified'
                genus = 'unclassified'
                
                if taxonomy != 'unclassified':
                    tax_parts = taxonomy.split(';')
                    if tax_parts:
                        species = tax_parts[-1]
                    if len(tax_parts) > 1:
                        genus = tax_parts[-2]
                
                assignment = TaxonomyAssignment(
                    read_id=read_id,
                    taxid=taxid,
                    species=species,
                    genus=genus,
                    confidence=1.0 if classified else 0.0,
                    alignment_length=seq_len,
                    classified=classified
                )
                
                assignments[read_id] = assignment
        
        return assignments
    
    def _align_with_minimap2(
        self,
        reads: List[ReadRecord]
    ) -> List[Tuple[ReadRecord, TaxonomyAssignment]]:
        """
        Align reads using minimap2 and assign taxonomy based on reference.
        """
        # Write reads to temporary FASTQ
        with tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.fastq',
            delete=False
        ) as tmp_fastq:
            tmp_fastq_path = tmp_fastq.name
            for read in reads:
                tmp_fastq.write(read.to_fastq())
        
        results = []
        
        try:
            # Run minimap2
            sam_path = tmp_fastq_path + '.sam'
            
            cmd = [
                'minimap2',
                '-ax', 'sr',
                '-t', str(self.threads),
                '--secondary=no',
                self.minimap2_index,
                tmp_fastq_path
            ]
            
            self.logger.info(f"Running minimap2: {' '.join(cmd)}")
            
            with open(sam_path, 'w') as sam_file:
                process = subprocess.run(
                    cmd,
                    stdout=sam_file,
                    stderr=subprocess.PIPE,
                    text=True
                )
            
            if process.returncode != 0:
                self.logger.error(f"minimap2 failed: {process.stderr}")
                raise RuntimeError(f"minimap2 alignment failed: {process.stderr}")
            
            # Parse SAM output and assign taxonomy
            assignments = self._parse_minimap2_sam(sam_path)
            
            # Match assignments with reads
            read_dict = {r.read_id: r for r in reads}
            
            for read_id, assignment in assignments.items():
                if read_id in read_dict:
                    read = read_dict[read_id]
                    results.append((read, assignment))
                    
                    # Update stats
                    self.stats['total_reads'] += 1
                    if assignment.classified:
                        self.stats['classified_reads'] += 1
                        self.stats['unique_species'].add(assignment.species)
                    else:
                        self.stats['unclassified_reads'] += 1
            
            self.logger.info(f"minimap2 aligned {self.stats['classified_reads']} reads")
            
        finally:
            # Cleanup
            for path in [tmp_fastq_path, sam_path]:
                if os.path.exists(path):
                    os.unlink(path)
        
        return results
    
    def _parse_minimap2_sam(self, sam_path: str) -> Dict[str, TaxonomyAssignment]:
        """Parse minimap2 SAM output and extract taxonomy."""
        assignments = {}
        
        with open(sam_path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                
                read_id = parts[0]
                flag = int(parts[1])
                ref_name = parts[2]
                mapq = int(parts[4])
                
                # Skip unmapped reads
                if flag & 0x4:
                    continue
                
                # Skip low MAPQ alignments
                if mapq < self.min_mapq:
                    continue
                
                # Extract species from reference name
                # Reference name format could be: species|strain|accession
                species = self._extract_species_from_ref(ref_name)
                
                assignment = TaxonomyAssignment(
                    read_id=read_id,
                    taxid='unknown',  # Would lookup from taxonomy DB
                    species=species,
                    genus='unknown',  # Would lookup from taxonomy DB
                    confidence=min(mapq / 60.0, 1.0),  # Normalize MAPQ
                    alignment_length=len(parts[9]),
                    classified=True
                )
                
                assignments[read_id] = assignment
        
        return assignments
    
    def _extract_species_from_ref(self, ref_name: str) -> str:
        """Extract species name from reference sequence name."""
        # This is a simplified version
        # In practice, you'd have a mapping file from reference to taxonomy
        
        # Common formats:
        # - NC_000913.3 Escherichia coli str. K-12 substr. MG1655
        # - E.coli_K12
        
        parts = ref_name.split()
        if len(parts) >= 2:
            # Try to extract genus and species
            return f"{parts[0]} {parts[1]}"
        
        return ref_name
    
    def get_stats(self) -> Dict:
        """Get alignment statistics."""
        stats = self.stats.copy()
        stats['unique_species'] = len(self.stats['unique_species'])
        return stats


def main():
    """Test the microbe aligner."""
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    # Create test reads
    test_reads = [
        ReadRecord(
            read_id='read1',
            sequence='ATCGATCGATCGATCGATCGATCGATCG',
            quality='IIIIIIIIIIIIIIIIIIIIIIIIIIII',
            barcode='AAAA-1',
            umi='TTTT',
            mapq=0,
            is_unmapped=True,
            is_secondary=False,
            is_supplementary=False
        ),
        ReadRecord(
            read_id='read2',
            sequence='GGGGGGGGGGGGGGGGGGGGGGGGGGGG',
            quality='IIIIIIIIIIIIIIIIIIIIIIIIIIII',
            barcode='BBBB-1',
            umi='CCCC',
            mapq=0,
            is_unmapped=True,
            is_secondary=False,
            is_supplementary=False
        )
    ]
    
    print("MicrobeAligner test - would require reference database")
    print(f"Test reads: {len(test_reads)}")
    print("To run actual alignment, provide --microbe-ref or --kraken-db")


if __name__ == '__main__':
    main()
