"""
BAM Parser module for scMicrobe.

Extracts cell barcodes, UMIs, and candidate reads from scRNA-seq BAM files.
"""

import logging
from typing import Iterator, Dict, List, Optional, Tuple
from dataclasses import dataclass
import pysam
from tqdm import tqdm


@dataclass
class ReadRecord:
    """Represents a single read record with metadata."""
    read_id: str
    sequence: str
    quality: str
    barcode: Optional[str]
    umi: Optional[str]
    mapq: int
    is_unmapped: bool
    is_secondary: bool
    is_supplementary: bool
    reference_name: Optional[str] = None
    reference_start: Optional[int] = None
    cigar: Optional[str] = None
    
    def to_fasta(self) -> str:
        """Convert to FASTA format."""
        return f">{self.read_id}\n{self.sequence}\n"
    
    def to_fastq(self) -> str:
        """Convert to FASTQ format."""
        return f"@{self.read_id}\n{self.sequence}\n+\n{self.quality}\n"


class BAMParser:
    """Parser for scRNA-seq BAM files."""
    
    def __init__(
        self,
        bam_path: str,
        barcode_tag: str = 'CB',
        umi_tag: str = 'UB',
        chunk_size: int = 10000
    ):
        """
        Initialize BAM parser.
        
        Args:
            bam_path: Path to input BAM file
            barcode_tag: BAM tag for cell barcode (default: CB)
            umi_tag: BAM tag for UMI (default: UB)
            chunk_size: Number of reads per chunk for processing
        """
        self.bam_path = bam_path
        self.barcode_tag = barcode_tag
        self.umi_tag = umi_tag
        self.chunk_size = chunk_size
        self.logger = logging.getLogger('BAMParser')
        
        # Statistics
        self.stats = {
            'total_reads': 0,
            'mapped_reads': 0,
            'unmapped_reads': 0,
            'low_mapq_reads': 0,
            'reads_with_barcode': 0,
            'reads_with_umi': 0
        }
    
    def _extract_barcode(self, read) -> Optional[str]:
        """Extract cell barcode from read."""
        if self.barcode_tag and read.has_tag(self.barcode_tag):
            barcode = read.get_tag(self.barcode_tag)
            # Remove '-1' suffix if present
            if barcode and '-' in barcode:
                barcode = barcode.split('-')[0]
            return barcode
        return None
    
    def _extract_umi(self, read) -> Optional[str]:
        """Extract UMI from read."""
        if self.umi_tag and read.has_tag(self.umi_tag):
            return read.get_tag(self.umi_tag)
        return None
    
    def _read_to_record(self, read) -> ReadRecord:
        """Convert pysam read to ReadRecord."""
        return ReadRecord(
            read_id=read.query_name,
            sequence=read.query_sequence or '',
            quality=pysam.qualities_to_qualitystring(read.query_qualities) if read.query_qualities else '',
            barcode=self._extract_barcode(read),
            umi=self._extract_umi(read),
            mapq=read.mapping_quality,
            is_unmapped=read.is_unmapped,
            is_secondary=read.is_secondary,
            is_supplementary=read.is_supplementary,
            reference_name=read.reference_name if not read.is_unmapped else None,
            reference_start=read.reference_start if not read.is_unmapped else None,
            cigar=read.cigarstring if not read.is_unmapped else None
        )
    
    def extract_candidate_reads(
        self,
        mapq_threshold: int = 30,
        include_unmapped: bool = True,
        include_low_mapq: bool = True
    ) -> List[ReadRecord]:
        """
        Extract candidate microbial reads from BAM.
        
        Candidate reads are those that:
        1. Are unmapped, OR
        2. Have low MAPQ (< threshold), OR
        3. Are not confidently mapped to host
        
        Args:
            mapq_threshold: MAPQ threshold for confident host alignment
            include_unmapped: Include unmapped reads
            include_low_mapq: Include low MAPQ reads
            
        Returns:
            List of ReadRecord objects
        """
        candidates = []
        
        self.logger.info(f"Parsing BAM file: {self.bam_path}")
        
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            # Try to get total reads, fallback to None if index not available
            try:
                total_reads = bam.mapped + bam.unmapped
            except ValueError:
                total_reads = None
            
            for read in tqdm(bam, total=total_reads, desc="Parsing BAM"):
                self.stats['total_reads'] += 1
                
                # Skip secondary and supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue
                
                # Track statistics
                if read.is_unmapped:
                    self.stats['unmapped_reads'] += 1
                else:
                    self.stats['mapped_reads'] += 1
                    if read.mapping_quality < mapq_threshold:
                        self.stats['low_mapq_reads'] += 1
                
                # Check if read has barcode
                barcode = self._extract_barcode(read)
                if barcode:
                    self.stats['reads_with_barcode'] += 1
                
                # Check if read has UMI
                if self._extract_umi(read):
                    self.stats['reads_with_umi'] += 1
                
                # Determine if this is a candidate read
                is_candidate = False
                
                if include_unmapped and read.is_unmapped:
                    is_candidate = True
                elif include_low_mapq and not read.is_unmapped and read.mapping_quality < mapq_threshold:
                    is_candidate = True
                
                if is_candidate:
                    record = self._read_to_record(read)
                    candidates.append(record)
        
        self._log_stats()
        return candidates
    
    def extract_reads_by_barcode(
        self,
        barcodes: List[str]
    ) -> Iterator[List[ReadRecord]]:
        """
        Extract reads for specific cell barcodes.
        
        Args:
            barcodes: List of cell barcodes to extract
            
        Yields:
            Chunks of ReadRecord objects
        """
        barcode_set = set(barcodes)
        chunk = []
        
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            for read in bam:
                barcode = self._extract_barcode(read)
                if barcode and barcode in barcode_set:
                    chunk.append(self._read_to_record(read))
                    
                    if len(chunk) >= self.chunk_size:
                        yield chunk
                        chunk = []
        
        if chunk:
            yield chunk
    
    def get_unique_barcodes(self) -> set:
        """Get all unique cell barcodes from BAM file."""
        barcodes = set()
        
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            for read in bam:
                barcode = self._extract_barcode(read)
                if barcode:
                    barcodes.add(barcode)
        
        return barcodes
    
    def _log_stats(self):
        """Log parsing statistics."""
        self.logger.info("BAM Parsing Statistics:")
        self.logger.info(f"  Total reads: {self.stats['total_reads']:,}")
        self.logger.info(f"  Mapped reads: {self.stats['mapped_reads']:,}")
        self.logger.info(f"  Unmapped reads: {self.stats['unmapped_reads']:,}")
        self.logger.info(f"  Low MAPQ reads: {self.stats['low_mapq_reads']:,}")
        self.logger.info(f"  Reads with barcode: {self.stats['reads_with_barcode']:,}")
        self.logger.info(f"  Reads with UMI: {self.stats['reads_with_umi']:,}")
    
    def write_to_fastq(
        self,
        reads: List[ReadRecord],
        output_path: str
    ):
        """Write reads to FASTQ file."""
        with open(output_path, 'w') as f:
            for read in reads:
                f.write(read.to_fastq())
        self.logger.info(f"Written {len(reads)} reads to {output_path}")


def main():
    """Test the BAM parser."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python bam_parser.py <bam_file>")
        sys.exit(1)
    
    logging.basicConfig(level=logging.INFO)
    
    parser = BAMParser(sys.argv[1])
    candidates = parser.extract_candidate_reads()
    
    print(f"Found {len(candidates)} candidate reads")
    
    if candidates:
        print("\nFirst 5 candidates:")
        for i, read in enumerate(candidates[:5]):
            print(f"  {i+1}. {read.read_id} | Barcode: {read.barcode} | UMI: {read.umi}")


if __name__ == '__main__':
    main()
