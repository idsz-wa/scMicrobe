"""
Host Filter module for scMicrobe.

Implements two-step host filtering:
1. Primary filter: Remove reads with MAPQ >= threshold
2. Secondary filter: Remove reads with high similarity to host genome
"""

import os
import logging
import tempfile
import subprocess
from typing import List, Dict, Optional
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor
import pysam

from bam_parser import ReadRecord


@dataclass
class AlignmentResult:
    """Result of alignment to host genome."""
    read_id: str
    identity: float
    aligned_length: int
    mismatches: int
    gaps: int
    is_host: bool


class HostFilter:
    """
    Two-step host filtering for microbial read identification.
    """
    
    def __init__(
        self,
        host_ref: str,
        mapq_threshold: int = 30,
        identity_threshold: float = 95.0,
        min_aligned_length: int = 50,
        threads: int = 4,
        aligner: str = 'minimap2'
    ):
        """
        Initialize host filter.
        
        Args:
            host_ref: Path to host reference genome (FASTA)
            mapq_threshold: MAPQ threshold for confident host alignment
            identity_threshold: Sequence identity threshold (%%) for host similarity
            min_aligned_length: Minimum aligned length to consider
            threads: Number of threads for alignment
            aligner: Alignment tool to use (minimap2 or bowtie2)
        """
        self.host_ref = host_ref
        self.mapq_threshold = mapq_threshold
        self.identity_threshold = identity_threshold
        self.min_aligned_length = min_aligned_length
        self.threads = threads
        self.aligner = aligner
        self.logger = logging.getLogger('HostFilter')
        
        # Statistics
        self.stats = {
            'primary_filtered': 0,
            'secondary_filtered': 0,
            'passed': 0
        }
        
        # Validate host reference
        if host_ref and not os.path.exists(host_ref):
            raise FileNotFoundError(f"Host reference not found: {host_ref}")
    
    def filter_reads(
        self,
        reads: List[ReadRecord],
        primary_only: bool = False
    ) -> List[ReadRecord]:
        """
        Apply two-step host filtering to reads.
        
        Args:
            reads: List of ReadRecord objects
            primary_only: Only apply primary filter (skip secondary alignment)
            
        Returns:
            List of non-host ReadRecord objects
        """
        self.logger.info(f"Starting host filtering on {len(reads)} reads")
        
        # Step 1: Primary filter (MAPQ-based)
        primary_passed = self._primary_filter(reads)
        self.stats['primary_filtered'] = len(reads) - len(primary_passed)
        self.logger.info(f"Primary filter: {self.stats['primary_filtered']} reads removed")
        
        if primary_only:
            return primary_passed
        
        # Step 2: Secondary filter (alignment-based)
        if self.host_ref and primary_passed:
            secondary_passed = self._secondary_filter(primary_passed)
            self.stats['secondary_filtered'] = len(primary_passed) - len(secondary_passed)
            self.logger.info(f"Secondary filter: {self.stats['secondary_filtered']} reads removed")
            return secondary_passed
        
        return primary_passed
    
    def _primary_filter(self, reads: List[ReadRecord]) -> List[ReadRecord]:
        """
        Primary filter: Remove reads confidently mapped to host.
        
        Reads with MAPQ >= threshold are considered host reads.
        """
        passed = []
        
        for read in reads:
            # If read is confidently mapped to host, filter it out
            if not read.is_unmapped and read.mapq >= self.mapq_threshold:
                continue
            passed.append(read)
        
        return passed
    
    def _secondary_filter(self, reads: List[ReadRecord]) -> List[ReadRecord]:
        """
        Secondary filter: Re-align to host and filter by identity.
        
        Reads showing high similarity (> threshold) to host are filtered out.
        """
        if not reads:
            return []
        
        # Write reads to temporary FASTQ
        with tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.fastq',
            delete=False
        ) as tmp_fastq:
            tmp_fastq_path = tmp_fastq.name
            for read in reads:
                tmp_fastq.write(read.to_fastq())
        
        try:
            # Align to host
            alignment_results = self._align_to_host(tmp_fastq_path)
            
            # Filter reads based on alignment results
            passed = []
            for read in reads:
                result = alignment_results.get(read.read_id)
                
                if result is None:
                    # No alignment to host - keep as candidate
                    passed.append(read)
                elif not result.is_host:
                    # Low similarity to host - keep as candidate
                    passed.append(read)
            
            return passed
            
        finally:
            # Cleanup
            if os.path.exists(tmp_fastq_path):
                os.unlink(tmp_fastq_path)
    
    def _align_to_host(self, fastq_path: str) -> Dict[str, AlignmentResult]:
        """
        Align reads to host genome and calculate identity.
        
        Returns:
            Dictionary mapping read_id to AlignmentResult
        """
        results = {}
        
        if self.aligner == 'minimap2':
            results = self._align_with_minimap2(fastq_path)
        elif self.aligner == 'bowtie2':
            results = self._align_with_bowtie2(fastq_path)
        else:
            raise ValueError(f"Unsupported aligner: {self.aligner}")
        
        return results
    
    def _align_with_minimap2(self, fastq_path: str) -> Dict[str, AlignmentResult]:
        """Align reads using minimap2."""
        results = {}
        
        # Build minimap2 command
        cmd = [
            'minimap2',
            '-ax', 'sr',  # Short reads
            '-t', str(self.threads),
            '--secondary=no',
            self.host_ref,
            fastq_path
        ]
        
        try:
            # Run minimap2 and capture SAM output
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            # Parse SAM output
            for line in process.stdout:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                
                read_id = parts[0]
                flag = int(parts[1])
                
                # Skip unmapped reads
                if flag & 0x4:  # Unmapped flag
                    continue
                
                # Calculate identity from CIGAR and NM tag
                cigar = parts[5]
                seq_len = len(parts[9])
                
                # Parse NM tag (edit distance)
                nm = 0
                for tag in parts[11:]:
                    if tag.startswith('NM:i:'):
                        nm = int(tag.split(':')[2])
                        break
                
                # Calculate identity
                matches = seq_len - nm
                identity = (matches / seq_len * 100) if seq_len > 0 else 0
                
                # Determine if host based on identity threshold
                is_host = identity >= self.identity_threshold
                
                results[read_id] = AlignmentResult(
                    read_id=read_id,
                    identity=identity,
                    aligned_length=seq_len,
                    mismatches=nm,
                    gaps=0,
                    is_host=is_host
                )
            
            process.wait()
            
        except FileNotFoundError:
            self.logger.error("minimap2 not found. Please install minimap2.")
            raise
        except Exception as e:
            self.logger.error(f"minimap2 alignment failed: {e}")
            raise
        
        return results
    
    def _align_with_bowtie2(self, fastq_path: str) -> Dict[str, AlignmentResult]:
        """Align reads using bowtie2."""
        # Similar implementation for bowtie2
        # This is a placeholder - full implementation would follow similar pattern
        self.logger.warning("bowtie2 alignment not fully implemented, using minimap2")
        return self._align_with_minimap2(fastq_path)
    
    def get_stats(self) -> Dict:
        """Get filtering statistics."""
        return self.stats.copy()


def main():
    """Test the host filter."""
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    # Create test reads
    test_reads = [
        ReadRecord(
            read_id='read1',
            sequence='ATCGATCGATCGATCGATCGATCG',
            quality='IIIIIIIIIIIIIIIIIIIIIIII',
            barcode='AAAA-1',
            umi='TTTT',
            mapq=40,  # High MAPQ - should be filtered
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False
        ),
        ReadRecord(
            read_id='read2',
            sequence='GGGGGGGGGGGGGGGGGGGGGGGG',
            quality='IIIIIIIIIIIIIIIIIIIIIIII',
            barcode='AAAA-1',
            umi='CCCC',
            mapq=1,  # Low MAPQ - should pass primary
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False
        ),
        ReadRecord(
            read_id='read3',
            sequence='TTTTTTTTTTTTTTTTTTTTTTTT',
            quality='IIIIIIIIIIIIIIIIIIIIIIII',
            barcode='BBBB-1',
            umi='AAAA',
            mapq=0,  # Unmapped
            is_unmapped=True,
            is_secondary=False,
            is_supplementary=False
        )
    ]
    
    # Test primary filter only
    filter_primary = HostFilter(
        host_ref=None,
        mapq_threshold=30,
        primary_only=True
    )
    
    passed = filter_primary.filter_reads(test_reads)
    print(f"\nPrimary filter test:")
    print(f"  Input: {len(test_reads)} reads")
    print(f"  Passed: {len(passed)} reads")
    print(f"  Stats: {filter_primary.get_stats()}")


if __name__ == '__main__':
    main()
