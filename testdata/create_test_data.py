#!/usr/bin/env python3
"""
Create test data for scMicrobe pipeline.
"""

import os
import random
import pysam
from typing import List, Tuple


def create_test_host_genome(output_path: str, length: int = 100000):
    """Create a test host reference genome."""
    bases = ['A', 'T', 'C', 'G']
    
    with open(output_path, 'w') as f:
        f.write(">chr1\n")
        for i in range(0, length, 60):
            seq = ''.join(random.choices(bases, k=min(60, length - i)))
            f.write(seq + "\n")
    
    print(f"Created host genome: {output_path} ({length} bp)")


def create_test_microbe_genomes(output_dir: str):
    """Create test microbial reference genomes."""
    os.makedirs(output_dir, exist_ok=True)
    
    microbes = {
        'E.coli': 50000,
        'S.aureus': 40000,
        'B.subtilis': 45000,
        'P.aeruginosa': 55000,
    }
    
    bases = ['A', 'T', 'C', 'G']
    
    for name, length in microbes.items():
        output_path = os.path.join(output_dir, f"{name}.fa")
        with open(output_path, 'w') as f:
            f.write(f">{name}|chromosome|{length}bp\n")
            for i in range(0, length, 60):
                seq = ''.join(random.choices(bases, k=min(60, length - i)))
                f.write(seq + "\n")
        print(f"Created microbe genome: {output_path} ({length} bp)")


def create_test_bam(
    output_path: str,
    host_ref: str,
    microbe_refs: str,
    n_cells: int = 10,
    n_reads_per_cell: int = 100
):
    """
    Create a test BAM file with host and microbial reads.
    
    Args:
        output_path: Output BAM file path
        host_ref: Host reference genome path
        microbe_refs: Directory containing microbe reference genomes
        n_cells: Number of cells to simulate
        n_reads_per_cell: Number of reads per cell
    """
    # Read host genome
    host_seq = ""
    with open(host_ref, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                host_seq += line.strip()
    
    # Read microbe genomes
    microbe_seqs = {}
    for fname in os.listdir(microbe_refs):
        if fname.endswith('.fa'):
            name = fname.replace('.fa', '')
            path = os.path.join(microbe_refs, fname)
            seq = ""
            with open(path, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        seq += line.strip()
            microbe_seqs[name] = seq
    
    # Create header
    header = {
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'LN': len(host_seq), 'SN': 'chr1'}] + 
              [{'LN': len(seq), 'SN': name} for name, seq in microbe_seqs.items()]
    }
    
    # Create BAM file
    with pysam.AlignmentFile(output_path, 'wb', header=header) as bam:
        read_id = 0
        
        for cell_idx in range(n_cells):
            barcode = f"CELL{cell_idx:04d}"
            
            for _ in range(n_reads_per_cell):
                read = pysam.AlignedSegment()
                read.query_name = f"read_{read_id:08d}"
                
                # 70% host reads, 30% microbial
                if random.random() < 0.7:
                    # Host read - confidently mapped
                    start = random.randint(0, len(host_seq) - 100)
                    seq = host_seq[start:start+100]
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigarstring = "100M"
                else:
                    # Microbial read
                    microbe_name = random.choice(list(microbe_seqs.keys()))
                    microbe_seq = microbe_seqs[microbe_name]
                    start = random.randint(0, len(microbe_seq) - 100)
                    seq = microbe_seq[start:start+100]
                    
                    # Reference ID for microbe
                    ref_id = list(microbe_seqs.keys()).index(microbe_name) + 1
                    read.reference_id = ref_id
                    read.reference_start = start
                    read.mapping_quality = random.choice([0, 1, 5, 10])  # Low MAPQ
                    read.cigarstring = "100M"
                
                read.query_sequence = seq
                read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                
                # Add tags
                read.set_tag('CB', barcode)
                read.set_tag('UB', f"UMI{random.randint(1, 1000):04d}")
                
                bam.write(read)
                read_id += 1
    
    # Sort and index
    sorted_path = output_path.replace('.bam', '.sorted.bam')
    pysam.sort('-o', sorted_path, output_path)
    pysam.index(sorted_path)
    os.rename(sorted_path, output_path)
    
    print(f"Created test BAM: {output_path}")
    print(f"  Cells: {n_cells}")
    print(f"  Reads per cell: {n_reads_per_cell}")
    print(f"  Total reads: {read_id}")


def create_test_data(output_dir: str = "testdata"):
    """Create all test data."""
    os.makedirs(output_dir, exist_ok=True)
    
    print("Creating test data...")
    print("=" * 50)
    
    # Create host genome
    host_ref = os.path.join(output_dir, "host.fa")
    create_test_host_genome(host_ref)
    
    # Create microbe genomes
    microbe_dir = os.path.join(output_dir, "microbes")
    create_test_microbe_genomes(microbe_dir)
    
    # Create test BAM
    bam_path = os.path.join(output_dir, "test_sample.bam")
    create_test_bam(bam_path, host_ref, microbe_dir)
    
    print("=" * 50)
    print("Test data creation complete!")
    print(f"Output directory: {output_dir}")


if __name__ == '__main__':
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "testdata"
    create_test_data(output_dir)
