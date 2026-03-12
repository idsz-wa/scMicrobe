#!/usr/bin/env python3
"""
scMicrobe: Single-Cell Microbial Quantification from Host scRNA-seq BAM

Main entry point for the scMicrobe pipeline.
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Optional
import yaml

from bam_parser import BAMParser
from host_filter import HostFilter
from microbe_aligner import MicrobeAligner
from quantifier import Quantifier
from contamination_filter import ContaminationFilter
from output_writer import OutputWriter


def setup_logging(log_level: str = "INFO"):
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('scmicro.log')
        ]
    )


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='scMicrobe: Single-Cell Microbial Quantification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with minimap2
  python scmicro.py --input sample.bam --output out/ --host-ref host.fa --microbe-ref microbes.fa
  
  # Using Kraken2 for classification
  python scmicro.py --input sample.bam --output out/ --aligner kraken2 --kraken-db /path/to/db
  
  # With custom barcode and UMI tags
  python scmicro.py --input sample.bam --output out/ --barcode-tag CR --umi-tag UR
        """
    )
    
    # Input/Output
    parser.add_argument('--input', '-i', required=True, help='Input BAM file')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--config', '-c', help='Configuration YAML file')
    
    # Reference files
    parser.add_argument('--host-ref', help='Host reference genome (FASTA)')
    parser.add_argument('--microbe-ref', help='Microbial reference genome (FASTA)')
    parser.add_argument('--kraken-db', help='Kraken2 database path')
    
    # Alignment options
    parser.add_argument('--aligner', choices=['minimap2', 'kraken2'], 
                       default='minimap2', help='Alignment method (default: minimap2)')
    parser.add_argument('--minimap2-index', help='Pre-built minimap2 index')
    
    # BAM tags
    parser.add_argument('--barcode-tag', default='CB', 
                       help='Cell barcode tag in BAM (default: CB)')
    parser.add_argument('--umi-tag', default='UB', 
                       help='UMI tag in BAM (default: UB)')
    
    # Filtering options
    parser.add_argument('--host-mapq', type=int, default=30,
                       help='Host MAPQ threshold for primary filtering (default: 30)')
    parser.add_argument('--host-identity', type=float, default=95.0,
                       help='Host sequence identity threshold %% (default: 95)')
    parser.add_argument('--min-microbe-mapq', type=int, default=1,
                       help='Minimum MAPQ for microbial alignment (default: 1)')
    
    # Quantification options
    parser.add_argument('--quant-method', choices=['reads', 'umi'], default='umi',
                       help='Quantification method (default: umi)')
    parser.add_argument('--umi-dedup-method', choices=['exact', 'hamming'], default='exact',
                       help='UMI deduplication method (default: exact)')
    
    # Contamination correction
    parser.add_argument('--ambient-cutoff', type=float, default=0.1,
                       help='Ambient RNA cutoff quantile (default: 0.1)')
    parser.add_argument('--min-cells', type=int, default=3,
                       help='Minimum cells for microbe detection (default: 3)')
    parser.add_argument('--skip-contamination-filter', action='store_true',
                       help='Skip contamination filtering')
    
    # Performance
    parser.add_argument('--threads', '-t', type=int, default=4,
                       help='Number of threads (default: 4)')
    parser.add_argument('--chunk-size', type=int, default=10000,
                       help='Reads per chunk for processing (default: 10000)')
    
    # Logging
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level (default: INFO)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    return parser.parse_args()


def load_config(config_path: Optional[str]) -> dict:
    """Load configuration from YAML file."""
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    return {}


def validate_args(args) -> bool:
    """Validate command line arguments."""
    if not os.path.exists(args.input):
        logging.error(f"Input BAM file not found: {args.input}")
        return False
    
    if args.aligner == 'minimap2' and not args.microbe_ref and not args.minimap2_index:
        logging.error("Minimap2 alignment requires --microbe-ref or --minimap2-index")
        return False
    
    if args.aligner == 'kraken2' and not args.kraken_db:
        logging.error("Kraken2 alignment requires --kraken-db")
        return False
    
    return True


def main():
    """Main pipeline execution."""
    args = parse_args()
    setup_logging(args.log_level)
    logger = logging.getLogger('scMicrobe')
    
    logger.info("Starting scMicrobe pipeline")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    # Load config if provided
    config = load_config(args.config)
    
    # Validate arguments
    if not validate_args(args):
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Step 1: Parse BAM and extract reads
    logger.info("Step 1: Parsing BAM file...")
    bam_parser = BAMParser(
        bam_path=args.input,
        barcode_tag=args.barcode_tag,
        umi_tag=args.umi_tag,
        chunk_size=args.chunk_size
    )
    
    # Step 2: Host filtering
    logger.info("Step 2: Host filtering...")
    host_filter = HostFilter(
        host_ref=args.host_ref,
        mapq_threshold=args.host_mapq,
        identity_threshold=args.host_identity,
        threads=args.threads
    )
    
    # Step 3: Microbial alignment
    logger.info(f"Step 3: Microbial alignment using {args.aligner}...")
    microbe_aligner = MicrobeAligner(
        aligner=args.aligner,
        microbe_ref=args.microbe_ref,
        kraken_db=args.kraken_db,
        minimap2_index=args.minimap2_index,
        min_mapq=args.min_microbe_mapq,
        threads=args.threads
    )
    
    # Step 4: Quantification
    logger.info("Step 4: Quantification...")
    quantifier = Quantifier(
        method=args.quant_method,
        umi_dedup_method=args.umi_dedup_method
    )
    
    # Step 5: Contamination filtering
    if not args.skip_contamination_filter:
        logger.info("Step 5: Contamination filtering...")
        contamination_filter = ContaminationFilter(
            ambient_cutoff=args.ambient_cutoff,
            min_cells=args.min_cells
        )
    else:
        contamination_filter = None
    
    # Step 6: Output writing
    logger.info("Step 6: Writing output...")
    output_writer = OutputWriter(output_dir=args.output)
    
    # Execute pipeline
    try:
        # Extract candidate reads from BAM
        candidate_reads = bam_parser.extract_candidate_reads(
            mapq_threshold=args.host_mapq
        )
        logger.info(f"Extracted {len(candidate_reads)} candidate reads")
        
        # Host filtering
        filtered_reads = host_filter.filter_reads(candidate_reads)
        logger.info(f"After host filtering: {len(filtered_reads)} reads")
        
        # Microbial alignment
        aligned_reads = microbe_aligner.align(filtered_reads)
        logger.info(f"After microbial alignment: {len(aligned_reads)} reads")
        
        # Quantification
        microbe_cell_matrix = quantifier.quantify(aligned_reads)
        logger.info(f"Matrix shape: {microbe_cell_matrix.shape}")
        
        # Contamination filtering
        if contamination_filter:
            microbe_cell_matrix = contamination_filter.filter(microbe_cell_matrix)
            logger.info(f"After contamination filter: {microbe_cell_matrix.shape}")
        
        # Write output
        output_writer.write_matrix(microbe_cell_matrix)
        output_writer.write_stats({
            'total_reads': len(candidate_reads),
            'host_filtered': len(filtered_reads),
            'microbial_aligned': len(aligned_reads),
            'final_matrix_shape': microbe_cell_matrix.shape
        })
        
        logger.info("Pipeline completed successfully!")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
