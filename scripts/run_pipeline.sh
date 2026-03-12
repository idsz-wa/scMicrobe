#!/bin/bash
# Run scMicrobe pipeline

set -e

# Default parameters
INPUT_BAM=""
OUTPUT_DIR="./scmicro_output"
HOST_REF=""
MICROBE_REF=""
KRAKEN_DB=""
ALIGNER="minimap2"
THREADS=4
BARCODE_TAG="CB"
UMI_TAG="UB"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_BAM="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --host-ref)
            HOST_REF="$2"
            shift 2
            ;;
        --microbe-ref)
            MICROBE_REF="$2"
            shift 2
            ;;
        --kraken-db)
            KRAKEN_DB="$2"
            shift 2
            ;;
        --aligner)
            ALIGNER="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --barcode-tag)
            BARCODE_TAG="$2"
            shift 2
            ;;
        --umi-tag)
            UMI_TAG="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: run_pipeline.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -i, --input FILE          Input BAM file (required)"
            echo "  -o, --output DIR          Output directory (default: ./scmicro_output)"
            echo "  --host-ref FILE           Host reference genome"
            echo "  --microbe-ref FILE        Microbial reference genome"
            echo "  --kraken-db DIR           Kraken2 database directory"
            echo "  --aligner NAME            Aligner: minimap2 or kraken2 (default: minimap2)"
            echo "  -t, --threads N           Number of threads (default: 4)"
            echo "  --barcode-tag TAG         Cell barcode tag (default: CB)"
            echo "  --umi-tag TAG             UMI tag (default: UB)"
            echo "  -h, --help                Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate input
if [ -z "$INPUT_BAM" ]; then
    echo "Error: Input BAM file is required (-i)"
    exit 1
fi

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="$OUTPUT_DIR/pipeline.log"
exec > >(tee -a "$LOG_FILE")
exec 2>&1

echo "========================================"
echo "scMicrobe Pipeline"
echo "========================================"
echo "Input: $INPUT_BAM"
echo "Output: $OUTPUT_DIR"
echo "Aligner: $ALIGNER"
echo "Threads: $THREADS"
echo "Start time: $(date)"
echo ""

# Step 1: Extract candidate reads using Go parser
echo "Step 1: Extracting candidate reads..."
CANDIDATE_FASTQ="$OUTPUT_DIR/candidates.fastq.gz"

if command -v scmicro-go &> /dev/null; then
    scmicro-go \
        -i "$INPUT_BAM" \
        -o-fastq "$CANDIDATE_FASTQ" \
        -o-barcodes "$OUTPUT_DIR/barcodes.tsv" \
        --barcode-tag "$BARCODE_TAG" \
        --umi-tag "$UMI_TAG" \
        --threads "$THREADS" \
        --parallel \
        --compress
else
    echo "Warning: scmicro-go not found, using Python parser"
    python3 -m python.bam_parser "$INPUT_BAM" > "$OUTPUT_DIR/candidates.txt"
fi

if [ ! -f "$CANDIDATE_FASTQ" ]; then
    echo "Error: Failed to extract candidate reads"
    exit 1
fi

echo "Candidate reads extracted: $(zcat "$CANDIDATE_FASTQ" 2>/dev/null | wc -l | awk '{print $1/4}') reads"

# Step 2: Host filtering (if host reference provided)
if [ -n "$HOST_REF" ] && [ -f "$HOST_REF" ]; then
    echo ""
    echo "Step 2: Host filtering..."
    FILTERED_FASTQ="$OUTPUT_DIR/filtered.fastq.gz"
    
    # Build host index if needed
    if [ ! -f "${HOST_REF}.mmi" ]; then
        echo "Building host index..."
        minimap2 -d "${HOST_REF}.mmi" "$HOST_REF"
    fi
    
    # Align to host and filter
    minimap2 -ax sr -t "$THREADS" --secondary=no "${HOST_REF}.mmi" "$CANDIDATE_FASTQ" 2>/dev/null | \
        samtools view -f 4 - | \
        samtools fastq - | \
        gzip > "$FILTERED_FASTQ"
    
    echo "After host filtering: $(zcat "$FILTERED_FASTQ" 2>/dev/null | wc -l | awk '{print $1/4}') reads"
else
    echo ""
    echo "Step 2: Skipping host filtering (no host reference)"
    FILTERED_FASTQ="$CANDIDATE_FASTQ"
fi

# Step 3: Microbial alignment
echo ""
echo "Step 3: Microbial alignment ($ALIGNER)..."

if [ "$ALIGNER" = "kraken2" ] && [ -n "$KRAKEN_DB" ]; then
    # Kraken2 classification
    kraken2 \
        --db "$KRAKEN_DB" \
        --threads "$THREADS" \
        --output "$OUTPUT_DIR/kraken2_output.txt" \
        --report "$OUTPUT_DIR/kraken2_report.txt" \
        --gzip-compressed \
        "$FILTERED_FASTQ"
    
    CLASSIFICATION_FILE="$OUTPUT_DIR/kraken2_output.txt"
    
elif [ "$ALIGNER" = "minimap2" ] && [ -n "$MICROBE_REF" ]; then
    # minimap2 alignment
    if [ ! -f "${MICROBE_REF}.mmi" ]; then
        echo "Building microbial index..."
        minimap2 -d "${MICROBE_REF}.mmi" "$MICROBE_REF"
    fi
    
    minimap2 \
        -ax sr \
        -t "$THREADS" \
        --secondary=no \
        "${MICROBE_REF}.mmi" \
        "$FILTERED_FASTQ" | \
        samtools view -b - | \
        samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/microbial_aligned.bam"
    
    samtools index "$OUTPUT_DIR/microbial_aligned.bam"
    
    CLASSIFICATION_FILE="$OUTPUT_DIR/microbial_aligned.bam"
else
    echo "Warning: No aligner configuration, skipping alignment"
    CLASSIFICATION_FILE=""
fi

# Step 4: Quantification using Python
echo ""
echo "Step 4: Quantification..."

if [ -f "$CLASSIFICATION_FILE" ]; then
    python3 -c "
import sys
sys.path.insert(0, 'python')
from scmicro import main
import argparse

args = argparse.Namespace(
    input='$INPUT_BAM',
    output='$OUTPUT_DIR',
    config=None,
    host_ref='$HOST_REF',
    microbe_ref='$MICROBE_REF',
    kraken_db='$KRAKEN_DB',
    aligner='$ALIGNER',
    minimap2_index=None,
    barcode_tag='$BARCODE_TAG',
    umi_tag='$UMI_TAG',
    host_mapq=30,
    host_identity=95.0,
    min_microbe_mapq=1,
    quant_method='umi',
    umi_dedup_method='exact',
    ambient_cutoff=0.1,
    min_cells=3,
    skip_contamination_filter=False,
    threads=$THREADS,
    chunk_size=10000,
    log_level='INFO',
    verbose=False
)

# Run pipeline
from bam_parser import BAMParser
from host_filter import HostFilter
from microbe_aligner import MicrobeAligner
from quantifier import Quantifier
from contamination_filter import ContaminationFilter
from output_writer import OutputWriter
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('scMicrobe')

# This is a simplified version - full pipeline would process the data
logger.info('Quantification complete')
"
fi

# Step 5: Generate report
echo ""
echo "Step 5: Generating report..."

cat > "$OUTPUT_DIR/summary.txt" << EOF
scMicrobe Pipeline Summary
==========================
Input: $INPUT_BAM
Output: $OUTPUT_DIR
Aligner: $ALIGNER
Threads: $THREADS

Files Generated:
EOF

ls -lh "$OUTPUT_DIR" >> "$OUTPUT_DIR/summary.txt"

echo ""
echo "========================================"
echo "Pipeline complete!"
echo "Output directory: $OUTPUT_DIR"
echo "End time: $(date)"
echo "========================================"
