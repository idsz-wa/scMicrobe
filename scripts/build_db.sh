#!/bin/bash
# Build microbial reference database for scMicrobe

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="${1:-./microbe_db}"
HOST_REF="${2:-}"

mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "Building microbial reference database in $DB_DIR"

# Step 1: Download NCBI RefSeq genomes
echo "Step 1: Downloading microbial genomes..."

# Bacteria
if [ ! -f "bacteria.complete.list" ]; then
    echo "Fetching bacteria assembly summary..."
    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O bacteria.assembly_summary.txt
    awk -F '\t' '$12=="Complete Genome" && $11=="latest" {print $20}' bacteria.assembly_summary.txt > bacteria.complete.list
fi

# Fungi
if [ ! -f "fungi.complete.list" ]; then
    echo "Fetching fungi assembly summary..."
    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt -O fungi.assembly_summary.txt
    awk -F '\t' '$12=="Complete Genome" && $11=="latest" {print $20}' fungi.assembly_summary.txt > fungi.complete.list
fi

# Viruses
if [ ! -f "viral.complete.list" ]; then
    echo "Fetching viral assembly summary..."
    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt -O viral.assembly_summary.txt
    awk -F '\t' '$12=="Complete Genome" && $11=="latest" {print $20}' viral.assembly_summary.txt > viral.complete.list
fi

# Step 2: Download genomes
echo "Step 2: Downloading genome sequences..."
mkdir -p genomes

for url in $(cat bacteria.complete.list fungi.complete.list viral.complete.list | head -n 1000); do
    acc=$(basename "$url")
    if [ ! -f "genomes/${acc}_genomic.fna.gz" ]; then
        wget -q -P genomes "${url}/${acc}_genomic.fna.gz" || true
    fi
done

# Step 3: Remove host-like sequences (if host reference provided)
if [ -n "$HOST_REF" ] && [ -f "$HOST_REF" ]; then
    echo "Step 3: Removing host-like sequences..."
    
    # Build minimap2 index for host
    if [ ! -f "${HOST_REF}.mmi" ]; then
        echo "Building host index..."
        minimap2 -d "${HOST_REF}.mmi" "$HOST_REF"
    fi
    
    # Filter each genome
    mkdir -p genomes_filtered
    for genome in genomes/*.fna.gz; do
        name=$(basename "$genome" .gz)
        if [ ! -f "genomes_filtered/${name}" ]; then
            echo "Filtering $name..."
            zcat "$genome" | minimap2 -t 4 "${HOST_REF}.mmi" - \
                | awk '$5 < 95 {print $1}' | sort -u > "genomes_filtered/${name}.host_contigs"
            
            # Extract non-host contigs
            zcat "$genome" | seqkit grep -v -f "genomes_filtered/${name}.host_contigs" \
                > "genomes_filtered/${name}" || cp <(zcat "$genome") "genomes_filtered/${name}"
        fi
    done
    
    GENOME_DIR="genomes_filtered"
else
    echo "Step 3: Skipping host filtering (no host reference provided)"
    GENOME_DIR="genomes"
fi

# Step 4: Build combined FASTA
echo "Step 4: Building combined reference..."
if [ ! -f "microbe_combined.fa" ]; then
    cat "$GENOME_DIR"/*.fna > microbe_combined.fa 2>/dev/null || true
fi

# Step 5: Build minimap2 index
echo "Step 5: Building minimap2 index..."
if [ ! -f "microbe_combined.mmi" ] && [ -f "microbe_combined.fa" ]; then
    minimap2 -d microbe_combined.mmi microbe_combined.fa
fi

# Step 6: Build Kraken2 database (optional)
if command -v kraken2-build &> /dev/null; then
    echo "Step 6: Building Kraken2 database..."
    mkdir -p kraken2_db
    kraken2-build --download-taxonomy --db kraken2_db || true
    kraken2-build --add-to-library microbe_combined.fa --db kraken2_db || true
    kraken2-build --build --db kraken2_db --threads 4 || true
fi

echo "Database build complete!"
echo "  Combined FASTA: $DB_DIR/microbe_combined.fa"
echo "  Minimap2 index: $DB_DIR/microbe_combined.mmi"
