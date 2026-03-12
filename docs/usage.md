# scMicrobe Usage Guide

## Installation

### Python Dependencies

```bash
cd scmicro/python
pip install -r requirements.txt
```

### Go Components (Optional, for better performance)

```bash
cd scmicro/go
go build -o scmicro-go .
# Or use make
make build
```

### External Tools

scMicrobe requires the following external tools:

- `minimap2` - For sequence alignment
- `samtools` - For BAM manipulation
- `kraken2` (optional) - For fast taxonomic classification

Install via conda:
```bash
conda install -c bioconda minimap2 samtools kraken2
```

## Quick Start

### 1. Build Microbial Reference Database

```bash
./scripts/build_db.sh ./microbe_db /path/to/host.fa
```

### 2. Run Pipeline

```bash
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --host-ref host.fa \
    --microbe-ref microbe_db/microbe_combined.fa \
    --aligner minimap2 \
    --threads 8
```

Or use the shell script:

```bash
./scripts/run_pipeline.sh \
    -i sample.bam \
    -o ./output \
    --host-ref host.fa \
    --microbe-ref microbe_db/microbe_combined.fa \
    -t 8
```

## Command Line Options

### Input/Output

| Option | Description | Default |
|--------|-------------|---------|
| `--input`, `-i` | Input BAM file (required) | - |
| `--output`, `-o` | Output directory | `./scmicro_output` |
| `--config`, `-c` | Configuration YAML file | - |

### Reference Files

| Option | Description |
|--------|-------------|
| `--host-ref` | Host reference genome (FASTA) |
| `--microbe-ref` | Microbial reference genome (FASTA) |
| `--kraken-db` | Kraken2 database path |

### Alignment Options

| Option | Description | Default |
|--------|-------------|---------|
| `--aligner` | Alignment method: `minimap2` or `kraken2` | `minimap2` |
| `--minimap2-index` | Pre-built minimap2 index | - |

### BAM Tags

| Option | Description | Default |
|--------|-------------|---------|
| `--barcode-tag` | Cell barcode tag in BAM | `CB` |
| `--umi-tag` | UMI tag in BAM | `UB` |

### Filtering Options

| Option | Description | Default |
|--------|-------------|---------|
| `--host-mapq` | Host MAPQ threshold | `30` |
| `--host-identity` | Host sequence identity threshold (%) | `95` |
| `--min-microbe-mapq` | Minimum MAPQ for microbial alignment | `1` |

### Quantification Options

| Option | Description | Default |
|--------|-------------|---------|
| `--quant-method` | Quantification: `reads` or `umi` | `umi` |
| `--umi-dedup-method` | UMI dedup: `exact` or `hamming` | `exact` |

### Contamination Correction

| Option | Description | Default |
|--------|-------------|---------|
| `--ambient-cutoff` | Ambient RNA cutoff quantile | `0.1` |
| `--min-cells` | Minimum cells for microbe detection | `3` |
| `--skip-contamination-filter` | Skip contamination filtering | `False` |

### Performance

| Option | Description | Default |
|--------|-------------|---------|
| `--threads`, `-t` | Number of threads | `4` |
| `--chunk-size` | Reads per chunk | `10000` |

## Output Files

The pipeline generates the following output files:

### Matrix Files

- `scmicro_matrix.tsv` - Dense matrix in TSV format (microbes × cells)
- `mtx/matrix.mtx` - Matrix Market format (for Seurat/Scanpy)
- `mtx/barcodes.tsv` - Cell barcodes
- `mtx/features.tsv` - Microbe features
- `scmicro.h5ad` - H5AD format (for Scanpy)

### Statistics

- `scmicro_stats.json` - Pipeline statistics
- `scmicro_qc.txt` - QC report

### Intermediate Files

- `candidates.fastq.gz` - Extracted candidate reads
- `filtered.fastq.gz` - Reads after host filtering
- `microbial_aligned.bam` - Microbial alignment (if using minimap2)

## Example Workflows

### Basic Analysis

```bash
# Run scMicrobe
python python/scmicro.py \
    --input pbmc_10x.bam \
    --output ./pbmc_output \
    --host-ref hg38.fa \
    --microbe-ref microbe_combined.fa \
    --threads 16

# Load in Scanpy
import scanpy as sc
adata = sc.read_h5ad('./pbmc_output/scmicro.h5ad')
sc.pp.filter_cells(adata, min_counts=1)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['species'])
```

### Using Kraken2 for Fast Classification

```bash
# Build Kraken2 database first
kraken2-build --download-taxonomy --db kraken2_db
kraken2-build --download-library bacteria --db kraken2_db
kraken2-build --build --db kraken2_db --threads 8

# Run with Kraken2
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --aligner kraken2 \
    --kraken-db kraken2_db \
    --threads 8
```

### Custom Configuration

Create a `config.yaml` file:

```yaml
host_filter:
  mapq_threshold: 30
  identity_threshold: 95

quantification:
  method: umi
  umi_dedup_method: exact

contamination:
  ambient_cutoff: 0.1
  min_cells: 3
  specificity_threshold: 2.0

output:
  formats:
    - mtx
    - h5ad
    - tsv
```

Run with config:

```bash
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --config config.yaml
```

## Troubleshooting

### Low microbial detection rate

- Check host filtering parameters (--host-mapq, --host-identity)
- Verify microbial reference database completeness
- Consider using Kraken2 for more sensitive detection

### High contamination

- Enable contamination filtering (default)
- Adjust --ambient-cutoff based on your data
- Check for reagent contamination in negative controls

### Memory issues

- Reduce --chunk-size
- Use Go-based BAM parser for lower memory footprint
- Process data in batches

### Performance optimization

- Use Go-based parser: `scmicro-go -i input.bam`
- Increase --threads
- Use pre-built indices
- Consider Kraken2 for large datasets
