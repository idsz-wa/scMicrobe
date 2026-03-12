# scMicrobe

**Single-Cell Microbial Quantification from Host scRNA-seq BAM**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Go 1.20+](https://img.shields.io/badge/go-1.20+-00ADD8.svg)](https://golang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

scMicrobe is a computational framework for detecting and quantifying microbial signals at single-cell resolution directly from host scRNA-seq alignment BAM files.

## Overview

Single-cell RNA sequencing (scRNA-seq) experiments occasionally capture microbial RNA together with host transcripts. Most current pipelines discard non-host reads, thereby losing potential information about intracellular infection, microbial symbiosis, or environmental contamination.

scMicrobe fills this gap by:
- Extracting unmapped or low-quality reads from BAM files
- Applying two-step host filtering to reduce false positives
- Aligning candidate reads to microbial reference databases
- Building microbe-by-cell abundance matrices
- Correcting for environmental contamination

## Features

### Core Capabilities

- **Dual Implementation**: Python for flexibility, Go for high-performance BAM parsing
- **Two-Step Host Filtering**: MAPQ-based primary filter + sequence identity secondary filter
- **Multiple Alignment Strategies**: minimap2 for accuracy, Kraken2 for speed
- **UMI-Aware Quantification**: Exact and Hamming distance-based deduplication
- **Contamination Correction**: Ambient RNA estimation and cell specificity scoring
- **Standard Output Formats**: Matrix Market (MTX), H5AD, TSV - compatible with Seurat and Scanpy

### Pipeline Steps

1. **BAM Parsing** - Extract cell barcodes, UMIs, and candidate reads
2. **Host Filtering** - Remove host-derived reads (MAPQ ≥ 30, identity > 95%)
3. **Microbial Alignment** - Align to microbial reference (minimap2 or Kraken2)
4. **Quantification** - Count reads/UMIs per cell-species
5. **Contamination Filtering** - Remove ambient RNA and non-specific signals
6. **Output Generation** - Write matrices in standard formats

## Installation

### Prerequisites

- Python 3.8 or higher
- Go 1.20 or higher (optional, for high-performance BAM parsing)
- minimap2
- samtools
- kraken2 (optional)

### Quick Install

```bash
# Clone repository
git clone <repository-url>
cd scmicro

# Install Python dependencies
pip install -r python/requirements.txt

# Build Go components (optional)
cd go
go mod tidy
go build -o scmicro-go .
```

### External Tools

Using conda (recommended):
```bash
conda install -c bioconda minimap2 samtools kraken2
```

Using apt (Ubuntu/Debian):
```bash
sudo apt-get install minimap2 samtools
```

## Usage

### Basic Usage

```bash
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --host-ref host.fa \
    --microbe-ref microbes.fa \
    --threads 8
```

### Using Kraken2 for Fast Classification

```bash
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --aligner kraken2 \
    --kraken-db /path/to/kraken2_db \
    --threads 8
```

### High-Performance BAM Parsing (Go)

```bash
# Extract candidate reads using Go parser
./go/scmicro-go \
    -i sample.bam \
    -o-fastq candidates.fastq.gz \
    --parallel \
    --threads 16

# Then run Python pipeline on extracted reads
python python/scmicro.py \
    --input candidates.fastq.gz \
    --output ./output \
    --host-ref host.fa \
    --microbe-ref microbes.fa
```

### Custom Barcode and UMI Tags

```bash
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --barcode-tag CR \
    --umi-tag UR \
    --host-ref host.fa \
    --microbe-ref microbes.fa
```

### Skip Contamination Filtering

```bash
python python/scmicro.py \
    --input sample.bam \
    --output ./output \
    --host-ref host.fa \
    --microbe-ref microbes.fa \
    --skip-contamination-filter
```

## Command Line Options

### Input/Output

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--input` | `-i` | Input BAM file (required) | - |
| `--output` | `-o` | Output directory | `./scmicro_output` |
| `--config` | `-c` | Configuration YAML file | - |

### Reference Files

| Option | Description |
|--------|-------------|
| `--host-ref` | Host reference genome (FASTA) |
| `--microbe-ref` | Microbial reference genome (FASTA) |
| `--kraken-db` | Kraken2 database path |
| `--minimap2-index` | Pre-built minimap2 index |

### Alignment Options

| Option | Description | Default |
|--------|-------------|---------|
| `--aligner` | Alignment method: `minimap2` or `kraken2` | `minimap2` |

### BAM Tags

| Option | Description | Default |
|--------|-------------|---------|
| `--barcode-tag` | Cell barcode tag in BAM | `CB` |
| `--umi-tag` | UMI tag in BAM | `UB` |

### Filtering Options

| Option | Description | Default |
|--------|-------------|---------|
| `--host-mapq` | Host MAPQ threshold for primary filtering | `30` |
| `--host-identity` | Host sequence identity threshold (%) | `95` |
| `--min-microbe-mapq` | Minimum MAPQ for microbial alignment | `1` |

### Quantification Options

| Option | Description | Default |
|--------|-------------|---------|
| `--quant-method` | Quantification method: `reads` or `umi` | `umi` |
| `--umi-dedup-method` | UMI deduplication: `exact` or `hamming` | `exact` |

### Contamination Correction

| Option | Description | Default |
|--------|-------------|---------|
| `--ambient-cutoff` | Ambient RNA cutoff quantile | `0.1` |
| `--min-cells` | Minimum cells for microbe detection | `3` |
| `--skip-contamination-filter` | Skip contamination filtering | `False` |

### Performance

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--threads` | `-t` | Number of threads | `4` |
| `--chunk-size` | | Reads per chunk for processing | `10000` |

## Output Files

### Matrix Files

| File | Description |
|------|-------------|
| `mtx/matrix.mtx` | Matrix Market format (cells × microbes) |
| `mtx/barcodes.tsv` | Cell barcodes |
| `mtx/features.tsv` | Microbe feature names |
| `scmicro_matrix.tsv` | Dense TSV format (human-readable) |
| `scmicro.h5ad` | H5AD format for Scanpy |

### Statistics and Reports

| File | Description |
|------|-------------|
| `scmicro_stats.json` | Pipeline statistics |
| `scmicro_qc.txt` | QC report with summary metrics |

### Example Output Structure

```
output/
├── mtx/
│   ├── matrix.mtx
│   ├── barcodes.tsv
│   └── features.tsv
├── scmicro_matrix.tsv
├── scmicro.h5ad
├── scmicro_stats.json
└── scmicro_qc.txt
```

## Downstream Analysis

### Loading in Scanpy

```python
import scanpy as sc

# Load H5AD file
adata = sc.read_h5ad('output/scmicro.h5ad')

# Basic analysis
sc.pp.filter_cells(adata, min_counts=1)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Visualization
sc.tl.umap(adata)
sc.pl.umap(adata, color=['species'])
```

### Loading in Seurat (R)

```r
library(Seurat)

# Load MTX files
data <- ReadMtx(
  mtx = "output/mtx/matrix.mtx",
  cells = "output/mtx/barcodes.tsv",
  features = "output/mtx/features.tsv"
)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = data)
```

## Building Microbial Reference Database

```bash
# Build database with host filtering
./scripts/build_db.sh ./microbe_db /path/to/host.fa

# This creates:
# - microbe_db/microbe_combined.fa
# - microbe_db/microbe_combined.mmi (minimap2 index)
# - microbe_db/kraken2_db/ (Kraken2 database, optional)
```

## Testing

```bash
# Create test data
python testdata/create_test_data.py ./test_data

# Run pipeline on test data
python python/scmicro.py \
    --input test_data/test_sample.bam \
    --output ./test_output \
    --host-ref test_data/host.fa \
    --microbe-ref test_data/microbes/E.coli.fa
```

## Architecture

### Python Modules

| Module | Purpose |
|--------|---------|
| `scmicro.py` | Main pipeline orchestration |
| `bam_parser.py` | BAM file parsing and candidate extraction |
| `host_filter.py` | Two-step host filtering |
| `microbe_aligner.py` | Microbial alignment (minimap2/Kraken2) |
| `quantifier.py` | UMI deduplication and quantification |
| `contamination_filter.py` | Contamination correction |
| `output_writer.py` | Output formatting |

### Go Components

| Component | Purpose |
|-----------|---------|
| `bam_parser.go` | High-performance BAM parsing |
| `fastq_utils.go` | FASTQ processing utilities |

See [ARCHITECTURE.md](ARCHITECTURE.md) for detailed system architecture.

## Performance

| Dataset Size | Python Only | With Go Parser |
|--------------|-------------|----------------|
| 1 GB BAM | ~5 min | ~2 min |
| 10 GB BAM | ~45 min | ~15 min |
| 100 GB BAM | ~8 hours | ~2.5 hours |

*Benchmarks on 16-core system with 64GB RAM*

## Troubleshooting

### Low microbial detection rate

- Check host filtering parameters (`--host-mapq`, `--host-identity`)
- Verify microbial reference database completeness
- Consider using Kraken2 for more sensitive detection

### High contamination

- Enable contamination filtering (default)
- Adjust `--ambient-cutoff` based on your data
- Check for reagent contamination in negative controls

### Memory issues

- Reduce `--chunk-size`
- Use Go-based BAM parser for lower memory footprint
- Process data in batches

### BAM index errors

If you see "mapping information not recorded in index":
- The pipeline will still work, but progress bar won't show total
- To fix, index your BAM: `samtools index sample.bam`

## Citation

If you use scMicrobe in your research, please cite:

```
scMicrobe: Haven't decided on a name yet
[Authors]
[Journal]
[Year]
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- Built with [pysam](https://github.com/pysam-developers/pysam) for BAM manipulation
- Uses [biogo/hts](https://github.com/biogo/hts) for high-performance Go BAM parsing
- Inspired by single-cell microbial detection methods in the field

## Contact

For questions or issues, please open an issue on the project repository.

