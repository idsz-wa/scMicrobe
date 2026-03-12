# Installation Guide

## Prerequisites

- Python 3.8+
- Go 1.21+ (optional, for high-performance BAM parsing)
- minimap2
- samtools
- kraken2 (optional)

## Quick Installation

### 1. Clone or Download

```bash
cd scmicro
```

### 2. Install Python Dependencies

```bash
pip install -r python/requirements.txt
```

### 3. Build Go Components (Optional)

```bash
cd go
go mod tidy
go build -o scmicro-go .
# Or use make:
make build
```

### 4. Install External Tools

#### Using Conda (Recommended)

```bash
conda install -c bioconda minimap2 samtools kraken2
```

#### Using apt (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install minimap2 samtools
```

#### From Source

```bash
# minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make

# samtools
git clone https://github.com/samtools/samtools.git
cd samtools && make
```

## Verify Installation

```bash
# Check Python dependencies
python -c "import pysam, numpy, scipy, pandas, scanpy; print('Python dependencies OK')"

# Check external tools
minimap2 --version
samtools --version
kraken2 --version  # optional

# Check Go component (if built)
./go/scmicro-go -h
```

## Build Microbial Reference Database

```bash
# Download and build database
./scripts/build_db.sh ./microbe_db /path/to/host.fa

# This will create:
# - microbe_db/microbe_combined.fa
# - microbe_db/microbe_combined.mmi (minimap2 index)
# - microbe_db/kraken2_db/ (Kraken2 database, optional)
```

## Test with Example Data

```bash
# Create test data
python testdata/create_test_data.py ./test_output

# Run pipeline on test data
python python/scmicro.py \
    --input test_output/test_sample.bam \
    --output ./test_results \
    --host-ref test_output/host.fa \
    --microbe-ref test_output/microbes/E.coli.fa \
    --threads 2
```

## Docker Installation (Optional)

Create a Dockerfile:

```dockerfile
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    minimap2 \
    samtools \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY python/requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt

# Copy scMicrobe
COPY . /opt/scmicro
WORKDIR /opt/scmicro

ENTRYPOINT ["python", "python/scmicro.py"]
```

Build and run:

```bash
docker build -t scmicro .
docker run -v $(pwd):/data scmicro --input /data/sample.bam --output /data/output
```

## Troubleshooting

### pysam installation fails

```bash
# Install build dependencies first
sudo apt-get install libbz2-dev liblzma-dev
pip install pysam
```

### minimap2 not found

Make sure minimap2 is in your PATH:

```bash
export PATH=$PATH:/path/to/minimap2
```

### Memory issues during database build

Reduce the number of genomes:

```bash
# Edit build_db.sh and change head -n 1000 to a smaller number
head -n 100 bacteria.complete.list > bacteria.subset.list
```

## Uninstallation

```bash
# Remove Python packages
pip uninstall pysam numpy scipy pandas scanpy anndata

# Remove scMicrobe
rm -rf scmicro

# Remove external tools (if installed manually)
rm /usr/local/bin/minimap2
rm /usr/local/bin/samtools
```
