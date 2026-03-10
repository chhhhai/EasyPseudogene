# EasyPseudogene Installation Guide

## Quick Start

1. **Install Dependencies** (using conda):
   ```bash
   conda create -n easypseudogene python=3.9
   conda activate easypseudogene
   conda install -c bioconda mmseqs2 miniprot samtools wise2
   ```

2. **Install EasyPseudogene**:
   ```bash
   cd EasyPseudogene
   chmod +x bin/easypseudogene scripts/*.sh
   ```

3. **Add to PATH** (optional):
   ```bash
   export PATH=$PATH:$(pwd)/bin
   ```

4. **Test Installation**:
   ```bash
   easypseudogene --help
   ```

## Detailed Installation

### Step 1: Install Dependencies

#### Option A: Using Conda (Recommended)

```bash
# Create environment
conda create -n easypseudogene python=3.9
conda activate easypseudogene

# Install tools
conda install -c bioconda mmseqs2
conda install -c bioconda miniprot
conda install -c bioconda samtools
conda install -c bioconda wise2
```

#### Option B: Manual Installation

1. **mmseqs2**: Download from https://github.com/soedinglab/MMseqs2
2. **miniprot**: Download from https://github.com/lh3/miniprot
3. **samtools**: Download from http://www.htslib.org/
4. **GeneWise**: Download from https://www.ebi.ac.uk/~birney/wise2/

### Step 2: Install EasyPseudogene

```bash
# Clone or copy EasyPseudogene directory
cd EasyPseudogene

# Make scripts executable
chmod +x bin/easypseudogene
chmod +x scripts/*.sh
chmod +x lib/*.py
```

### Step 3: Verify Installation

```bash
# Check dependencies
easypseudogene --help

# Test with a small example (if you have test data)
easypseudogene \
  --proteins test_proteins.fa \
  --genome test_genome.fa \
  --threads 4 \
  --visualize
```

## System Requirements

- **OS**: Linux or macOS
- **Memory**: 8GB minimum, 32GB+ recommended
- **CPU**: Multi-core recommended (32-128 threads)
- **Disk**: 10GB+ free space (depends on genome size)

## Troubleshooting

### Issue: "command not found: mmseqs"

**Solution**: Install mmseqs2 via conda or add to PATH

### Issue: "GeneWise not found"

**Solution**: 
- Install GeneWise via conda, or
- Set `--wiseconfig` and `--wisepath` options

### Issue: "Permission denied"

**Solution**: 
```bash
chmod +x bin/easypseudogene scripts/*.sh lib/*.py
```

### Issue: "Python module not found"

**Solution**: Ensure Python 3.6+ is installed and in PATH

## Next Steps

After installation, see `README.md` for usage examples.
