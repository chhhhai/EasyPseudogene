# EasyPseudogene

EasyPseudogene is an automated pipeline for identifying pseudogenes in target genomes by comparing them with reference protein sequences. The pipeline uses mmseqs2, miniprot, and GeneWise to detect genes with premature stop codons and/or frameshifts.

## Features

- **Automated Pipeline**: Single command to run the entire pseudogene identification process
- **Comprehensive Filtering**: Multiple filtering steps to remove false positives
- **Flexible Configuration**: Customizable filtering parameters via config file
- **Whitelist Support**: Preserve known reference pseudogenes through all filtering stages
- **High Recall**: Achieves >90% recall on reference datasets

## Requirements

### Dependencies

- **mmseqs2**: Fast protein sequence searching
- **miniprot**: Spliced protein-to-genome alignment
- **GeneWise**: Gene prediction and disruption detection
- **samtools**: FASTA indexing and extraction
- **Python 3**: All filtering scripts
- **Standard Unix tools**: awk, sort, comm, bash

### Installation

1. Install dependencies:
   ```bash
   # mmseqs2
   conda install -c bioconda mmseqs2
   
   # miniprot
   conda install -c bioconda miniprot
   
   # samtools
   conda install -c bioconda samtools
   
   # GeneWise
   conda install -c bioconda wise2
   ```

2. Install EasyPseudogene:
   ```bash
   cd EasyPseudogene
   chmod +x bin/easypseudogene
   chmod +x scripts/*.sh
   ```

3. Add to PATH (optional):
   ```bash
   export PATH=$PATH:$(pwd)/bin
   ```

## Usage

### Basic Usage

```bash
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --visualize
```

### With Whitelist

```bash
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --whitelist known_pseudogenes.txt \
  --visualize
```

### With Custom Configuration

```bash
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --config config/custom.conf
```

### With Visualization

```bash
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --visualize
```

### Advanced Options

```bash
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --species dolphin \
  --output ./results \
  --wiseconfig /path/to/wisecfg \
  --wisepath /path/to/genewise/bin \
  --whitelist ref_list.txt \
  --visualize
```

## Input Files

### Required

- **--proteins**: Reference protein FASTA file
  - Format: Standard FASTA with protein sequences
  - Headers should contain transcript IDs (e.g., `>ENST00000123456.1`)

- **--genome**: Target genome FASTA file
  - Format: Standard FASTA with genomic sequences
  - Can be multi-chromosome or single sequence

- **--threads**: Number of threads to use
  - Recommended: 32-128 depending on system

### Optional

- **--whitelist**: File with reference transcript IDs (one per line)
  - These transcripts will be preserved through all filtering stages
  - Useful for maintaining high recall on known pseudogenes

## Output

The pipeline creates the following directory structure:

```
output_dir/
└── species_name/
    ├── mmseqs/              # mmseqs2 search results
    ├── miniprot/            # miniprot alignment results
    ├── regions/             # Extracted genomic regions
    ├── genewise/            # GeneWise predictions
    └── post_genewise_filter/ # Filtered results
        ├── putative_pseudogene_ids.txt
        ├── unitary_pseudogene_ids.txt
        ├── high_conf_final_pseudogene_ids.txt
        ├── high_confidence_pseudogene_ids.txt   # Final high-confidence results
```

### Key Output Files

- **high_confidence_pseudogene_ids.txt**: Final high-confidence pseudogenes (recommended for downstream analysis)
- **putative_pseudogene_ids.txt**: All putative pseudogenes before filtering
- **multi_mutation_for_manual_review.txt**: Candidates with multiple disruptions requiring manual review
- **genewise_visualization.html**: Interactive HTML visualization (if --visualize option used)

## Configuration

Edit `config/default.conf` to customize filtering parameters:

```bash
# Disable family filter
ENABLE_FAMILY_FILTER=0

# Adjust global alignment quality thresholds
GLOBAL_MIN_COV=0.8
GLOBAL_MIN_ID=0.4
GLOBAL_MIN_ALN_LEN=250
```

## Pipeline Steps

1. **Pre-GeneWise Screening** (mmseqs2 + miniprot)
   - Fast protein-to-genome mapping
   - Spliced alignment identification
   - Candidate region extraction

2. **GeneWise Analysis**
   - Gene prediction on candidate regions
   - Detection of stop codons and frameshifts

3. **Post-GeneWise Filtering**
   - Prefilter: Keep only transcripts with gene_symbol
   - Unitary filters: Remove large families, predicted genes, redundancy
   - False-positive filters: Remove low-quality alignments
   - Quality filters: Multi-disruption and global alignment quality

## Examples

### Example : Basic Run (Take the comparison between humans and cetaceans as an example)

```bash
# abstract ensembl human longest pep
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz

python ./examples/ensembl_select_longest_cds_then_protein.py --cds ./Homo_sapiens.GRCh38.cds.all.fa --pep ./Homo_sapiens.GRCh38.pep.all.fa  --out-prefix human # human.longest_cds.proteins.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/595/GCF_011762595.2_mTurTru1.mat.Y/GCF_011762595.2_mTurTru1.mat.Y_genomic.fna.gz

mv GCF_011762595.2_mTurTru1.mat.Y_genomic.fna dolphin.fa

easypseudogene \
  --proteins human.longest_cds.proteins.fa \
  --genome dolphin.fa \
  --threads 96 \
  --output result \
  --visualize
```

### 

## Troubleshooting

### GeneWise Not Found

If GeneWise is not in PATH, specify paths:

```bash
easypseudogene \
  --proteins proteins.fa \
  --genome genome.fa \
  --threads 96 \
  --wiseconfig /path/to/wisecfg \
  --wisepath /path/to/genewise/bin
```

### Low Memory

Reduce threads or adjust mmseqs2 parameters:

```bash
THREADS=32 MMSEQS_MAX_SEQS=3 easypseudogene \
  --proteins proteins.fa \
  --genome genome.fa \
  --threads 32
```

### Skip Steps

Skip GeneWise if already run:

```bash
easypseudogene \
  --proteins proteins.fa \
  --genome genome.fa \
  --threads 96 \
  --skip-genewise
```

## Visualization

EasyPseudogene includes an interactive HTML visualization tool to explore GeneWise results.

### Generate Visualization After Pipeline

Add `--visualize` flag to automatically generate visualization:

```bash
easypseudogene \
  --proteins human.proteins.fa \
  --genome dolphin.fa \
  --threads 96 \
  --visualize
```

### Standalone Visualization

Use the standalone visualization tool on existing results:

```bash
easypseudogene-visualize \
  --genewise output/dolphin/genewise/pseudogene-genewise.txt \
  --output visualization.html
```

Or auto-detect from results directory:

```bash
easypseudogene-visualize \
  --species dolphin \
  --results-dir ./easypseudogene_output \
  --output visualization.html
```

### Visualization Features

The HTML report includes:
- **Statistics Dashboard**: Overview of total pseudogenes, mutation types, and scores
- **Interactive Search**: Filter pseudogenes by transcript ID
- **Detailed Alignments**: View query vs target protein alignments with highlighted mutations
- **Mutation Highlights**: Stop codons and frameshifts are color-coded in the alignment
- **Mutation Summary**: List of all mutations per pseudogene

Open the generated HTML file in any web browser to explore the results interactively.

## Citation

If you use EasyPseudogene in your research, please cite:

- mmseqs2: Steinegger & Söding (2017) Nature Biotechnology
- miniprot: Li (2023) Bioinformatics
- GeneWise: Birney et al. (2004) Genome Research

## License

This software is provided as-is for research purposes.

## Support

For issues and questions, please check the documentation or contact the developers.
