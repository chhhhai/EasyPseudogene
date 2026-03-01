# EasyPseudogene Quick Start

## Installation

```bash
# Install dependencies
conda create -n easypseudogene python=3.9
conda activate easypseudogene
conda install -c bioconda mmseqs2 miniprot samtools genewise

# Install EasyPseudogene
cd EasyPseudogene
chmod +x bin/easypseudogene scripts/*.sh
export PATH=$PATH:$(pwd)/bin
```

## Basic Usage

```bash
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96
```

## Output

Results will be in `./easypseudogene_output/species_name/post_genewise_filter/`:

- `high_confidence_pseudogene_ids.txt` - Final high-confidence pseudogenes (recommended)
- `putative_pseudogene_ids.txt` - All putative pseudogenes

## More Examples

See `README.md` for detailed documentation and examples.
