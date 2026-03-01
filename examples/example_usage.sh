#!/bin/bash
# Example usage of EasyPseudogene

# Example 1: Basic usage
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96

# Example 2: With whitelist
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --whitelist known_pseudogenes.txt

# Example 3: With custom output directory
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --output ./my_results

# Example 4: With GeneWise configuration
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --wiseconfig /path/to/wisecfg \
  --wisepath /path/to/genewise/bin

# Example 5: Skip steps (if already run)
easypseudogene \
  --proteins reference.proteins.fa \
  --genome target.genome.fa \
  --threads 96 \
  --skip-genewise
