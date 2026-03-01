#!/usr/bin/env bash
# Example: Generate visualization for GeneWise results

# Example 1: Visualize during pipeline run
easypseudogene \
  --proteins human.proteins.fa \
  --genome dolphin.fa \
  --threads 96 \
  --visualize

# Example 2: Standalone visualization from existing results
easypseudogene-visualize \
  --genewise easypseudogene_output/dolphin/genewise/pseudogene-genewise.txt \
  --output dolphin_pseudogenes.html

# Example 3: Auto-detect paths from results directory
easypseudogene-visualize \
  --species dolphin \
  --results-dir ./easypseudogene_output \
  --output visualization.html

# Example 4: Visualize all results (no filtering)
easypseudogene-visualize \
  --genewise easypseudogene_output/dolphin/genewise/pseudogene-genewise.txt \
  --output all_pseudogenes.html
