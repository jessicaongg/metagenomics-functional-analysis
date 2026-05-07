# Methodology

This repository demonstrates a general metagenomics and microbiome functional profiling workflow using R.

## Analysis Workflow

1. **Input data preparation**
   - Taxonomic count tables and sample metadata are loaded into R.
   - Sample IDs are matched between feature tables and metadata.

2. **Alpha diversity analysis**
   - Shannon, Observed, Simpson and Chao diversity metrics are calculated using `phyloseq`.
   - Diversity metrics are summarised by sample group.

3. **Taxonomic profiling**
   - Count data are transformed into relative abundance.
   - Dominant taxa and top taxa by group are identified.

4. **Shared and unique taxa analysis**
   - Taxa are classified as shared, unique or partially shared across groups.

5. **Functional profiling**
   - KEGG Ortholog counts are summarised by group.
   - KEGG pathways are mapped using `KEGGREST`.
   - Functional pathways are classified as core, unique or partially shared.

6. **Visualisation**
   - Diversity and taxonomic composition plots are generated using `ggplot2`.

## Notes

Raw sequencing data and large intermediate files are not included in this repository.
Example input tables should be placed in the `data/` directory if running the workflow locally.
