# HumanEmbryo_scRNA

This repository includes codes used in the bioinformatic analysis of scRNA-seq data for paper "Single-cell transcriptomics of post-implantation human embryos reveal lineage dynamics and ploidy-dependent signatures".

"cell_demultiplexing.txt" contains codes for processing of fastq files via CellRanger and single cell demultiplexing via Freemuxlet.

"embryo_scRNA_data_processing" contains codes for scRNA-seq data QC and filtering, integration, lineage-specific markers identification, and copy number variations (CNV) calling for aneuploid and euploid cell classification.

"gene_expression_analysis" contains codes for Pseudobulk, differential expression analysis for identification of enriched markers, gene set enrichment analysis (GSEA).

"interactome_analysis.R" contains codes for analysis of ligand-receptor interactions and crosstalk strength between aneuploid and euploid populations of the mosaic embryos.
