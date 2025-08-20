# HumanEmbryo_scRNA

This repository includes codes used in the bioinformatic analysis of scRNA-seq data for the paper  
*"Single-cell transcriptomics of post-implantation human embryos reveal lineage dynamics and ploidy-dependent signatures"*.

- **cell_demultiplexing.txt**  
  Codes for processing of FASTQ files via CellRanger and single cell demultiplexing via Freemuxlet.

- **embryo_scRNA_data_processing**  
  Codes for scRNA-seq data QC and filtering, integration, lineage-specific markers identification, and copy number variations (CNV) calling for aneuploid and euploid cell classification.

- **gene_expression_analysis**  
  Codes for pseudobulk generation and differential expression analysis to identify enriched markers, as well as gene set enrichment analysis (GSEA).

- **interactome_analysis.R**  
  Codes for analysis of ligand–receptor interactions and crosstalk strength between aneuploid and euploid populations of mosaic embryos.

- **01_monocle_lineage_trajectories.R**  
  Codes for lineage-resolved trajectory inference using Monocle3. Includes conversion of Seurat objects into Monocle3 CellDataSets, clustering, graph learning, pseudotime ordering, identification of trajectory-dependent genes, and marker validation for Epiblast (EPI), Hypoblast (HYPO), and Trophoblast (TB) lineages.

- **02_scenic_core.R**  
  Codes for gene regulatory network inference using SCENIC. Includes gene filtering, co-expression network construction (GENIE3), motif enrichment to build regulons (RcisTarget), AUCell-based regulon activity scoring, and calculation of Regulon Specificity Scores (RSS) across embryo annotations.

