# Embryo scRNA-seq data processing from cellranger output
# Identify lineage-specific and ploidy-specific expression patterns 

# Load requried packages
library(tidyverse)
library(dplyr)
library(Seurat)
library(harmony)
library(infercnv)

# Create seurat objects with cellranger count output from all 25 embryos 
lib_all <- Read10X(data.dir = "../filtered_feature_bc_matrix")
seurat_all <- CreateSeuratObject(counts = lib_all, 
                                    min.cells = 3, 
                                    min.features = 200)

# Add freemuxlet cell demultiplexing results
demuxafy_all <- read.table("../combined_results.tsv")
seurat_all <- AddMetaData(seurat_all, demuxafy_all)
seurat_all <- subset(seurat_all, 
                     subset = Freemuxlet_DropletType == "singlet") # Remove doublets

# Remove cells with mitochondrial gene content >25% or UMI count <1000
lib_all[["percent.mito"]] <- PercentageFeatureSet(lib_all, pattern = "^MT-")
lib_all <- subset(lib_all, subset = percent.mito < 25 & nCount_RNA >= 1000)

# Keep embryos with >10 cells after filtering
lib_all <- subset(lib_all, subset == id %in% 
           c("0","1","2","3","4","7","8","9","10","11","12","13",
             "14","15","16","18","19","20"))

# Run the standard workflow for visualization and clustering
lib_all <- NormalizeData(lib_all)
lib_all <- FindVariableFeatures(lib_all,selection.method="vst", nfeatures = 2000)
lib_all <- ScaleData(lib_all)
lib_all <- RunPCA(lib_all)
ElbowPlot(lib_all, ndims=100)

lib_all <- FindNeighbors(lib_all, reduction = "pca", dims = 1:20)
lib_all <- FindClusters(lib_all, resolution = 0.1)
lib_all <- RunUMAP(lib_all, dims = 1:20)

# Batch correction with Harmony
seurat_all.harmony <- seurat_all %>% 
  RunHarmony(group.by.vars = "batch", plot_convergence = FALSE)
seurat_all.harmony <- seurat_all.harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.1)

cluster_markers <- FindAllMarkers(lib_all, 
                                  logfc.threshold = 1) # Identify cluster-specific markers

# Create UMAP figure after Harmonry integration
DimPlot(seurat_all.harmony, reduction = "umap", label = T)

# Create UMAP labeled by embryo origins
DimPlot(seurat_all.harmony, reduction = "umap", group.by = "id")

# Create heatmap for lineage-specific markers
top_genes <- c("KRT19","PEG10","LGALS3","TPM1","FABP5","GCSH", # CTB markers
               "CGB2","CGB3","CGB5","CGB7","ATF3","ERVW-1","HOPX","ANXA1", # STB markers
               "MMP2","EPAS1","ITGA1", # EVT markers
               "POU5F1","DPPA5","TDGF1","MT1X","MT1G","MT1H","KHDC3L", # EPI markers
               "APOA1","APOA2","APOE","APOC1","MARCKS","S100A4","FN1","COL4A1" # HYPO markers
)

DoHeatmap(seurat_all.harmony, group.by = "seurat_clusters", features = top_genes, label = F)+
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                       mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Violin plot of lineage-specific markers 
CTB.vln <- VlnPlot(seurat_all.harmony, features = c("KRT19","PEG10","LGALS3"), group.by = "seurat_clusters")
STB.vln <- VlnPlot(seurat_all.harmony, features = c("CGB3","ERVW-1","HOPX"), group.by = "seurat_clusters")
EVT.vln <- VlnPlot(seurat_all.harmony, features = c("MMP2","EPAS1","ITGA1"), group.by = "seurat_clusters")
EPI.vln <- VlnPlot(seurat_all.harmony, features = c("POU5F1","DPPA5","ITGA1"), group.by = "seurat_clusters")
HYPO.vln <- VlnPlot(seurat_all.harmony, features = c("APOA1", "MARCKS", "FN1", "S100A4"), group.by = "seurat_clusters")

# Unsupervised inferCNV run with no euploid reference
raw_matrix1 <- seurat_all.harmony@assays$RNA@counts
anno1 <- seurat_all.harmony@meta.data %>% select(id)
anno1$id <- factor(anno1$id)

infercnv_1 <- CreateInfercnvObject(raw_counts_matrix = raw_matrix1,
                                 annotations_file = anno1,
                                 gene_order_file = "../gencode_v19_gene_pos_infercnv.txt",
                                 ref_group_names = NULL, 
                                 chr_exclude = c("chrM"))

infercnv_1 = infercnv::run(infercnv_1,
                         cutoff=0.1,
                         out_dir=paste0("../infercnv_results1"),
                         cluster_by_groups = T,
                         cluster_references = T,
                         denoise = T,
                         HMM = F,
                         num_threads=4,
                         output_format = "png",
                         useRaster = T,
                         per_chr_hmm_subclusters = F)

# Calculate aberration sum for each gene in each cell
cnv_raw <- readRDS("../infercnv_results1/run.final.infercnv_obj")
cnv_raw <- cnv_raw@expr.data %>% as.matrix()
cnv_raw <- apply(cnv_raw, c(1,2), function(x) (x-1)) # Calculate gene-level deviation from global average (1)
cnv_raw <- abs(cnv_raw) 
cnv_raw <- as.data.frame(colSums(cnv_raw)) # Sum up aberration score across all genes for each cell
colnames(cnv_raw)= "cnv_raw"
seurat_all.harmony <- AddMetaData(seurat_all.harmony, cnv_raw) 

# Find cells with lowest aberration sum in each lineage
ab_sum <- seurat_all.harmony@meta.data
ab_sum$lineage_rank <- paste0(ab_sum$seurat_clusters, ab_sum$cnv_raw) 

eup_ref_list <- c("lib230524_GGTTCTCCAAGATGTA-1",
                  "lib230331_CAGGCCAAGCGTGAAC-1",
                  "lib230331_CTCCAACGTTCGGACC-1",
                  "lib230517_AAACCCAAGGTTGACG-1",
                  "lib230517_GAGGGATCAGCTTCGG-1",
                  "lib230517_ACACAGTGTAAGCAAT-1",
                  "lib230517_GCAGCCAGTCACCGCA-1",
                  "lib230517_CCGATGGCATGACGTT-1",
                  "lib230517_TGGGCTGTCGGCTGGT-1",
                  "lib230517_TTGCTGCGTACGAGTG-1") # Create reference list

# Create distribution plot for aberration sum across all cells
VlnPlot(seurat_all.harmony, features = "cnv_raw", cols = "gray")
VlnPlot(seurat_all.harmony, features = "cnv_raw", group.by = "seurat_clusters")

# Label euploid reference cells on UMAP
DimPlot(seurat_all.harmony, reduction = "umap", label = T, label.size = 3, 
        cells.highlight = list(eup_ref_list),sizes.highlight = 0.6)+
  scale_color_manual(labels = c("Unselected", "Euploid reference"), 
                     values = c("#D3D3D3","red"))

# Add reference label to cells for next inferCNV run 
x <- c(rep(21,10)) # set euploid reference cells as cluster 21
eup_cells <- data.frame(eup_ref_list,x)
colnames(eup_cells) <- c("rownames", "infercnv_label")
eup_cells <- column_to_rownames(eup_cells, var = "rownames")

metadata <- seurat_all.harmony@meta.data
metadata <- select(metadata, id)
colnames(metadata) <- "infercnv_label"

metadata <- metadata[!(row.names(metadata) %in% eup_ref_list), ]
metadata <- rbind(metadata, eup_cells)
seurat_all.harmony <- AddMetaData(seurat_all.harmony, metadata)

# Set up parameters for new inferCNV run with euploid reference cells
raw_matrix2 <- seurat_all.harmony@assays$RNA@counts
anno2 <- seurat_all.harmony@meta.data %>% select(infercnv_label)
anno2$infercnv_label <- factor(anno2$infercnv_label)

infercnv2 <- CreateInfercnvObject(raw_counts_matrix = raw_matrix2,
                                 annotations_file = anno2,
                                 gene_order_file = "../gencode_v19_gene_pos_infercnv.txt",
                                 ref_group_names = c("21"), # Annotation for euploid reference cells 
                                 chr_exclude = c("chrM"))

infercnv2 = infercnv::run(infercnv,
                         cutoff=0.1,
                         out_dir=paste0("../infercnv_results2"),
                         cluster_by_groups = T,
                         cluster_references = T,
                         denoise = T,
                         HMM = T,
                         num_threads=4,
                         output_format = "png",
                         useRaster = T,
                         per_chr_hmm_subclusters = F)


# Extract CNV prediction information
cnv = readRDS("../run.final.infercnv_obj")
expr_comb <- cnv@expr.data %>% as.matrix()
expr_comb <- apply(expr_comb, c(1,2), function(x) (x-3)) # Quantifying CNV in each cell from HMM output

# Find gene range for each chromosome
chr <- as.data.frame(cnv@gene_order)

# Create master list of CNV proportion for each chromosome in each cell
# 1. Create an empty list to store the dataframes for each chromosome (excluding chr19, X, Y)
chr_list <- list()

# 2. Loop through chromosome numbers 
for (chr_num in c(1:18, 20:22)) {
  chr_name <- paste0("chr",chr_num)  
  chr_subset <- chr[chr$chr == chr_name, ]  
  chr_list[[chr_name]] <- chr_subset  
}


# Identify gain proportion, loss proportion, and CNV proportion (gain+loss) for each chromosome
# 1. Specify chr numbers to be included
chr_nums <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
              "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
              "chr18", "chr20", "chr21", "chr22")

# 2. Create an empty list to store the results
cnv_expr <- list()

for (chr_num in chr_nums) {
  gene_names <- rownames(chr_list[[as.character(chr_num)]])
  cnv_subset <- expr_comb[rownames(expr_comb) %in% gene_names, ]
  
  if (is.null(dim(cnv_subset)) || any(dim(cnv_subset) < 2)) {
    cnv_result <- data.frame(
      "loss" = NA,
      "gain" = NA
    )
  } else {
    # Create a DataFrame with loss and gain columns
    column_names <- c(paste0("loss", chr_num), paste0("gain", chr_num))
    cnv_result <- data.frame(
      "loss" = colSums(cnv_subset < 0) / nrow(cnv_subset),
      "gain" = colSums(cnv_subset > 0) / nrow(cnv_subset)
    )
  }
  
  # Set the column names
  colnames(cnv_result) <- column_names
  
  # Store the result in the cnv_expr list
  cnv_expr[[as.character(chr_num)]] <- cnv_result
}

# 3. Combine columns in each list within cnv_expr
combined_cnv_expr <- bind_cols(cnv_expr)

# 4. Calculate the sum column for each chromosome
for (chr_num in chr_nums) {
  combined_cnv_expr[[paste0("sum", chr_num)]] <- 
    combined_cnv_expr[[paste0("gain", chr_num)]] +
    combined_cnv_expr[[paste0("loss", chr_num)]]
}


# Find the chromosome with highest CNV proportion in each cell
cnv_summary <- combined_cnv_expr[,-c(1:42)]
cnv_summary$max_cnv_prop <- apply(cnv_summary[1:21], 1, max, na.rm=TRUE)
cnv_summary$max_chr <- names(cnv_summary)[max.col(cnv_summary[1:21], "first")]

# Set threshold to identify aneuploid population 
cnv_summary <- cnv_summary %>% mutate(cnv_status = if_else(
  max_cnv_prop <= 0.7, "Euploid", "Aneuploid"))
cnv_summary <- select(cnv_summary, cnv_status, max_cnv_prop, max_chr)
seurat_all.harmony <- AddMetaData(seurat_all.harmony, cnv_summary)

# Identify embryo types by % Aneuploidy
seurat_all.harmony$samples <- paste0(seurat_all.harmony$cnv_status, seurat_all.harmony$id)
table(seurat_all.harmony$samples) # Find #aneuploid cells and euploid cells in each embryo

meta <- seurat_all.harmony@meta.data
meta$emb_type <- ifelse(meta$id %in% c(0, 4, 12), "Mosaic",
                       ifelse(lin$id %in% c(1, 5, 6, 7, 8, 10), "Aneuploid",
                              ifelse(lin$id %in% c(2, 3, 9, 11, 13, 14, 15, 16, 18, 19, 20), "Euploid", NA)))
meta <- select(meta,emb_type)
seurat_all.harmony <- AddMetaData(seurat_all.harmony,meta)

# Reorganize embryo order based on embryo types
paper_id <- c("12","15","1","2","13","16","17","18","19","3","20","4",
              "14","5","6","7","8","9","10","11")
meta <- seurat_all.harmony@meta.data
meta$paper_id <- paste0(meta$id)
meta$paper_id <- factor(meta$id)
levels(meta$paper_id) <- paper_id
meta <- select(meta, paper_id)
seurat_all.harmony <- AddMetaData(seurat_all.harmony, meta)

# Create distribution plot of max chromosome-CNV proportion for all cells 
VlnPlot(seurat_all.harmony, features = "max_cnv_prop")
VlnPlot(seurat_all.harmony, features = "max_cnv_prop", group.by = "emb_type")

# Label UMAP with cnv_status and emb_type
DimPlot(seurat_all.harmony, reduction = "umap", group.by = "cnv_status", cols = c("red","green"))
DimPlot(seurat_all.harmony, reduction = "umap", group.by = "emb_type")

# Label UMAP with cnv_status in mosaic embryos
mos_eup <- subset(seurat_all.harmony, subset = emb_type == "Mosaic" 
                  & cnv_status == "Euploid") %>% WhichCells()
mos_anp <- subset(seurat_all.harmony, subset = emb_type == "Mosaic" 
                  & cnv_status == "Aneuploid") %>% WhichCells()

mosaic.plot <- DimPlot(seurat_all.harmony, reduction = "umap",
                        cells.highlight = list(mos_anp,mos_eup), sizes.highlight = 0.3)+ 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8))+
  scale_color_manual(labels = c("Other embryos","Euploid cells","Aneuploid cells"), 
                     values = alpha(c("gray","#0CB702","#F8766D"),0.7))

LabelClusters(mosaic.plot, id = "ident", size = 3.5, repel = T, box.padding = 1)












