# Differential expression (DE) analysis between cells from four comparison groups

# 1. Load requried packages for DE analysis
library(EnhancedVolcano)
library(EdgeR)

# 2. Run pseudobulk and DE analysis using EdgeR
balance_group$barcode <- rownames(balance_group)
balance_group$samples <- paste0(balance_group$emb_type, balance_group$emb_barcode
aggr <- AggregateExpression(balance_group, group.by = "samples",
                            assays = "RNA",
                            slot = "counts",
                            return.seurat = F)
aggr <- aggr$RNA
aggr <- as.data.frame(aggr)
anno <- data.frame(samples = colnames(aggr))
anno <- anno %>% mutate(condition = ifelse(grepl("Aneuploid", samples), "Aneuploid", "Euploid"))
anno <- column_to_rownames(anno, var = "samples")
conditions <- factor(anno$condition)
conditions <- relevel(conditions, ref="Euploid")
model.desig <- model.matrix(~0 + conditions)
colnames(model.design) <- make.names(gsub("conditions", "", colnames(model.design)))
rownames(model.design) <- rownames(anno)
contrast <- makeContrasts("Aneuploid-Euploid", levels = colnames(model.design))

y <- DGEList(aggr)
y <- calcNormFactors(y)
y <- estimateDisp(y, design = model.design)
fit <- glmQLFit(y, model.design)
y_LRT <- glmQLFTest(fit, contrast = contrast)
edger1 <- as.data.frame(topTags(y_LRT, 1000000))
edger1$gene <- rownames(edger1)
edger1$rank <- sign(edger1$logFC)*-log10(edger1$FDR) # For GSEA
edger1$minuslogFDR <- -log10(edger1$FDR)

# 3. Volcano plot showing DEGs
EnhancedVolcano(edger1, 
                rownames(edger1),
                x ="logFC", 
                y ="FDR",
                FCcutoff = 0,
                pCutoff = 0.05,
                labSize = 4,
                gridlines.major = F,
                gridlines.minor = F,
                axisLabSize = 13,
                labFace = 'bold',
                boxedLabels = TRUE,
                xlim = c(-4,4),
                drawConnectors = TRUE,
                max.overlaps = 200


# Find common DEGs between comparison groups that contain aneuploid cells
# 1. Subset upregulated genes from each comparison group
UP_AE <- subset(edger1, subset = FDR < 0.05 & logFC > 1)
UP_AEmos <- subset(edger2, subset = FDR < 0.05 & logFC > 1)
UP_AA <- subset(edger3, subset = FDR < 0.05 & logFC > 1)

# 2. Subset downregulated genes from each comparison group
DN_AE <- subset(edger1, subset = FDR < 0.05 & logFC < -1)
DN_AEmos <- subset(edger2, subset = FDR < 0.05 & logFC < -1)
DN_AA <- subset(edger3, subset = FDR < 0.05 & logFC < -1)

# 3. Find common DEGs in all 3 groups
all3_up <- UP_AE[rownames(UP_AE) %in% rownames(UP_AEmos) & rownames(UP_AE) %in% rownames(UP_AA),]
all3_dn <- DN_AE[rownames(DN_AE) %in% rownames(DN_AEmos) & rownames(DN_AE) %in% rownames(DN_AA),]

# 4. Find common DEGs between A vs. E and Amosaic vs. Emosaic
AEAEmos_UP <- UP_AE[rownames(UP_AE) %in% rownames(UP_AEmos),]
AEAEmos_DN <- DN_AE[rownames(DN_AE) %in% rownames(DN_AEmos),]

# 5. Find common DEGs between A vs. E and Amosaic vs. A
AEAA_UP <- UP_AE[rownames(UP_AE) %in% rownames(UP_AA),]
AEAA_DN <- DN_AE[rownames(DN_AE) %in% rownames(DN_AA),]

# 6. Find common DEGs between Amosaic vs. Emosaic and Amosaic vs. A
AEmosAA_UP <- UP_AEmos[rownames(UP_AEmos) %in% rownames(UP_AA),]
AEmosAA_DN <- DN_AEmos[rownames(DN_AEmos) %in% rownames(DN_AA),]


# Run GSEA on DEGs of each comparison group
# Load requried packages
library(fgsea)
library(AnnotationDbi)
pathways <- fgsea::gmtPathways("../Human_GOBP_AllPathways_no_GO_iea_August_08_2023_symbol.gmt")

# gsea for DEGs in A vs. E
rank1 <- edger1$rank
names(rank1) <- rownames(edger1)
gsea1 <- fgsea(pathways = pathways,
               stats = rank1,
               minSize = 15,
               maxSize = 500)

# gsea for DEGs in A mosaic vs. E mosaic
rank2 <- edger2$rank
names(rank2) <- rownames(edger2)
gsea2 <- fgsea(pathways = pathways,
               stats = rank2,
               minSize = 15,
               maxSize = 500)

# gsea for DEGs in A mosaic vs. A
rank3 <- edger3$rank
names(rank3) <- rownames(edger3)
gsea3 <- fgsea(pathways = pathways,
               stats = rank3,
               minSize = 15,
               maxSize = 500)

# gsea for DEGs in E mosaic vs. E
rank4 <- edger4$rank
names(rank4) <- rownames(edger4)
gsea4 <- fgsea(pathways = pathways,
               stats = rank4,
               minSize = 15,
               maxSize = 500)

# Enrichment plot for selected pathways in A mosaic vs. E mosaic (rank2)
plotEnrichment(pathways[["HALLMARK_DNA_REPAIR%MSIGDBHALLMARK%HALLMARK"]], 
               rank2, ticksSize = 0.2)
plotEnrichment(pathways[["HALLMARK_GLYCOLYSIS%MSIGDBHALLMARK%HALLMARK"]], 
               rank2, ticksSize = 0.2)
plotEnrichment(pathways[["UBIQUITIN-DEPENDENT PROTEIN CATABOLIC PROCESS%GOBP%GO:0006511"]], 
               rank2, ticksSize = 0.2)
plotEnrichment(pathways[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE%MSIGDBHALLMARK%HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]], 
               rank2, ticksSize = 0.2)







