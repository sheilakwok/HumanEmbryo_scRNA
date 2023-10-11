# Differential expression (DE) analysis between cells from four comparison groups
# 1) A vs. E; 2) Amosaic vs. Emosaic; 3) Amosaic vs. A; 4) Emosaic vs. E

# Load requried packages for DE analysis
library(EnhancedVolcano)
library(EdgeR)

# Subset cells from each embryo type
aneuploid <- subset(seurat_all.harmony, subset = emb_type == "Aneuploid")
euploid <- subset(seurat_all.harmony, subset = emb_type == "Euploid")
mos_anp <- subset(seurat_all.harmony, subset = emb_type == "Mosaic" & cnv_status == "Aneuploid")
mos_eup <- subset(seurat_all.harmony, subset = emb_type == "Mosaic" & cnv_status == "Euploid")

# Find DEGs between A vs. E 
# 1. Balance lineage composition between comparison groups
meta <- aneuploid@meta.data 
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(627, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(27, replace = FALSE)
sampled_HYPO <- meta %>% filter(seurat_clusters == "HYPO") %>% sample_n(11, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(182, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(311, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_HYPO, sampled_EVT, sampled_STB), rownames))
balance1 <- subset(aneuploid, cells = list)

meta <- euploid@meta.data 
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(627, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(27, replace = FALSE)
sampled_HYPO <- meta %>% filter(seurat_clusters == "HYPO") %>% sample_n(11, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(182, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(311, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_HYPO, sampled_EVT, sampled_STB), rownames))
balance2 <- subset(euploid, cells = list)
balance_group1 <- merge(balance1, y=balance2)

# 2. Run pseudobulk and DE analysis using EdgeR
balance_group1$barcode <- rownames(balance_group1)
balance_group1$samples <- paste0(balance_group1$emb_type, balance_group1$barcode
aggr <- AggregateExpression(balance_group1, group.by = "samples",
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
                max.overlaps = 200,
                selectLab = c('CCND1',"DDIT4","HIF1A","RPS28","MCM5",
                              "ATP5MEP2","NLRP7","ATF6")) 


# Find DEGs between A mosaic vs. E mosaic
# 1. Balance lineage composition between comparison groups
meta <- mos_eup@meta.data 
meta <- select(meta, seurat_clusters)
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(29, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(3, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(21, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(121, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_EVT, sampled_STB), rownames))
balance1 <- subset(mos_eup, cells = list)

meta <- mos_anp@meta.data 
meta <- select(meta, seurat_clusters)
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(29, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(3, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(21, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(121, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_EVT, sampled_STB), rownames))
balance2 <- subset(mos_anp, cells = list)
balance_group2 <- merge(balance1, y=balance2)

# 2. Run pseudobulk and DE analysis using EdgeR
balance_group2$barcode <- rownames(balance_group2)
balance_group2$samples <- paste0(balance_group2$cnv_status, balance_group2$barcode
aggr <- AggregateExpression(balance_group2, group.by = "samples",
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
model.design <- model.matrix(~0 + conditions)
colnames(model.design) <- make.names(gsub("conditions", "", colnames(model.design)))
rownames(model.design) <- rownames(anno)
contrast <- makeContrasts("Aneuploid-Euploid", levels = colnames(model.design))

y <- DGEList(aggr)
y <- calcNormFactors(y)
y <- estimateDisp(y, design = model.design)
fit <- glmQLFit(y, model.design)
y_LRT <- glmQLFTest(fit, contrast = contrast)
edger2 <- as.data.frame(topTags(y_LRT, 1000000))
edger2$gene <- rownames(edger2)
edger2$rank <- sign(edger2$logFC)*-log10(edger2$FDR) # For GSEA
edger2$minuslogFDR <- -log10(edger2$FDR)

# 3. Volcano plot showing DEGs
EnhancedVolcano(edger2, 
                rownames(edger2),
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
                ylim = c(0,40),
                xlim = c(-4,4),
                drawConnectors = TRUE,
                max.overlaps = 200,
                selectLab = c("DUT","PHB2","DDIT4","ATF4","UBE2J1","MYC","TUBA1C",
                              "TP63","CCNE1","ULK1","CLU")) 


# Find DEGs between A mosaic vs. A (100% CNV)
# 1. Balance lineage composition between comparison groups
meta <- mos_anp@meta.data 
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(29, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(3, replace = FALSE)
sampled_HYPO <- meta %>% filter(seurat_clusters == "HYPO") %>% sample_n(3, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(21, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(85, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_EVT, sampled_STB, sampled_HYPO), rownames))
balance1 <- subset(mos_anp, cells = list)

aneuploid_100 <- subset(aneuploid, subset = max_cnv_prop == 100) # Subset aneuploid cells with 100% CNV prop.
meta <- aneuploid_100@meta.data 
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(29, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(3, replace = FALSE)
sampled_HYPO <- meta %>% filter(seurat_clusters == "HYPO") %>% sample_n(3, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(21, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(85, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_EVT, sampled_STB,sampled_HYPO), rownames))
balance2 <- subset(aneuploid_100, cells = list)
balance_group3 <- merge(balance1, y=balance2)

# 2. Run pseudobulk and DE analysis using EdgeR
balance_group3$barcode <- rownames(balance_group3)
balance_group3$samples <- paste0(balance_group3$emb_type, balance_group3$barcode
aggr <- AggregateExpression(balance_group3, group.by = "samples",
                          assays = "RNA",
                            slot = "counts",
                            return.seurat = F)
aggr <- aggr$RNA
aggr <- as.data.frame(aggr)
anno <- data.frame(samples = colnames(aggr))
anno <- anno %>% mutate(condition = ifelse(grepl("Mosaic", samples), "Mosaic", "Aneuploid"))
anno <- column_to_rownames(anno, var = "samples")
conditions <- factor(anno$condition)
conditions <- relevel(conditions, ref="Aneuploid")
model.desig <- model.matrix(~0 + conditions)
colnames(model.design) <- make.names(gsub("conditions", "", colnames(model.design)))
rownames(model.design) <- rownames(anno)
contrast <- makeContrasts("Mosaic-Aneuploid", levels = colnames(model.design))

y <- DGEList(aggr)
y <- calcNormFactors(y)
y <- estimateDisp(y, design = model.design)
fit <- glmQLFit(y, model.design)
y_LRT <- glmQLFTest(fit, contrast = contrast)
edger3 <- as.data.frame(topTags(y_LRT, 1000000))
edger3$gene <- rownames(edger3)
edger3$rank <- sign(edger3$logFC)*-log10(edger3$FDR) # For GSEA
edger3$minuslogFDR <- -log10(edger3$FDR)

# 3. Volcano plot showing DEGs
EnhancedVolcano(edger3, 
                rownames(edger3),
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
                max.overlaps = 200,
                selectLab = c("CASP4","DDIT4","HIF1A","CCND1","TGM2","SEMA4D",
                              "CHMP2A","HSP90B1","ATF4","TP63")) 


# Find DEGs between E mosaic vs. E (0% CNV)
# 1. Balance lineage composition between comparison groups
meta <- mos_anp@meta.data 
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(115, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(12, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(69, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(102, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_EVT, sampled_STB), rownames))
balance1 <- subset(mos_anp, cells = list)

euploid_0 <- subset(euploid, subset = max_cnv_prop == 0) # Subset euploid cells with 0% CNV prop.
meta <- euploid_0@meta.data 
sampled_CTB <- meta %>% filter(seurat_clusters == "CTB") %>% sample_n(115, replace = FALSE)
sampled_EPI <- meta %>% filter(seurat_clusters == "EPI") %>% sample_n(12, replace = FALSE)
sampled_EVT <- meta %>% filter(seurat_clusters == "EVT") %>% sample_n(69, replace = FALSE)
sampled_STB <- meta %>% filter(seurat_clusters == "STB") %>% sample_n(102, replace = FALSE)
list <- unlist(lapply(list(sampled_CTB, sampled_EPI, sampled_EVT, sampled_STB), rownames))
balance2 <- subset(euploid_0, cells = list)
balance_group4 <- merge(balance1, y=balance2)

# 2. Run pseudobulk and DE analysis using EdgeR
balance_group4$barcode <- rownames(balance_group4)
balance_group4$samples <- paste0(balance_group4$emb_type, balance_group4$barcode
aggr <- AggregateExpression(balance_group4, group.by = "samples",
                          assays = "RNA",
                            slot = "counts",
                            return.seurat = F)
aggr <- aggr$RNA
aggr <- as.data.frame(aggr)
anno <- data.frame(samples = colnames(aggr))
anno <- anno %>% mutate(condition = ifelse(grepl("Mosaic", samples), "Mosaic", "Euploid"))
anno <- column_to_rownames(anno, var = "samples")
conditions <- factor(anno$condition)
conditions <- relevel(conditions, ref="Euploid")
model.desig <- model.matrix(~0 + conditions)
colnames(model.design) <- make.names(gsub("conditions", "", colnames(model.design)))
rownames(model.design) <- rownames(anno)
contrast <- makeContrasts("Mosaic-Euploid", levels = colnames(model.design))

y <- DGEList(aggr)
y <- calcNormFactors(y)
y <- estimateDisp(y, design = model.design)
fit <- glmQLFit(y, model.design)
y_LRT <- glmQLFTest(fit, contrast = contrast)
edger4 <- as.data.frame(topTags(y_LRT, 1000000))
edger4$gene <- rownames(edger4)
edger4$rank <- sign(edger4$logFC)*-log10(edger4$FDR) # For GSEA
edger4$minuslogFDR <- -log10(edger4$FDR)

# 3. Volcano plot showing DEGs
EnhancedVolcano(edger4, 
                rownames(edger4),
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
                max.overlaps = 200,
                selectLab = c("CEBPD","RPS20","ULK1","PSMD9","CASP4",
                              "HSPA1A","RPS29","ATP5ME")) 


# Create boxplot of logFC values in the 4 comparison groups
edger1_deg <- subset(edger1, subset = FDR < 0.05 & abs(logFC) > 1)
edger2_deg <- subset(edger2, subset = FDR < 0.05 & abs(logFC) > 1)
edger3_deg <- subset(edger3, subset = FDR < 0.05 & abs(logFC) > 1)
edger4_deg <- subset(edger4, subset = FDR < 0.05 & abs(logFC) > 1)

FC <- cbind(edger1_deg$logFC, edger2_deg$logFC, edger3_deg$logFC, edger4_deg$logFC)
colnames(FC) <- c("A vs. E","A mosaic vs. E mosaic","A mosaic vs. A","E mosaic vs. E")
boxplot(FC, ylab = "logFC",
        col = c("lightblue", "lightgreen", "orange", "lightgrey"),
        space = 5)


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







