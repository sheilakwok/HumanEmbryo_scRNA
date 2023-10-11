
# Load requried packages
library(CellChat)
library(reticulate)
library(patchwork)
library(NMF)
library(ggalluvial)

# Subset cells from mosaic embryos 
mosaic_emb <- subset(seurat_all.harmony, subset = emb_type == "Mosaic")

# Standard workflow of CellChat
cellchat.DB <- CellChatDB.human
data.input <- GetAssayData(mosaic_emb, assay = "RNA", slot = "data")
meta <- mosaic_emb@meta.data
meta <- select(meta, cnv_status)
cellchat.mosaic <- createCellChat(object = data.input, meta = meta, group.by = "cnv_status")

cellchat.mosaic@DB <- cellchat.DB
cellchat.mosaic <- subsetData(cellchat.mosaic)
cellchat.mosaic <- identifyOverExpressedGenes(cellchat.mosaic)
cellchat.mosaic <- identifyOverExpressedInteractions(cellchat.mosaic)
cellchat.mosaic <- projectData(cellchat.mosaic, PPI.human)
cellchat.mosaic <- computeCommunProb(cellchat.mosaic,raw.use = FALSE)
cellchat.mosaic <- filterCommunication(cellchat.mosaic, min.cells = 10)
cellchat.mosaic <- computeCommunProbPathway(cellchat.mosaic)
cellchat.mosaic <- aggregateNet(cellchat.mosaic)
groupSize.aneuploid <- as.numeric(table(cellchat.mosaic@idents))

# Create bubbleplot of L-R interactions between aneuploid and euploid cells
netVisual_bubble(cellchat.mosiac, sources.use = 1, targets.use = 2, # Euploid (source) -> Aneuploid (target)
                 remove.isolate = FALSE) 
netVisual_bubble(cellchat.mosiac, sources.use = 2, targets.use = 1, # Aneuploid (source) -> Euploid (target)
                 remove.isolate = FALSE) 

# Compute the network centrality scores
cellchat.mosaic <- netAnalysis_computeCentrality(cellchat.mosaic, slot.name = "netP")

# Scatter plot to visualize aggregated communication networks for each cell type
netAnalysis_signalingRole_scatter(cellchat.mosaic,
                                  label.size = 4,
                                  dot.size = c(4, 6))

# Heatmap to visualize dominant cell types for each signaling pathway
p1 <- netAnalysis_signalingRole_heatmap(cellchat.mosaic, pattern = "outgoing", height = 10)
p2 <- netAnalysis_signalingRole_heatmap(cellchat.mosaic, pattern = "incoming", height = 10)
p1+p2



 




