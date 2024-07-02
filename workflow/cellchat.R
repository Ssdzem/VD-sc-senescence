######### Instructions ###########
######## Variables ###############

input.seurat = "/datos/sensence/emilio/VD/output/seurat/aem_GRN_pseudobulk.rds"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"

sample.name = "aem_cellchat"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat"
#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/aem/cellchat"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"

######## Libraries ##############

library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)

####### Input prep ############

seurat_obj <- readRDS(input.seurat)
Idents(seurat_obj) <- seurat_obj@meta.data[["SingleR.labels"]]
data.input <- seurat_obj[["RNA3"]]@data # normalized data matrix
labels <- Idents(seurat_obj)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = seurat_obj, group.by = "ident", assay = "RNA3")

######## setup interaction db ###########

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat@DB <- CellChatDB.use

####### Cellchat obj prepros ###########

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)