######### Instructions ###########
#Manual inputs involve the variable adjustement (line 5), if you want to check specific pathway, check out lines 112-122
######## Variables ###############

input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_hanai_ES.rds"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"

sample.name = "hanai_cellchat"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat"
#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/scall/hanai"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"

######## Libraries ##############

library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

####### Input prep ############


## get the pseudobulk
seurat_obj <- readRDS(input.seurat)

## change the labels of the top 10% enriched senmayo cells to "senescent_cells"
senescent_cells <- which(seurat_obj@meta.data[["SENmayo"]] > quantile(seurat_obj@meta.data[["SENmayo"]], probs = c(0.90)))
seurat_obj$SingleR.labels[senescent_cells] <- "senescent_cells"

## subset for vitamin D and vehicle
seurat_vehicle <- subset(x = seurat_obj, subset = experiment == "vehicle")
seurat_vd <- subset(x = seurat_obj, subset = experiment == "vd")

#remove the  mismatched cell types since cellchat is very picky with it 
common_cell_types <- intersect(unique(seurat_vehicle@meta.data[["SingleR.labels"]]), unique(seurat_vd@meta.data[["SingleR.labels"]]))
seurat_vehicle <- subset(seurat_vehicle, subset = SingleR.labels %in% common_cell_types)
seurat_vd <- subset(seurat_vd, subset = SingleR.labels %in% common_cell_types)

#arrange the cell types in the same order since cellchat is also picky with it
seurat_vehicle@meta.data[["SingleR.labels"]] <- factor(seurat_vehicle@meta.data[["SingleR.labels"]], levels = unique(seurat_vehicle@meta.data[["SingleR.labels"]]))
seurat_vd@meta.data[["SingleR.labels"]] <- factor(seurat_vd@meta.data[["SingleR.labels"]], levels = unique(seurat_vehicle@meta.data[["SingleR.labels"]]))

#get the cellchat objects for vehicle and vd
Idents(seurat_vehicle) <- seurat_vehicle@meta.data[["SingleR.labels"]]
labels <- Idents(seurat_vehicle)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat_vehicle <- createCellChat(object = seurat_vehicle, group.by = "ident", assay = "RNA")

Idents(seurat_vd) <- seurat_vd@meta.data[["SingleR.labels"]]
labels <- Idents(seurat_vd)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat_vd <- createCellChat(object = seurat_vd, group.by = "ident", assay = "RNA")

######## setup interaction db ###########

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

####### Cellchat obj prepros ###########

#process intervention
cellchat_vd <- subsetData(cellchat_vd) # This step is necessary even if using the whole database
future::plan("multicore", workers = 50) # do parallel
cellchat_vd <- identifyOverExpressedGenes(cellchat_vd)
cellchat_vd <- identifyOverExpressedInteractions(cellchat_vd)
# cellchat_vd <- projectData(cellchat_vd, PPI.mouse) use this if shallow sequencing
cellchat_vd <- computeCommunProb(cellchat_vd)
cellchat_vd <- computeCommunProbPathway(cellchat_vd)
df.net_vd <- subsetCommunication(cellchat_vd)
cellchat_vd <- aggregateNet(cellchat_vd)
cellchat_vd <- netAnalysis_computeCentrality(cellchat_vd)
gc()
#process vehicle
cellchat_vehicle <- subsetData(cellchat_vehicle) # This step is necessary even if using the whole database
future::plan("multicore", workers = 50) # do parallel
cellchat_vehicle <- identifyOverExpressedGenes(cellchat_vehicle)
cellchat_vehicle <- identifyOverExpressedInteractions(cellchat_vehicle)
#cellchat_vehicle <- projectData(cellchat_vehicle, PPI.mouse) use this if shallow sequencing (currently function not working)
cellchat_vehicle <- computeCommunProb(cellchat_vehicle)
cellchat_vehicle <- computeCommunProbPathway(cellchat_vehicle)
df.net_vehicle <- subsetCommunication(cellchat_vehicle)
cellchat_vehicle <- aggregateNet(cellchat_vehicle)
cellchat_vehicle <- netAnalysis_computeCentrality(cellchat_vehicle)

#merge it boi
object.list <- list(vehicle = cellchat_vehicle, vd = cellchat_vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

rm(seurat_obj, seurat_vd, seurat_vehicle, meta, df.net_vd, df.net_vehicle, CellChatDB.use, CellChatDB)
####### Plots ##########

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
png(filename = file.path(plots.dir, paste0(sample.name, "_celltype_interactions.png")),
    width = 20, height = 11.75, units = "cm", res=300)
gg1 + gg2
dev.off()
## 2d representation change for each cell
gg1 <- netAnalysis_signalingRole_scatter(cellchat_vehicle) 
gg2 <- netAnalysis_signalingRole_scatter(cellchat_vd)
gg1 + gg2


#### here should be the pathway analysis of the most differentiated cell types so watch out, this process is MANUAL

#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Fibroblast_Fbn1 high")
#gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
#patchwork::wrap_plots(plots = list(gg1,gg2))

### Here we make the functional analysis, first we make the 2d plot
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

#here is the sauce, check wich pathways change more IN THE JOINT MANIFOLD THAT MAKES THE CLUSTERS YOU NEED TO CHECK THE FUCKING PAPER TO UNDERSTAND
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, paste0(sample.name, "_path_functionalan.png")),
    width = 30, height = 18, units = "cm", res=300)
gg1 + gg2
dev.off()
### here it actually checks pathways boi
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
png(filename = file.path(plots.dir, paste0(sample.name, "_path_comm.png")),
    width = 12.5, height = 10, units = "in", res=300)
gg1
dev.off()
### heatmap all (lets see if it runs)

# combining all the identified signaling pathways from different datasets 
pathway.union <- union(cellchat_vehicle@netP$pathways, cellchat_vd@netP$pathways)
########### Overall
#vehicle
centr <- slot(cellchat_vehicle, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat_vehicle@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)
for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]]$outdeg
  incoming[, i] <- centr[[i]]$indeg
}
mat <- t(outgoing + incoming)
mat1 <- mat[rownames(mat) %in% pathway.union, , drop = FALSE]
mat <- matrix(0, nrow = length(pathway.union), ncol = ncol(mat))
idx <- match(rownames(mat1), pathway.union)
mat[idx[!is.na(idx)], ] <- mat1
dimnames(mat) <- list(pathway.union, colnames(mat1))
mat.ori_vehicle <- mat

#vd
centr <- slot(cellchat_vd, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat_vehicle@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)
for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]]$outdeg
  incoming[, i] <- centr[[i]]$indeg
}
mat <- t(outgoing + incoming)
mat1 <- mat[rownames(mat) %in% pathway.union, , drop = FALSE]
mat <- matrix(0, nrow = length(pathway.union), ncol = ncol(mat))
idx <- match(rownames(mat1), pathway.union)
mat[idx[!is.na(idx)], ] <- mat1
dimnames(mat) <- list(pathway.union, colnames(mat1))
mat.ori_vd <- mat

### here i want to see the statistically significant differentiated pathways

p_values <- numeric(nrow(mat.ori_vd))
test_results <- list()

# Perform Wilcoxon signed-rank test for each row pair
for (i in 1:nrow(mat.ori_vd)) {
  vd_row <- mat.ori_vd[i, ]
  vehicle_row <- mat.ori_vehicle[i, ]
  test_result <- wilcox.test(vd_row, vehicle_row, paired = TRUE)
  p_values[i] <- test_result$p.value
  test_results[[i]] <- test_result
}

# Display results
result_df_all <- data.frame(
  Row = rownames(mat.ori_vd),
  P_Value = p_values,
  Significant = ifelse(p_values < 0.05, "Yes", "No")
)
# Compute difference between means
mean_diff <- rowMeans(mat.ori_vd) - rowMeans(mat.ori_vehicle)

# Add difference between medians to result_df
result_df_all$Mean_Difference <- mean_diff
result_df_all <- result_df_all %>% filter(Significant == "Yes") %>% arrange(desc(Mean_Difference))
#Add red color for upregulated by VD and blue for downregulated
result_df_all$Color <- ifelse(result_df_all$Mean_Difference > 0, "red", "blue")
#Heatmap altered function
netAnalysis_signalingRole_heatmap_ESD <- function (object, signaling = NULL, pattern = c("outgoing", "incoming", 
                                                                                         "all"), slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
                                                   title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
                                                   cluster.rows = FALSE, cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  print(unique(is.na(mat)))
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[is.na(mat)] <- 0
  mat[is.infinite(mat)] <- 0
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat)#.ori)
  pSum.original <- pSum
  #  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  print(unique(is.na(mat)))
  ht1 = Heatmap(mat, col = color.heatmap.use, 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.cols, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size, col = result_df_all$Color), 
                column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                           "cm"), height = unit(height, "cm"), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                 "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                                                          "mm")))
  return(ht1)
}

#heatmap
ht1 <- netAnalysis_signalingRole_heatmap_ESD(cellchat_vehicle, pattern = "all", signaling = result_df_all$Row, title = "vehicle", width = 5, height = 12, cluster.cols = T, cluster.rows = F)
ht2 <- netAnalysis_signalingRole_heatmap_ESD(cellchat_vd, pattern = "all", signaling = result_df_all$Row, title = "vd", width = 5, height = 12, cluster.cols = T, cluster.rows = F)
png(filename = file.path(plots.dir, paste0(sample.name, "_ht_diffpath.png")),
    width = 7, height = 8.25, units = "in", res=300)
draw(ht1 + ht2)
dev.off()
#### to check the freaking L-R as well jajaja

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in vehicle", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in vd", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
#png(filename = file.path(plots.dir, paste0(sample.name, "_bubbleL-R.png")),
#    width = 15, height = 10, units = "in", res=300)
gg1 + gg2
#dev.off()

########## Save objects ###############

saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_scall.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_scall.rds")))
