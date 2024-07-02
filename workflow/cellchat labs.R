############ Outgoing pathways #############
#vehicle
centr <- slot(cellchat_vehicle, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat_vehicle@idents), names(centr))
for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]]$outdeg}
mat <- t(outgoing)
mat1 <- mat[rownames(mat) %in% pathway.union, , drop = FALSE]
mat <- matrix(0, nrow = length(pathway.union), ncol = ncol(mat))
idx <- match(rownames(mat1), pathway.union)
mat[idx[!is.na(idx)], ] <- mat1
dimnames(mat) <- list(pathway.union, colnames(mat1))
mat.ori_vehicle <- mat

#vd
centr <- slot(cellchat_vd, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat_vehicle@idents), names(centr))
for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]]$outdeg}
mat <- t(outgoing)
mat1 <- mat[rownames(mat) %in% pathway.union, , drop = FALSE]
mat <- matrix(0, nrow = length(pathway.union), ncol = ncol(mat))
idx <- match(rownames(mat1), pathway.union)
mat[idx[!is.na(idx)], ] <- mat1
dimnames(mat) <- list(pathway.union, colnames(mat1))
mat.ori_vd <- mat
#JAJAJAJAJ this shit is for the column annotation of the heatmap
Strength = anno_barplot(colSums(mat.ori), 
                        border = FALSE, gp = gpar(fill = color.use, col = color.use))

#here is the actual code for the row anno
pSum <- rowSums(mat.ori)
pSum.original <- pSum
pSum <- -1/log(pSum)
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
result_df_outgoing <- data.frame(
  Row = rownames(mat.ori_vd),
  P_Value = p_values,
  Significant = ifelse(p_values < 0.05, "Yes", "No")
)
# Compute difference between means
median_diff <- rowMeans(mat.ori_vd) - rowMeans(mat.ori_vehicle)

# Add difference between medians to result_df
result_df_outgoing$Mean_Difference <- mean_diff
print(result_df_outgoing %>% filter(Significant == "Yes"))
             
########## Ingoing pathways ###############

#vehicle
centr <- slot(cellchat_vehicle, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat_vehicle@idents), names(centr))
for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]]$indeg}
mat <- t(outgoing)
mat1 <- mat[rownames(mat) %in% pathway.union, , drop = FALSE]
mat <- matrix(0, nrow = length(pathway.union), ncol = ncol(mat))
idx <- match(rownames(mat1), pathway.union)
mat[idx[!is.na(idx)], ] <- mat1
dimnames(mat) <- list(pathway.union, colnames(mat1))
mat.ori_vehicle <- mat

#vd
centr <- slot(cellchat_vd, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(cellchat_vehicle@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat_vehicle@idents), names(centr))
for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]]$outdeg}
mat <- t(outgoing)
mat1 <- mat[rownames(mat) %in% pathway.union, , drop = FALSE]
mat <- matrix(0, nrow = length(pathway.union), ncol = ncol(mat))
idx <- match(rownames(mat1), pathway.union)
mat[idx[!is.na(idx)], ] <- mat1
dimnames(mat) <- list(pathway.union, colnames(mat1))
mat.ori_vd <- mat
#JAJAJAJAJ this shit is for the column annotation of the heatmap
Strength = anno_barplot(colSums(mat.ori), 
                        border = FALSE, gp = gpar(fill = color.use, col = color.use))

#here is the actual code for the row anno
pSum <- rowSums(mat.ori)
pSum.original <- pSum
pSum <- -1/log(pSum)
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
result_df_ingoing <- data.frame(
  Row = rownames(mat.ori_vd),
  P_Value = p_values,
  Significant = ifelse(p_values < 0.05, "Yes", "No")
)
# Compute difference between means
median_diff <- rowMeans(mat.ori_vd) - rowMeans(mat.ori_vehicle)

# Add difference between medians to result_df
result_df_ingoing$Mean_Difference <- mean_diff
print(result_df_ingoing %>% filter(Significant == "Yes"))

########### Overall ####################
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

#JAJAJAJAJ this shit is for the column annotation of the heatmap
#Strength = anno_barplot(colSums(mat.ori), 
#                       border = FALSE, gp = gpar(fill = color.use, col = color.use))

#here is the actual code for the row anno
pSum <- rowSums(mat.ori)
pSum.original <- pSum
pSum <- -1/log(pSum)
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
print(result_df_all %>% filter(Significant == "Yes"))

########## just to see why the values so crazy vehicle vs vitamin D overrall, response:bad log #####
#vehicle
pSum <- rowSums(mat.ori_vehicle)
pSum.original <- pSum
pSum <- -1/log(pSum)
pSum[is.na(pSum)] <- 0
idx1 <- which(is.infinite(pSum) | pSum < 0)
if (length(idx1) > 0) {
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                       length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
}
pSumvehicle <- pSum

#vd
pSum <- rowSums(mat.ori_vd)
pSum.original <- pSum
pSum <- -1/log(pSum)
pSum[is.na(pSum)] <- 0
idx1 <- which(is.infinite(pSum) | pSum < 0)
if (length(idx1) > 0) {
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                       length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
}
pSumvd <- pSum



########## why i cant cluster fucking NA ###############

df <- data.frame(group = colnames(mat))
heatmap(mat)
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
htprove1 <- tpm_marge(cellchat_vehicle, pattern = "all", signaling = result_df_all$Row, title = "vehicle", width = 5, height = 12, cluster.cols = T, cluster.rows = F)
htprove2 <- tpm_marge(cellchat_vd, pattern = "all", signaling = result_df_all$Row, title = "vd", width = 5, height = 12, cluster.cols = T, cluster.rows = F)
draw(htprove1 + htprove2)




#### and to see the differential significance of the cell types R: NO SIGNIFICANT DIFF WITH ALL MAT SUBSETS (normalized, all pathways, diff pathways...)

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
mat.vehicle_celltype <- mat
mat.vehicle_celltype <- mat.vehicle_celltype[result_df_all$Row, ]
mat.vehicle_celltype <- sweep(mat.vehicle_celltype, 1L, apply(mat.vehicle_celltype, 1, max), "/", check.margin = FALSE)
mat.vehicle_celltype[is.na(mat.vehicle_celltype)] <- 0
mat.vehicle_celltype[is.infinite(mat.vehicle_celltype)] <- 0


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
mat.vd_celltype <- mat
mat.vd_celltype <- mat.vd_celltype[result_df_all$Row, ]
mat.vd_celltype <- sweep(mat.vd_celltype, 1L, apply(mat.vd_celltype, 1, max), "/", check.margin = FALSE)
mat.vd_celltype[is.na(mat.vd_celltype)] <- 0
mat.vd_celltype[is.infinite(mat.vd_celltype)] <- 0


# Initialize vectors and lists
p_values <- numeric(ncol(mat.vd_celltype))
test_results <- list()

# Perform Wilcoxon signed-rank test for each column pair
for (i in 1:ncol(mat.vd_celltype)) {
  vd_col <- mat.vd_celltype[, i]
  vehicle_col <- mat.vehicle_celltype[, i]
  test_result <- wilcox.test(vd_col, vehicle_col, paired = T)
  p_values[i] <- test_result$p.value
  test_results[[i]] <- test_result
}

# Display results
result_df_all_celltype <- data.frame(
  Column = colnames(mat.vd_celltype),
  P_Value = p_values,
  Significant = ifelse(p_values < 0.05, "Yes", "No")
)
# Compute difference between means
median_diff <- colMeans(mat.vd_celltype) - colMeans(mat.vehicle_celltype)

# Add difference between medians to result_df
result_df_all_celltype$Mean_Difference <- median_diff

# Print significant results
print(result_df_all_celltype)
