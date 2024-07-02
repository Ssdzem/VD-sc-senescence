################# instructions #########
 # need to adjust lines 35-50 and 495-507 for input directories (data) and annotation data (see ann_codes)
########################################
# The standard Seurat workflow takes raw single-cell expression data and
# aims to find clusters within the data.
#
# This process consists of:
# data normalization and variable feature selection
# data scaling
# cell cycle scoring and regression
# a PCA on variable features
# construction of a shared-nearest-neighbors graph
# clustering using a modularity optimizer
# Finally, visualization of clusters in a two-dimensional space
#
#################### import libraries and set options ##################

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(sva))
suppressMessages(library(SeuratData))
suppressMessages(library(ggpubr)) # apparently you need the latest dev version
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
suppressMessages(library(SingleR))
suppressMessages(library(celldex))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(biomaRt))
set.seed(12345)

################################  set up variables #####################
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/hanai/vd"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"
#In lin_2022  input it does not count since i have already the gene matrix count, so i developed another script called "readmatrix_lin.R"

sample.name = "hanai_vd"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"

#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/hanai/vd/"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"

########################## functions ###################################
# it's dangerous to go alone! take this.
#
# function to match normal gene names to my weird gene names
#get_gene_id <- function(q, ids=rownames(pbmc.seurat)) {
#  paste0("hrg", q, "$") %>% grepl(ids) %>%
#    which() -> i
#  
#  return(ids[i])
#}

# create DimPlot for cell cycle anaylsis
plot_dp_cc <- function(s=pbmc.seurat.cc, t) {
  
  DimPlot(s, pt.size = 1.55, group.by = "Phase", order = c("G2M", "S", "G1")) + 
    ggtitle(t) +
    scale_color_manual(values = c("G1" = "#0075DC", "S"="#9DCC00","G2M"="#FF5005")) +
    theme(legend.position="top")
}

# source multibar heatmap function
#  source("workflow/scripts/multibar_heatmap_function.R") #!!

##########################  read in matrix ################################
options(Seurat.object.assay.version = "v3") #SIngleR uses V3
pbmc.seurat <- Read10X(input.seurat)
pbmc.seurat <- CreateSeuratObject(counts = pbmc.seurat, min.cells = 3, min.features = 200, project = 'vd')
class(pbmc.seurat[["RNA"]])
######################### Quality control #################################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

###########################################################################
########################## Normalization ##################################
###########################################################################
# "The number of useful reads obtained from a sequencing experiment will
# vary between cells, and one must correct for this difference"
#
# "Normalization ensures that any observed heterogeneity or differential expression 
# within the cell population are driven by biology and not technical biases."
#
# LogNormalize for each cell:
# divide feature expression measurements by the total expression
# multiply the result by a scale factor (default is 10,000)
# log-transform the result 
#
# The log-transformation is useful because 
# differences in the log-values represent log-fold changes in expression
#
# Assumption:
# any cell-specific bias (capture or amplification efficiency) 
# affects all genes equally via scaling of the expected mean count for that cell
#
########## obtain log-transformed normalized expression values
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)


### SPLIT HERE

########## Feature selection
# Feature selection identifies genes with the strongest biological signal relative to the technical noise
# cell clustering is dependent on the selection of these genes (to compare the gene expression of cells, 
# we have to aggregate all of the gene's differences into a single (dis)similarity metric between pairs of cells)
#
# we have to select genes that contain useful info about the biology of the system 
# we have to ignore the genes whose expression differences is actually random noise
# we have to reduce the size of the data to improve computational efficiency
#
# it is common to identify highly variable features (genes) for this.
# we select the most variable genes (based on their expression across the population)
# because true biological differences manifest as increased variation 
# as opposed to genes affected by technical noise or a baseline uninteresting biological variation

# obtain genes with a higher-than-expected variance
# (Seurat fits empirically the relationship between variance and mean expression)
pbmc.seurat <- FindVariableFeatures(pbmc.seurat, selection.method = "vst", nfeatures = 2000)

# variable features plot
top50 <- head(VariableFeatures(pbmc.seurat), 50)
top50.colors <- rep("black", 50)
top50.colors[!grepl("hrg", top50)] <- "lightslateblue"

p1 <- VariableFeaturePlot(pbmc.seurat, pt.size = 1.5)

p2 <- LabelPoints(plot = p1,
                  points = top50, labels = gsub(pattern = ".+hrg", "", top50),
                  color	= top50.colors,
                  repel = TRUE,
                  segment.color =	"grey",
                  segment.size = 0.25,
                  segment.alpha = 0.85) +
  
  theme(legend.text = element_text(size = 15),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.position="bottom",
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        title = element_text(size=22)) +
  
  guides(colour = guide_legend(override.aes = list(size=4))) +
  
  ggtitle(paste("Variable Features (transcripts)", sample.name))

# I could change the axes to have more ticks
#max(pbmc.seurat@assays$RNA@meta.features$vst.variance.standardized)

png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()

rm(top50, top50.colors, p1,p2)

###########################################################################
########################## Data Scaling ###################################
###########################################################################
# Scales the expression of each gene to have a variance of 1
# so all genes have equal contributions in downstream analyses
#
# All gene expression will be centered to 0

# this would be to scale only highly variable genes (the ones used for the PCA)
#pbmc.seurat <- ScaleData(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

# but since we are planning to do a heatmap, that is also built on scaled genes
# we scale all genes 
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))

###########################################################################
################## cell cycle scoring and regression ######################
###########################################################################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# We can segregate this list into markers of G2/M phase and markers of S phase
# and match them to my IDs
#cc.genes.updated.2019$s.genes -> sh.genes

#cc.genes.updated.2019$g2m.genes -> g2h.genes 

#If working with other species than human, you need to use biomart to change the names of the genes

#we use archived databases because biomart is kind of yucky
#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#g2hm <- getLDS(attributes = c("external_gene_name"),
#       filters = "external_gene_name", values = g2h.genes, mart = human,
#       attributesL = c("external_gene_name"), martL = mouse)

#g2m.genes <- g2hm[,2]

#shm.genes <-  getLDS(attributes = c("external_gene_name"),
#       filters = "external_gene_name", values = sh.genes, mart = human,
#       attributesL = c("external_gene_name"), martL = mouse)

#sm.genes <- shm.genes[,2]

# After scoring each gene for cell cycle phase, we can perform PCA using the expression of cell cycle genes. 
# If the cells group by cell cycle in the PCA, then we would want to regress out cell cycle variation,
# unless cells are differentiating

########## Assign Cell-Cycle Scores
# assign each cell a score, based on its expression of G2/M and S phase markers
# the marker sets should be anti correlated in their expression levels
# cells expressing neither are likely not cycling and in G1 phase.
#pbmc.seurat <- CellCycleScoring(pbmc.seurat, s.features = sm.genes, g2m.features = g2m.genes, set.ident = TRUE)

##### Visualize the distribution of cell cycle markers across
#c("Pcna", "Top2a", "Mcm6", "Mki67") -> cc.features

#p <- RidgePlot(pbmc.seurat, features = cc.features, combine = TRUE, ncol=2)

#p <- lapply(seq(1:4), function(x) {
#  p[[x]] <- p[[x]] + labs(title = gsub(pattern = ".+hrg", "", cc.features)[x]) + 
#    ylab("") +
#    scale_fill_manual(values = c("G1" = "#0075DC", "S"="#9DCC00","G2M"="#FF5005")) 
#})

#png(filename=file.path(plots.dir, paste0(sample.name, "_CC_markers_ridge_plot.png")),width=20, height=16, units="in", res=300)
# print((p[[1]] | p[[2]]) / (p[[3]] | p[[4]]))
# dev.off()

##### BEFORE plots
# Run a PCA on variable features
#pbmc.seurat.cc <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))
#a <- plot_dp_cc(t = paste0("PCA on Variable Features\n", sample.name))
#a <- plot_dp_cc(t = "PCA on Variable Features")
#a
# Run a PCA on cell cycle genes to see if cells separate by phase
#pbmc.seurat.cc <- RunPCA(pbmc.seurat.cc, features = c(sm.genes, g2m.genes))
#b <- plot_dp_cc(t = paste0("PCA on Cell Cycle genes\n", sample.name))
#b <- plot_dp_cc(t = "PCA on Cell Cycle genes")
#b
########## Regress out cell cycle scores during data scaling
#pbmc.seurat.cc <- ScaleData(pbmc.seurat.cc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc.seurat.cc))

##### AFTER plots
# Now, a PCA on the variable genes should not  return components associated with cell cycle
#pbmc.seurat.cc <- RunPCA(pbmc.seurat.cc, features = VariableFeatures(object = pbmc.seurat.cc))
#c <- plot_dp_cc(t = paste0("PCA on Variable Features (Regressed cell cycle scores)\n", sample.name))
#c <- plot_dp_cc(t = "PCA on Variable Features (Regressed cell cycle scores)")
#c
# In a PCA on only cell cycle genes, cells should no longer separate by cell-cycle phase
#pbmc.seurat.cc <- RunPCA(pbmc.seurat.cc, features = c(sm.genes, g2m.genes))
#d <- plot_dp_cc(t = paste0("PCA on Cell Cycle genes (Regressed cell cycle scores)\n", sample.name))
#d <- plot_dp_cc(t = "PCA on Cell Cycle genes (Regressed cell cycle scores)")
#d
# create plot
# png(filename=file.path(plots.dir, paste0(sample.name, "_CellCycle_dim_plots.png")), width=24, height=18, units="in", res=300)
# print((a | b ) / (c | d ) )
# dev.off()

##### barplot of cells in each phase
#pbmc.md <- pbmc.seurat@meta.data[c("Phase")]

#cc.bar <- ggplot(data = pbmc.md, aes(Phase, fill=Phase)) + geom_bar() +
#  scale_fill_manual(values = c("G1" = "#0075DC", "S"="#9DCC00","G2M"="#FF5005")) +
#  theme_minimal() + xlab("Cell Cycle Phase") + ylab("Number of Cells")
#cc.bar
#up <- (a | b | (  ( p[[1]] | p[[2]] ) / ( p[[3]] | p[[4]] ) )) 
#down <- (c | d | (cc.bar))

##### one big plot
#png(filename=file.path(plots.dir, paste0(sample.name, "_CellCycle_plots.png")),
#    width=22, height=13, units="in", res=300)

#print(up/down + 
        #         plot_layout( guides = "collect") + 
#        plot_annotation(title = paste0("\n",sample.name,"\n"), 
#                        theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5))) 
#)

#dev.off()

# clean up
#rm(s.genes, g2m.genes, p, a, b, c, d, cc.features, pbmc.seurat.cc, cc.bar, pbmc.md, up, down,human,mouse, shm.genes, sh.genes, sm.genes, g2hm)

###########################################################################
################## linear dimensional reduction ###########################
###########################################################################
# Dimensionality reduction aims to reduce the number of separate dimensions in the data
# This is possible because different genes are correlated if they are affected by the same biological process
#
# We compress multiple features into a single dimension (an “eigengene”)
# this reduces computational work (perform calculations only for a few dimensions rather than thousands of genes
# this reduces noise by averaging across multiple genes to obtain a more precise representation of the patterns in the data
# this enables effective plotting of the data (humans can't see more than 3 dimensions)
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

# Principal components analysis (PCA) discovers axes in high-dimensional space that capture the largest amount of variation

##### dim plot
dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

##### elbow plot
eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

##### Visualize top genes associated with reduction components
vd.plots <- VizDimLoadings(pbmc.seurat, dims = 1:6, reduction = "pca", combine = FALSE)

vd.plots <- lapply(X = vd.plots, FUN = function(x) {
  s <- x$data$feature
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_features_contribution.png")),
    width = 22, height = 14, units = "in", res=300)

patchwork::wrap_plots(vd.plots, ncol=3) +
  plot_annotation(title = paste0("\n", sample.name, " - Contribution of features (transcripts) to Principal Components\n"),
                  theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)) )

dev.off()

rm(vd.plots)

##### Dimensional reduction heat map
# Both cells and features are ordered according to their PCA scores. 
# Though clearly a supervised analysis, it is a valuable tool for exploring correlated feature sets.
hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)

# might turn this into a function
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
# you can't know how much I hate donig this indexing trick
hm.plots <- lapply(seq(1:15), FUN = function(x) {
  hm.plots[[x]] <- hm.plots[[x]] + ggtitle(paste0("PC_", x))
  
})

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_heatmap.png")),
    width = 26, height = 32, units = "in", res=300)

patchwork::wrap_plots(hm.plots, ncol = 3) + 
  plot_annotation(title = paste0("\n", sample.name, " - Contribution of features (transcripts) to Principal Components\n"),
                  theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)) )

dev.off()
rm(hm.plots)

##### Select PCs
# By definition, the top PCs capture the dominant factors of heterogeneity in the data set.
# In the context of scRNA-seq, our assumption is that biological processes affect multiple genes in a coordinated manner.
#
# This means that the earlier PCs are likely to represent biological structure
# as more variation can be captured by considering the correlated behavior of many genes. 
#
# By comparison, random technical or biological noise is expected to affect each gene independently. 
# There is unlikely to be an axis that can capture random variation across many genes, 
# meaning that noise should mostly be concentrated in the later PCs. 
#
# This motivates the use of the earlier PCs in our downstream analyses, 
# which concentrates the biological signal to simultaneously reduce computational work and remove noise.
#
# I will select 10 PCs, but I am plotting the JackStraw analysis to assess this decision 
#
# randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, 
# and repeat this procedure. Identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

# Significant PCs should show a p-value distribution (black curve) that is strongly skewed to the left 
# compared to the null distribution (dashed line).
# The p-value for each PC is based on a proportion test comparing the number of genes with a p-value 
# below a particular threshold (score.thresh),
# compared with the proportion of genes expected under a uniform distribution of p-values.

###########################################################################
########## Shared nearest neighbor graph  and clustering ##################
###########################################################################
# Seurat first constructs a KNN graph based on the euclidean distance in PCA space. 
# Each node is a cell (connected to its nearest neighbors by edges)
# Edges are weighted based on the similarity between the cells involved, 
# with higher weight given to cells that are more closely related. 
# A refinement step (Jaccard similarity) is used to refine the edge weights between any two cells
# based on the shared overlap in their local neighborhoods
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)

# partition the graph into highly interconnected quasi-cliques/communities using Louvain.
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
# the resolution parameter sets the ‘granularity’ of the downstream clustering, 
# with increased values leading to a greater number of clusters. 

# "clustering is an unsupervised learning procedure that is used to empirically define groups of cells 
# with similar expression profiles. 
# Its primary purpose is to summarize complex scRNA-seq data into a digestible format
# for human interpretation."

###########################################################################
######################## Visualization in 2D space ########################
###########################################################################
pbmc.seurat <- RunUMAP(pbmc.seurat, dims = 1:10)

p.umap <- DimPlot(pbmc.seurat, reduction = "umap", label = TRUE, label.size = 5, group.by = pbmc.seurat@meta.data[["SingleR.labels"]]) +
  ggtitle(paste0("Clustering UMAP\n", sample.name))
print(p.umap)
p.pca <- DimPlot(pbmc.seurat, reduction = "pca", label = TRUE) +
  ggtitle(paste0("Clustering PCA\n", sample.name))

png(filename = file.path(plots.dir, paste0(sample.name, "_clusters_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(p.pca + p.umap)

dev.off()

# In scRNA-seq clusters are groups of cells with similar transcriptomes
###########################################################################
############################## Gene Markers ############################### 
###########################################################################
# find all markers means for each cluster, find markers that make it different
# to all of the other cells.
# But markers may bleed over between closely-related clusters because
# they are not forced to be specific to only one cluster.
#
# Identify the transcripts that drive separation of the clusters
# By default, FindAllMarkers uses Wilcoxon rank-sum (Mann-Whitney-U) test
# to find Differentially Expressed transcripts. 
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)

# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells
# logfc.threshold: limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
#
#A positive marker is a gene that is overexpressed in a cell population compared to other populations.
#A negative marker is a gene that is underexpressed in a cell population compared to other populations.

#ggplot(data = pbmc.markers, aes(p_val)) + geom_histogram()+ theme_minimal()
#ggplot(data = pbmc.markers, aes(p_val_adj)) + geom_histogram()+ theme_minimal()

# select significant markers 
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)

# add column to identify gene type
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene),"Protein Coding", "Transposable Element"))
# export markers
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

# Beware that comparing each cluster to the average of all other cells is sensitive to
# the population composition, which introduces
# an element of unpredictability to the marker sets due to variation in cell type abundances
#
# Consider using pairwise comparisons instead
###########################################################################
################## Assign Cell type identity ##############################
###########################################################################
ref <- pbmc.boneM@assays[["RNA"]]@counts
test <- pbmc.seurat@assays[["RNA"]]@counts
#run SingleR using the expression matrix of the seurat object
test_singleR <- SingleR( test = test, ref = ref, labels = pbmc.boneM@meta.data[["CellType"]], num.threads = 60)

#move the results to the seurat obj
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels

#rm remaining data
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)

###########################################################################
############################# Cell annotation #############################
###########################################################################
#Print it, hell yeah
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

###########################################################################
############################# Explore markers #############################
###########################################################################
# Cells by cluster

uno <- DimPlot(pbmc.seurat, reduction = "umap", label = TRUE) +
  ggtitle("Seurat Clusters")

dos <- ggplot(data = pbmc.markers, aes(cluster, fill=gene_type)) +
  geom_bar(position=position_dodge()) +
  theme_minimal() +
  ggtitle("Significant markers by cluster") + 
  xlab("Cluster") + ylab("Number of Markers") 
  

tres <- ggplot(data = pbmc.seurat@meta.data, aes(seurat_clusters, fill=SingleR.labels)) + geom_bar() +
  ylab("Number of cells") + xlab("Seurat Clusters") + theme_minimal() +
  scale_fill_manual(values = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  ggtitle("Number of Cells and Cell Type by Cluster")

png(filename = file.path(plots.dir, paste0(sample.name, "_markers_barplot.png")),
    width = 20, height = 14, units = "in", res=300)

print((uno | dos) / (tres) +
        plot_annotation(title = paste0("\n",sample.name,"\n"),
                        theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)))
)

dev.off()

rm(uno, dos, tres)

##### top 25 markers by cluster
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC) %>% 
  as.data.frame() -> top25.markers

## Feature plots
wow <-  lapply(split.data.frame(top25.markers, f = top25.markers$cluster), function(f) {
  
  f <- f[order(f$avg_log2FC, decreasing = TRUE), ]
  
  t25 <- f$gene
  t25.colors <- rep("black", 25)
  t25.colors[!grepl("hrg", t25)] <- "lightslateblue"
  
  l <- FeaturePlot(pbmc.seurat, features = t25, combine = FALSE) 
  
  l <- lapply(seq(1:25), function(i) {
    l[[i]] <- l[[i]] + ggtitle(gsub(".+hrg", "", t25[[i]])) +
      theme(plot.title = element_text(colour = t25.colors[i], face = "bold")) 
  })
  
  mine <- patchwork::wrap_plots(l, ncol = 5) +
    plot_annotation(title = paste0("\n", sample.name, "   Cluster ", unique(f$cluster), " - Top 25 (log2FC) Markers\n"),
                    theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)) )
  
  return(mine)
  
})

pdf(file.path(plots.dir, paste0(sample.name, "_markers_by_cluster.pdf")), width = 30, height = 28)
dev.off()

# I really called an object wow because it is 1am wow
rm(wow)

## Markers Heatmap

#h <- DoHeatmap(pbmc.seurat, features = top25.markers$gene) 

#h <- DoMultiBarHeatmap(pbmc.seurat, features = top25.markers$gene, label = FALSE, group.by = "seurat_clusters")

#h.labs = ggplot_build(h)$layout$panel_params[[1]]$y$get_labels()

#h.cols = rep("black", length(h.labs))
#h.cols[!grepl("hrg", h.labs)] <- "lightslateblue"

#h.labs = gsub(".+hrg", "", h.labs)

#h <- h + scale_y_discrete(labels=h.labs) + theme(axis.text.y = element_text(colour = h.cols)) +
#ggtitle(paste0(sample.name, "  -  markers by cluster")) +
#theme(title = element_text(size=20))

#png(filename = file.path(plots.dir, paste0(sample.name, "_markers_heatmap.png")),
#width = 25, height = 25, units = "in", res=300)

#print(h)

#dev.off()

#rm(h, h.labs, h.cols)



###########################################################################
############################### save object ############################### 
###########################################################################
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))


