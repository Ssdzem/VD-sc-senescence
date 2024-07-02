# 1. Instructions ####
# this is a code to run everything in order to make the experiment the most reproducible possible
#this will be on read only, major changes will be on individual code of workflow
# 2. Prepros ####
#  2.1 libraries ####
suppressMessages(library(Seurat))
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))
suppressMessages(library(tidyverse))
suppressMessages(library(vroom))
suppressMessages(library(patchwork))
suppressMessages(library(sva))
suppressMessages(library(SeuratData))
suppressMessages(library(ggpubr))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
set.seed(12345)

#  2.2 QC,processing and annotation ####
#   2.2.1 aem ####
#    2.2.1.1 aem_vd ####
#     2.2.1.1.1 Values ####
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/aem/vd"
sample.name = "aem_vd"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/aem/vd/"

#     2.2.1.1.2 load ####
options(Seurat.object.assay.version = "v3") #SIngleR uses V3
pbmc.seurat <- Read10X(input.seurat)
pbmc.seurat <- CreateSeuratObject(counts = pbmc.seurat, min.cells = 3, min.features = 200, project = 'vd')
class(pbmc.seurat[["RNA"]])
#     2.2.1.1.3 Quality control ####
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

#     2.2.1.1.4 Normalization ####
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#     2.2.1.1.5 Feature selection ####
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
png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()
rm(top50, top50.colors, p1,p2)

#     2.2.1.1.6 
#     2.2.1.1.6 Data scaling ####
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))
#     2.2.1.1.7 Dimensional reduction ####
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

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

hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
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

pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

#     2.2.1.1.8 SNN graph and clustering ####
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
#     2.2.1.1.9 UMAP visualization ####
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
#     2.2.1.1.10 Markers of each cluster ####
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene), "Protein Coding", "Transposable Element"))
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

#     2.2.1.1.11 Annotation ####
barcodes_annotation <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_barcodes_anno.csv")
dge <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_dge.csv")
gene_list <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_gene.csv")
ref <- merge(gene_list, dge, by = "...1")
ref <- ref[,2:8289]
genes <- ref[,1]
rownames(ref) <- genes
ref <- ref[,2:8288]
ref <- as.matrix(ref)
test <- pbmc.seurat@assays[["RNA"]]@counts
test_singleR <- SingleR( test = test, ref = ref, labels =  barcodes_annotation$Idents.pbmc., num.threads = 60)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

#     2.2.1.1.12 Explore markers ####
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

#     2.2.1.1.13 Save ####
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))
#     2.2.1.1.14 clean ####
rm(list=ls())
gc()
#    2.2.1.2 aem_vehicle ####
#     2.2.1.2.1 Values ####
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/aem/vehicle"
sample.name = "aem_vehicle"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/aem/vehicle/"

#     2.2.1.2.2 load ####
options(Seurat.object.assay.version = "v3") #SIngleR uses V3
pbmc.seurat <- Read10X(input.seurat)
pbmc.seurat <- CreateSeuratObject(counts = pbmc.seurat, min.cells = 3, min.features = 200, project = 'vd')
class(pbmc.seurat[["RNA"]])
#     2.2.1.2.3 Quality control ####
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

#     2.2.1.2.4 Normalization ####
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#     2.2.1.2.5 Feature selection ####
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
png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()
rm(top50, top50.colors, p1,p2)

#     2.2.1.1.6 
#     2.2.1.2.6 Data scaling ####
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))
#     2.2.1.2.7 Dimensional reduction ####
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

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

hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
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

pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

#     2.2.1.2.8 SNN graph and clustering ####
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
#     2.2.1.2.9 UMAP visualization ####
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
#     2.2.1.2.10 Markers of each cluster ####
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene), "Protein Coding", "Transposable Element"))
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

#     2.2.1.2.11 Annotation ####
barcodes_annotation <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_barcodes_anno.csv")
dge <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_dge.csv")
gene_list <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_gene.csv")
ref <- merge(gene_list, dge, by = "...1")
ref <- ref[,2:8289]
genes <- ref[,1]
rownames(ref) <- genes
ref <- ref[,2:8288]
ref <- as.matrix(ref)
test <- pbmc.seurat@assays[["RNA"]]@counts
test_singleR <- SingleR( test = test, ref = ref, labels = barcodes_annotation$Idents.pbmc., num.threads = 60)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

#     2.2.1.2.12 Explore markers ####
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

#     2.2.1.2.13 Save ####
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))
#     2.2.1.2.14 clean ####
rm(list=ls())
gc()
#   2.2.2 Hanai ####
#    2.2.2.1 hanai_vd ####
#     2.2.2.1.1 Values ####
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/hanai/vd"
sample.name = "hanai_vd"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/hanai/vd/"

#     2.2.2.1.2 load ####
options(Seurat.object.assay.version = "v3") #SIngleR uses V3
pbmc.seurat <- Read10X(input.seurat)
pbmc.seurat <- CreateSeuratObject(counts = pbmc.seurat, min.cells = 3, min.features = 200, project = 'vd')
class(pbmc.seurat[["RNA"]])
#     2.2.2.1.3 Quality control ####
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

#     2.2.2.1.4 Normalization ####
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#     2.2.2.1.5 Feature selection ####
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
png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()
rm(top50, top50.colors, p1,p2)

#     2.2.1.1.6 
#     2.2.2.1.6 Data scaling ####
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))
#     2.2.2.1.7 Dimensional reduction ####
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

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

hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
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

pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

#     2.2.2.1.8 SNN graph and clustering ####
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
#     2.2.2.1.9 UMAP visualization ####
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
#     2.2.2.1.10 Markers of each cluster ####
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene), "Protein Coding", "Transposable Element"))
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

#     2.2.2.1.11 Annotation ####
pbmc.bone <- Connect(filename = "/datos/sensence/emilio/VD/input/ref/hanai/SCA_v.0.5.0_preprint.loom", mode = "r")
pbmc.bone <- as.Seurat(pbmc.bone)
selected_cell_types <- c("Osteoclasts", "Osteoblasts","Chondrocytes", "Endothelium", "LepR+ BM-MSCs", "Fibroblasts", "Osteocytes", "Pericytes","Periosteum", "Adipocytes", "Immune cells", "Vascular smooth muscle", "Tenocytes")
pbmc.bone <- subset(pbmc.bone, subset = CellType %in% selected_cell_types)
ref <- pbmc.bone@assays[["RNA"]]@counts
test <- pbmc.seurat@assays[["RNA"]]@counts
test_singleR <- SingleR( test = test, ref = ref, labels = pbmc.bone@meta.data[["CellType"]], num.threads = 60)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

#     2.2.2.1.12 Explore markers ####
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

#     2.2.2.1.13 Save ####
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))
#     2.2.2.1.14 clean ####
rm(list=ls())
gc()
#    2.2.2.2 hanai_vehicle ####
#     2.2.2.2.1 Values ####
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/hanai/vehicle"
sample.name = "hanai_vehicle"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/hanai/vehicle/"

#     2.2.2.2.2 load ####
options(Seurat.object.assay.version = "v3") #SIngleR uses V3
pbmc.seurat <- Read10X(input.seurat)
pbmc.seurat <- CreateSeuratObject(counts = pbmc.seurat, min.cells = 3, min.features = 200, project = 'vd')
class(pbmc.seurat[["RNA"]])
#     2.2.2.2.3 Quality control ####
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

#     2.2.2.2.4 Normalization ####
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#     2.2.2.2.5 Feature selection ####
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
png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()
rm(top50, top50.colors, p1,p2)

#     2.2.1.1.6 
#     2.2.2.2.6 Data scaling ####
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))
#     2.2.2.2.7 Dimensional reduction ####
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

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

hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
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

pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

#     2.2.2.2.8 SNN graph and clustering ####
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
#     2.2.2.2.9 UMAP visualization ####
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
#     2.2.2.2.10 Markers of each cluster ####
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene), "Protein Coding", "Transposable Element"))
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

#     2.2.2.2.11 Annotation ####
pbmc.bone <- Connect(filename = "/datos/sensence/emilio/VD/input/ref/hanai/SCA_v.0.5.0_preprint.loom", mode = "r")
pbmc.bone <- as.Seurat(pbmc.bone)
selected_cell_types <- c("Osteoclasts", "Osteoblasts","Chondrocytes", "Endothelium", "LepR+ BM-MSCs", "Fibroblasts", "Osteocytes", "Pericytes","Periosteum", "Adipocytes", "Immune cells", "Vascular smooth muscle", "Tenocytes")
pbmc.bone <- subset(pbmc.bone, subset = CellType %in% selected_cell_types)
ref <- pbmc.bone@assays[["RNA"]]@counts
test <- pbmc.seurat@assays[["RNA"]]@counts
test_singleR <- SingleR( test = test, ref = ref, labels = pbmc.bone@meta.data[["CellType"]], num.threads = 60)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

#     2.2.2.2.12 Explore markers ####
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

#     2.2.2.2.13 Save ####
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))
#     2.2.2.2.14 clean ####
rm(list=ls())
gc()
#   2.2.3 lin ####
#    2.2.3.1 lin_vd ####
#     2.2.3.1.1 Values ####
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/lin/vd"
sample.name = "lin_vd"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/lin/vd/"

#     2.2.3.1.2 load ####
VD_path = "/datos/sensence/emilio/VD/input/sc_data/lin/vd/GSM5266944_VD_matrix.tsv"
data_VD <- read.delim(VD_path, header = TRUE, sep = '\t', row.names = "X")
counts_matrix = Matrix(as.matrix(data_VD) , sparse=TRUE)
options(Seurat.object.assay.version = "v3")
pbmc.seurat <- CreateSeuratObject(counts = counts_matrix, min.cells = 3, min.features = 200)
class(pbmc.seurat[["RNA"]])
#     2.2.3.1.3 Quality control ####
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

#     2.2.3.1.4 Normalization ####
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#     2.2.3.1.5 Feature selection ####
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
png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()
rm(top50, top50.colors, p1,p2)

#     2.2.1.1.6 
#     2.2.3.1.6 Data scaling ####
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))
#     2.2.3.1.7 Dimensional reduction ####
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

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

hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
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

pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

#     2.2.3.1.8 SNN graph and clustering ####
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
#     2.2.3.1.9 UMAP visualization ####
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
#     2.2.3.1.10 Markers of each cluster ####
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene), "Protein Coding", "Transposable Element"))
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

#     2.2.3.1.11 Annotation ####
refskin = "/datos/sensence/emilio/VD/input/ref/lin"
options(Seurat.object.assay.version = "v3")
pbmc.skin <- Read10X(refskin)
pbmc.skin <- CreateSeuratObject(counts = pbmc.skin, min.cells = 3, min.features = 200)
first_colnames <- c("barcodes", "first_level_ann")
first_level <- vroom("https://github.com/kasperlab/Joost_et_al_2020_Cell_Stem_Cell/raw/master/celltypes_1st_level.txt", 
                     col_names = first_colnames)
second_colnames <- c("barcodes", "second_level_ann")
second_level <- vroom("https://github.com/kasperlab/Joost_et_al_2020_Cell_Stem_Cell/raw/master/celltypes_2nd_level.txt",
                      col_names = second_colnames,
                      delim = "\t")
second_cellnames <- c("second_level_ann", "cell_type")
second_names <- vroom("https://github.com/kasperlab/Joost_et_al_2020_Cell_Stem_Cell/raw/master/cell_types_cluster_names.txt",
                      col_names = second_cellnames,
                      delim = ":")
for (col in names(second_names)) {
  second_names[[col]] <- gsub("'", "", second_names[[col]])
  second_names[[col]] <- gsub(",", "", second_names[[col]])
}
cell_ann <- merge(first_level, second_level, by.x="barcodes")
cell_ann <- merge(cell_ann, second_names, by.x="second_level_ann")
cell_ann <- cell_ann[,2:4]
cell_ann <- subset(cell_ann, !grepl("5w$", barcodes))
cell_ann$barcodes <- sub("-9w$", "-1", cell_ann$barcodes)
pbmc.skin@meta.data$barcodes <- rownames(pbmc.skin@meta.data)
merged_data <- merge(x = pbmc.skin@meta.data, y = cell_ann, by.x = "barcodes", by.y = "barcodes", all.x = TRUE)
pbmc.skin[["cell.type"]] <- merged_data$cell_type
rows_with_na <- merged_data[apply(merged_data, 1, function(row) any(is.na(row))), ]
rows_with_na <- rows_with_na[,1]
pbmc.skin <- pbmc.skin[,!colnames(pbmc.skin) %in% rows_with_na]
rm(cell_ann, first_level, merged_data, second_level, second_names, col, first_colnames, refskin, rows_with_na,second_cellnames, second_colnames)

ref <- pbmc.skin@assays[["RNA"]]@counts
test <- pbmc.seurat@assays[["RNA"]]@counts
test_singleR <- SingleR( test = test, ref = ref, labels = pbmc.skin@meta.data[["cell.type"]], num.threads = 60)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

#     2.2.3.1.12 Explore markers ####
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

#     2.2.3.1.13 Save ####
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))
#     2.2.3.1.14 clean ####
rm(list=ls())
gc()
#    2.2.3.2 lin_vehicle ####
#     2.2.3.2.1 Values ####
input.seurat = "/datos/sensence/emilio/VD/input/sc_data/lin/vehicle"
sample.name = "lin_vehicle"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/lin/vehicle/"

#     2.2.3.2.2 load ####
control_path = "/datos/sensence/emilio/VD/input/sc_data/lin/vehicle/GSM5266942_C5_matrix.tsv"
data_control <- read.delim(control_path, header = TRUE, sep = '\t', row.names = "X")
counts_matrix = Matrix(as.matrix(data_control) , sparse=TRUE)
options(Seurat.object.assay.version = "v3")
pbmc.seurat <- CreateSeuratObject(counts = counts_matrix, min.cells = 3, min.features = 200, project = "vehicle")
class(pbmc.seurat[["RNA"]])
#     2.2.3.2.3 Quality control ####
pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
png(filename = file.path(plots.dir, paste0(sample.name, "_QC.png")),
    width = 20, height = 16, units = "in", res=300)
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

pbmc.seurat <- subset(pbmc.seurat, subset = nFeature_RNA > 200 & percent.mt < 7.5)

#     2.2.3.2.4 Normalization ####
pbmc.seurat <- NormalizeData(pbmc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#     2.2.3.2.5 Feature selection ####
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
png(filename = file.path(plots.dir, paste0(sample.name, "_variableFeatures_scatter_plot.png")),
    width = 20, height = 16, units = "in", res=300)
print(p2)
dev.off()
rm(top50, top50.colors, p1,p2)

#     2.2.1.1.6 
#     2.2.3.2.6 Data scaling ####
pbmc.seurat <- ScaleData(pbmc.seurat, features = rownames(pbmc.seurat))
#     2.2.3.2.7 Dimensional reduction ####
pbmc.seurat <- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

dplot <- DimPlot(pbmc.seurat, reduction = "pca", group.by = "orig.ident", cols = "darkgoldenrod1") + 
  ggtitle(paste0("PCA on Variable Features\n", sample.name)) +
  theme(legend.position="bottom")

eplot <- ElbowPlot(pbmc.seurat) + ggtitle("Elbow plot")

png(filename = file.path(plots.dir, paste0(sample.name, "_PCA_dim_plot.png")),
    width = 18, height = 9, units = "in", res=300)

print(dplot + eplot)

dev.off()

rm(dplot, eplot)

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

hm.plots <- DimHeatmap(pbmc.seurat, dims = 1:15, reduction = "pca", balanced = TRUE, combine = FALSE,  fast = FALSE)
hm.plots <- lapply(X = hm.plots, FUN = function(x) {
  s <- x$plot_env$feature.order
  s.colors <- rep("black", length(s))
  s.colors[!grepl("hrg", s)] <- "lightslateblue"
  
  s <- gsub(".+hrg", "", s)
  x <- x + scale_y_discrete(labels = s) + theme(axis.text.y = element_text(face = "bold", color = s.colors))
  
  return(x)
})
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

pbmc.seurat <- JackStraw(pbmc.seurat, num.replicate = 100)
pbmc.seurat <- ScoreJackStraw(pbmc.seurat, dims = 1:20)

png(filename = file.path(plots.dir, paste0(sample.name, "_JackStraw_plot.png")),
    width = 20, height = 16, units = "in", res=300)
JackStrawPlot(pbmc.seurat, dims = 1:20) + 
  ggtitle(paste0(sample.name, "  -  JackStraw p-values distribution by principal component"))
dev.off()

#     2.2.3.2.8 SNN graph and clustering ####
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:10)
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 2)
#     2.2.3.2.9 UMAP visualization ####
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
#     2.2.3.2.10 Markers of each cluster ####
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
pbmc.markers <- filter(pbmc.markers, p_val_adj < 0.05)
pbmc.markers <- data.frame(pbmc.markers, gene_type = if_else(grepl("hrg", pbmc.markers$gene), "Protein Coding", "Transposable Element"))
write.csv(x = pbmc.markers, 
          file = file.path(output.dir, paste0(sample.name, "_clusters_markers.csv")))

#     2.2.3.2.11 Annotation ####
refskin = "/datos/sensence/emilio/VD/input/ref/lin"
options(Seurat.object.assay.version = "v3")
pbmc.skin <- Read10X(refskin)
pbmc.skin <- CreateSeuratObject(counts = pbmc.skin, min.cells = 3, min.features = 200)
first_colnames <- c("barcodes", "first_level_ann")
first_level <- vroom("https://github.com/kasperlab/Joost_et_al_2020_Cell_Stem_Cell/raw/master/celltypes_1st_level.txt", 
                     col_names = first_colnames)
second_colnames <- c("barcodes", "second_level_ann")
second_level <- vroom("https://github.com/kasperlab/Joost_et_al_2020_Cell_Stem_Cell/raw/master/celltypes_2nd_level.txt",
                      col_names = second_colnames,
                      delim = "\t")
second_cellnames <- c("second_level_ann", "cell_type")
second_names <- vroom("https://github.com/kasperlab/Joost_et_al_2020_Cell_Stem_Cell/raw/master/cell_types_cluster_names.txt",
                      col_names = second_cellnames,
                      delim = ":")
for (col in names(second_names)) {
  second_names[[col]] <- gsub("'", "", second_names[[col]])
  second_names[[col]] <- gsub(",", "", second_names[[col]])
}
cell_ann <- merge(first_level, second_level, by.x="barcodes")
cell_ann <- merge(cell_ann, second_names, by.x="second_level_ann")
cell_ann <- cell_ann[,2:4]
cell_ann <- subset(cell_ann, !grepl("5w$", barcodes))
cell_ann$barcodes <- sub("-9w$", "-1", cell_ann$barcodes)
pbmc.skin@meta.data$barcodes <- rownames(pbmc.skin@meta.data)
merged_data <- merge(x = pbmc.skin@meta.data, y = cell_ann, by.x = "barcodes", by.y = "barcodes", all.x = TRUE)
pbmc.skin[["cell.type"]] <- merged_data$cell_type
rows_with_na <- merged_data[apply(merged_data, 1, function(row) any(is.na(row))), ]
rows_with_na <- rows_with_na[,1]
pbmc.skin <- pbmc.skin[,!colnames(pbmc.skin) %in% rows_with_na]
rm(cell_ann, first_level, merged_data, second_level, second_names, col, first_colnames, refskin, rows_with_na,second_cellnames, second_colnames)

ref <- pbmc.skin@assays[["RNA"]]@counts
test <- pbmc.seurat@assays[["RNA"]]@counts
test_singleR <- SingleR( test = test, ref = ref, labels = pbmc.skin@meta.data[["cell.type"]], num.threads = 60)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
rm(barcodes_annotation, dge, gene_list, ref, test, test_singleR)
pct.p <- DimPlot(pbmc.seurat, label = F , group.by = "SingleR.labels", repel = T, label.size = 7,
                 cols = DiscretePalette(length(unique(pbmc.seurat$SingleR.labels)))) +
  theme(legend.position="bottom")+
  ggtitle(sample.name)

png(filename = file.path(plots.dir, paste0(sample.name, "_cluster by cell type_UMAP.png")),
    width = 19, height = 9, units = "in", res=300)

print(pct.p)

dev.off()

#     2.2.3.2.12 Explore markers ####
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

#     2.2.3.2.13 Save ####
saveRDS(object = pbmc.seurat, file = file.path(output.dir, paste0(sample.name, ".rds")))
#     2.2.2.1.14 clean ####
rm(list=ls())
gc()

#   2.3 Integration ####
#    2.3.1 aem ####
#     2.3.1.1 libraries ####

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)

#     2.3.1.2 Values ####

input.seurat.vd = "/datos/sensence/emilio/VD/output/seurat/prepros/aem_vd.rds"
input.seurat.vehicle = "/datos/sensence/emilio/VD/output/seurat/prepros/aem_vehicle.rds"
sample.name = "aem_combined"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/aem/"
set.seed(12345)

#     2.3.1.3 load ########

vehicle <- readRDS(input.seurat.vehicle)
vd <- readRDS(input.seurat.vd)

#     2.3.1.4 Set idents #######

Idents(vehicle) <- 'vehicle'
df <- data.frame(Idents(vehicle))
vehicle <- AddMetaData(object = vehicle, metadata = df, col.name = 'experiment')
Idents(vehicle)

Idents(vd) <- 'vd'
df <- data.frame(Idents(vd))
vd <- AddMetaData(object = vd, metadata = df, col.name = 'experiment')
Idents(vd)


#     2.3.1.5 Merge/integration/prepros ####

combined <- merge(vehicle, y = vd, project = "combined")
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$experiment)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- IntegrateLayers(object = combined, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                            verbose = FALSE)

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 1)

#     2.3.1.6 Plotting individually #####
combined <- RunUMAP(combined, dims = 1:30, reduction = "harmony")

DimPlot(combined, reduction = "umap", group.by = c("experiment", "SingleR.labels"))
umap <- DimPlot(combined, reduction = "umap", split.by = "experiment", group.by = ("SingleR.labels"))
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP_annotated_integrated.png")),
    width = 20, height = 16, units = "in", res=300)
umap + labs(title= paste0(sample.name, "_UMAP_annnotated"))
dev.off()

#     2.3.1.7 Save #####

saveRDS(object = combined, file = file.path(output.dir, paste0(sample.name, ".rds")))


#     2.3.1.8 Clean ####
rm(list=ls())
gc()

#    2.3.1 hanai ####
#     2.3.1.2 Values ####

input.seurat.vd = "/datos/sensence/emilio/VD/output/seurat/prepros/hanai_vd.rds"
input.seurat.vehicle = "/datos/sensence/emilio/VD/output/seurat/prepros/hanai_vehicle.rds"
sample.name = "hanai_combined"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/hanai/"
set.seed(12345)

#     2.3.1.3 load ########

vehicle <- readRDS(input.seurat.vehicle)
vd <- readRDS(input.seurat.vd)

# Here we will start the purging of the immune cells since they are from bone marrow and adds much noise 
vehicle <- subset(vehicle, SingleR.labels != 'Immune cells')
saveRDS(vehicle, "/datos/sensence/emilio/VD/output/seurat/hanai_control.rds")
vd <- subset (vd, SingleR.labels != 'Immune cells')
saveRDS(vd, "/datos/sensence/emilio/VD/output/seurat/hanai_vd.rds")

#     2.3.1.4 Set idents #######

Idents(vehicle) <- 'vehicle'
df <- data.frame(Idents(vehicle))
vehicle <- AddMetaData(object = vehicle, metadata = df, col.name = 'experiment')
Idents(vehicle)

Idents(vd) <- 'vd'
df <- data.frame(Idents(vd))
vd <- AddMetaData(object = vd, metadata = df, col.name = 'experiment')
Idents(vd)


#     2.3.1.5 Merge/integration/prepros ####

combined <- merge(vehicle, y = vd, project = "combined")
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$experiment)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- IntegrateLayers(object = combined, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                            verbose = FALSE)

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 1)

#     2.3.1.6 Plotting individually #####
combined <- RunUMAP(combined, dims = 1:30, reduction = "harmony")

DimPlot(combined, reduction = "umap", group.by = c("experiment", "SingleR.labels"))
umap <- DimPlot(combined, reduction = "umap", split.by = "experiment", group.by = ("SingleR.labels"))
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP_annotated_integrated.png")),
    width = 20, height = 16, units = "in", res=300)
umap + labs(title= paste0(sample.name, "_UMAP_annnotated"))
dev.off()

#     2.3.1.7 Save #####

saveRDS(object = combined, file = file.path(output.dir, paste0(sample.name, ".rds")))



#     2.3.1.8 Clean ####
rm(list=ls())
gc()

#    2.3.2 lin ####
#     2.3.2.2 Values ####

input.seurat.vd = "/datos/sensence/emilio/VD/output/seurat/prepros/lin_vd.rds"
input.seurat.vehicle = "/datos/sensence/emilio/VD/output/seurat/prepros/lin_vehicle.rds"
sample.name = "lin_combined"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/lin/"
set.seed(12345)

#     2.3.2.3 load ########

vehicle <- readRDS(input.seurat.vehicle)
vd <- readRDS(input.seurat.vd)

#     2.3.2.4 Set idents #######

Idents(vehicle) <- 'vehicle'
df <- data.frame(Idents(vehicle))
vehicle <- AddMetaData(object = vehicle, metadata = df, col.name = 'experiment')
Idents(vehicle)

Idents(vd) <- 'vd'
df <- data.frame(Idents(vd))
vd <- AddMetaData(object = vd, metadata = df, col.name = 'experiment')
Idents(vd)


#     2.3.2.5 Merge/integration/prepros ####

combined <- merge(vehicle, y = vd, project = "combined")
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$experiment)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- IntegrateLayers(object = combined, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                            verbose = FALSE)

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 1)

#     2.3.2.6 Plotting individually #####
combined <- RunUMAP(combined, dims = 1:30, reduction = "harmony")

DimPlot(combined, reduction = "umap", group.by = c("experiment", "SingleR.labels"))
umap <- DimPlot(combined, reduction = "umap", split.by = "experiment", group.by = ("SingleR.labels"))
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP_annotated_integrated.png")),
    width = 20, height = 16, units = "in", res=300)
umap + labs(title= paste0(sample.name, "_UMAP_annnotated"))
dev.off()

#     2.3.2.7 Save #####

saveRDS(object = combined, file = file.path(output.dir, paste0(sample.name, ".rds")))



#     2.3.2.8 Clean ####
rm(list=ls())
gc()

# 3. Enrichment ####

#  3.1 libraries ####
library(Seurat)
library(GSEABase)
library(escape)
library(AUCell)
library(tidyverse)
library(dittoSeq)
library(gridExtra)
library(patchwork)
library(biomaRt)
library(ggrepel)
library(ggstatsplot)
library(rstatix)
library(BiocParallel)
set.seed(12345)
#  3.2 aem ####
#   3.2.1 Set values ##########

input.seurat = "/datos/sensence/emilio/VD/output/seurat/prepros/aem_combined.rds"
sample.name = "enrich_aem"
output.dir = "/datos/sensence/emilio/VD/output/seurat/enrichment/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/enrichment/aem/"

#   3.2.2 import datasets ####

### Senescence scores
cellage <- read_csv("input/scores/cellage.csv")

genage <- read_tsv("input/scores/geneage.tsv")

senmayo_mouse <- read_csv("input/scores/senmayo_mouse.csv")

#### Seurat data
aem_combined <- readRDS(input.seurat)

#   3.2.3 Clean data  #####

#### Convert to respective species

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gs_cellage <- getLDS(attributes = c("external_gene_name"),
                     filters = "external_gene_name", values = unique(cellage$`Gene symbol`), mart = human,
                     attributesL = c("external_gene_name"), martL = mouse)
rm(human,mouse, cellage)

#### Prepros integration

aem_combined <- NormalizeData(aem_combined)
aem_combined <- FindVariableFeatures(aem_combined)
aem_combined <- ScaleData(aem_combined)
aem_combined <- RunPCA(aem_combined)

aem_combined <- FindNeighbors(aem_combined, reduction = "pca", dims = 1:30)
aem_combined <- FindClusters(aem_combined, resolution = 1)
aem_combined <- RunUMAP(aem_combined, dims = 1:30, reduction = "pca")

#   3.2.4 Gene sets #####

gs_cellage <- gs_cellage[,2]
gs_geneage <- unique(genage$`Gene Symbol`)
gs_senmayo <- unique(senmayo_mouse$`Gene(murine)`)

gs_cellage <- GeneSet(gs_cellage)
gs_geneage <- GeneSet(gs_geneage)
gs_senmayo <- GeneSet(gs_senmayo)

geneSets <- list(geneSet1=gs_cellage, geneSet2=gs_geneage, geneSet3=gs_senmayo)

geneSets$geneSet1@setName <- "Cellage"
geneSets$geneSet2@setName <- "GenAge"
geneSets$geneSet3@setName <- "SENmayo"
geneSets <- GeneSetCollection(geneSets)

#   3.2.5 Enrichment ####

ES <- escape.matrix(aem_combined, gene.sets = geneSets, method = "AUCell", groups = 1000, BPPARAM = SnowParam(workers = 40))
aem_combined <- Seurat::AddMetaData(aem_combined, ES)

#   3.2.6 Plot #####

### Vln plot
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot SENMayo by subtype and experiment.png")),
    width = 20, height = 16, units = "in", res=300)
vlnplot_aem <- VlnPlot(
  object = aem_combined,
  features = "SENmayo",
  group.by = "SingleR.labels",
  pt.size = 0.5,
  log = FALSE,
  split.by = "experiment",
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  ncol = NULL,
  slot = "data",
  split.plot = FALSE,
  stack = FALSE,
  combine = TRUE,
  fill.by = "feature",
  flip = FALSE,
  add.noise = TRUE
) + ggtitle(paste0(sample.name,"Violin_plot_SENmayo_by_subtype_and_experiment"))
plot(vlnplot_aem)
dev.off()

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$cellage <- aem_combined@meta.data$Cellage
umap_data$geneage <- aem_combined@meta.data$GenAge
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment
label_counts <- table(umap_data$SingleR.labels)
valid_labels <- names(label_counts[label_counts >= 80])
umap_data <- umap_data %>%
  filter(SingleR.labels %in% valid_labels)

#see the overall violin plot 
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot overrall SENmayo.png")),
    width = 20, height = 16, units = "in", res=300)
overal_plot <- ggbetweenstats(
  data  = umap_data,
  x     = experiment,
  y     = SENmayo,
  title = paste0(sample.name, "__overall_SENmayo score values"),
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 7, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
  #pairwise.display = "none"
)
overal_plot + 
  theme(axis.text.x=element_text(size = 17.5))  + theme(plot.title=element_text(size=25))
dev.off()
### Significance 

significance <- data.frame(
  Variable = character(),
  Wilcoxon_Result = character(),
  P_Value = numeric(),
  Median_Vehicle = numeric(),
  Median_VD = numeric(),
  stringsAsFactors = FALSE
)

# List of variables to compare
variables_to_compare <- c("Cellage", "GenAge", "SENmayo")

# Iterate through each variable
for (variable in variables_to_compare) {
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(get(variable) ~ experiment, data = aem_combined@meta.data)
  
  # Check if the Wilcoxon statistic is numeric
  wilcox_statistic <- ifelse(is.numeric(wilcox_result$statistic), as.character(wilcox_result$statistic), NA)
  
  # Calculate median for each group
  median_vehicle <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vehicle", variable], na.rm = TRUE)
  median_vd <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vd", variable], na.rm = TRUE)
  
  # Create a row for the result dataframe
  result_row <- data.frame(
    Variable = variable,
    Wilcoxon_Result = wilcox_statistic,
    P_Value = wilcox_result$p.value,
    Median_Vehicle = median_vehicle,
    Median_VD = median_vd,
    stringsAsFactors = FALSE
  )
  
  # Add the row to the result dataframe
  significance <- bind_rows(significance, result_row)
}

# Print the resulting dataframe
print(significance)
write.csv(x = significance, 
          file = file.path(output.dir, paste0(sample.name, "_significance.csv")))

### UMAP

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment

# Choose the specific labels you want to label
labels_to_label <- unique(umap_data$SingleR.labels)

# Create UMAP plot using ggplot2 with labels for each unique cell subtype and with certain experiment
base_vehicle <- ggplot(umap_data %>% filter(experiment == "vehicle"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )   +
  ggtitle("AEM_Vehicle")

# Display the plot
print(base_vehicle)

#again but with different experiment
base_vd <- ggplot(umap_data %>% filter(experiment == "vd"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )  +
  ggtitle("AEM_VD")

# Display the plot
print(base_vd)

#print them together
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP SENmayo experiments.png")),
    width = 20, height = 16, units = "in", res=300)
print((base_vehicle | base_vd) +
        plot_annotation(title = paste0("\n",sample.name,"\n"),
                        theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)))
)
dev.off()

#### Pairwise wilcoxon top 3 cell subtypes

pairwise.wilcox.test(umap_data$SENmayo, umap_data$SingleR.labels,
                     p.adjust.method = "BH")
# Subset the data for "vehicle" and "vd" experiments
vehicle_data <- aem_combined@meta.data %>%
  filter(experiment == "vehicle")

vd_data <- aem_combined@meta.data %>%
  filter(experiment == "vd")

# Get unique labels
unique_labels <- unique(aem_combined@meta.data$SingleR.labels)
# Initialize a data frame to store the results
wilcox_results <- data.frame(
  SingleR_label = character(),
  wilcox_value = numeric(),
  p_value = numeric(),
  experiment_comparison = character(),
  VD_minus_vehicle = numeric(), # New column for the mean difference
  stringsAsFactors = FALSE
)

# Perform Wilcoxon tests for each label
for (label in unique_labels) {
  # Subset data for the current label
  vehicle_label_data <- vehicle_data %>%
    filter(SingleR.labels == label)
  
  vd_label_data <- vd_data %>%
    filter(SingleR.labels == label)
  
  # Check if both datasets have enough observations
  if (nrow(vehicle_label_data) >= 3 && nrow(vd_label_data) >= 3) {
    # Calculate means for each group
    vehicle_mean <- mean(vehicle_label_data$SENmayo)
    vd_mean <- mean(vd_label_data$SENmayo)
    
    # Calculate the difference in means
    mean_difference <- vd_mean - vehicle_mean
    
    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(vehicle_label_data$SENmayo, vd_label_data$SENmayo)
    
    # Store the results along with the calculated mean difference
    wilcox_results <- rbind(
      wilcox_results,
      data.frame(
        SingleR_label = label,
        wilcox_value = wilcox_result$statistic,
        p_value = wilcox_result$p.value,
        experiment_comparison = "vehicle vs vd",
        VD_minus_vehicle = mean_difference,  # Store the mean difference
        stringsAsFactors = FALSE
      )
    )
  } else {
    # Print a message or handle cases where there are not enough observations
    cat("Skipping label", label, "due to insufficient observations.\n")
  }
}

# Order the results by wilcox_value
wilcox_results <- wilcox_results %>% arrange(p_value)

# Get the top 3 results
top_results <- head(wilcox_results, 3)

# Print the top results
print(top_results)
write.csv(x = top_results, 
          file = file.path(output.dir, paste0(sample.name, "_top_3_wilcox_cellsubtypes.csv")))

#   3.2.7 Save ####

saveRDS(object = aem_combined, file = paste0(output.dir,sample.name, "_ES.rds"))

#   3.2.8 clean ####
rm(list=ls())
gc()

#  3.3 hanai ####
#   3.3.1 Set values ##########

input.seurat = "/datos/sensence/emilio/VD/output/seurat/prepros/hanai_combined.rds"
sample.name = "enrich_hanai"
output.dir = "/datos/sensence/emilio/VD/output/seurat/enrichment/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/enrichment/hanai/"

#   3.3.2 import datasets ####

### Senescence scores
cellage <- read_csv("input/scores/cellage.csv")

genage <- read_tsv("input/scores/geneage.tsv")

senmayo_mouse <- read_csv("input/scores/senmayo_mouse.csv")

#### Seurat data
aem_combined <- readRDS(input.seurat)

#   3.3.3 Clean data  #####

#### Convert to respective species

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gs_cellage <- getLDS(attributes = c("external_gene_name"),
                     filters = "external_gene_name", values = unique(cellage$`Gene symbol`), mart = human,
                     attributesL = c("external_gene_name"), martL = mouse)
rm(human,mouse, cellage)

#### Prepros integration

aem_combined <- NormalizeData(aem_combined)
aem_combined <- FindVariableFeatures(aem_combined)
aem_combined <- ScaleData(aem_combined)
aem_combined <- RunPCA(aem_combined)

aem_combined <- FindNeighbors(aem_combined, reduction = "pca", dims = 1:30)
aem_combined <- FindClusters(aem_combined, resolution = 1)
aem_combined <- RunUMAP(aem_combined, dims = 1:30, reduction = "pca")

#   3.3.4 Gene sets #####

gs_cellage <- gs_cellage[,2]
gs_geneage <- unique(genage$`Gene Symbol`)
gs_senmayo <- unique(senmayo_mouse$`Gene(murine)`)

gs_cellage <- GeneSet(gs_cellage)
gs_geneage <- GeneSet(gs_geneage)
gs_senmayo <- GeneSet(gs_senmayo)

geneSets <- list(geneSet1=gs_cellage, geneSet2=gs_geneage, geneSet3=gs_senmayo)

geneSets$geneSet1@setName <- "Cellage"
geneSets$geneSet2@setName <- "GenAge"
geneSets$geneSet3@setName <- "SENmayo"
geneSets <- GeneSetCollection(geneSets)

#   3.3.5 Enrichment ####

ES <- escape.matrix(aem_combined, gene.sets = geneSets, method = "AUCell", groups = 1000, BPPARAM = SnowParam(workers = 40))
aem_combined <- Seurat::AddMetaData(aem_combined, ES)

#   3.3.6 Plot #####

### Vln plot
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot SENMayo by subtype and experiment.png")),
    width = 20, height = 16, units = "in", res=300)
vlnplot_aem <- VlnPlot(
  object = aem_combined,
  features = "SENmayo",
  group.by = "SingleR.labels",
  pt.size = 0.5,
  log = FALSE,
  split.by = "experiment",
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  ncol = NULL,
  slot = "data",
  split.plot = FALSE,
  stack = FALSE,
  combine = TRUE,
  fill.by = "feature",
  flip = FALSE,
  add.noise = TRUE
) + ggtitle(paste0(sample.name,"Violin_plot_SENmayo_by_subtype_and_experiment"))
plot(vlnplot_aem)
dev.off()

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$cellage <- aem_combined@meta.data$Cellage
umap_data$geneage <- aem_combined@meta.data$GenAge
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment
label_counts <- table(umap_data$SingleR.labels)
valid_labels <- names(label_counts[label_counts >= 80])
umap_data <- umap_data %>%
  filter(SingleR.labels %in% valid_labels)

#see the overall violin plot 
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot overrall SENmayo.png")),
    width = 20, height = 16, units = "in", res=300)
overal_plot <- ggbetweenstats(
  data  = umap_data,
  x     = experiment,
  y     = SENmayo,
  title = paste0(sample.name, "__overall_SENmayo score values"),
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 7, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
  #pairwise.display = "none"
)
overal_plot + 
  theme(axis.text.x=element_text(size = 17.5))  + theme(plot.title=element_text(size=25))
dev.off()
### Significance 

significance <- data.frame(
  Variable = character(),
  Wilcoxon_Result = character(),
  P_Value = numeric(),
  Median_Vehicle = numeric(),
  Median_VD = numeric(),
  stringsAsFactors = FALSE
)

# List of variables to compare
variables_to_compare <- c("Cellage", "GenAge", "SENmayo")

# Iterate through each variable
for (variable in variables_to_compare) {
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(get(variable) ~ experiment, data = aem_combined@meta.data)
  
  # Check if the Wilcoxon statistic is numeric
  wilcox_statistic <- ifelse(is.numeric(wilcox_result$statistic), as.character(wilcox_result$statistic), NA)
  
  # Calculate median for each group
  median_vehicle <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vehicle", variable], na.rm = TRUE)
  median_vd <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vd", variable], na.rm = TRUE)
  
  # Create a row for the result dataframe
  result_row <- data.frame(
    Variable = variable,
    Wilcoxon_Result = wilcox_statistic,
    P_Value = wilcox_result$p.value,
    Median_Vehicle = median_vehicle,
    Median_VD = median_vd,
    stringsAsFactors = FALSE
  )
  
  # Add the row to the result dataframe
  significance <- bind_rows(significance, result_row)
}

# Print the resulting dataframe
print(significance)
write.csv(x = significance, 
          file = file.path(output.dir, paste0(sample.name, "_significance.csv")))

### UMAP

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment

# Choose the specific labels you want to label
labels_to_label <- unique(umap_data$SingleR.labels)

# Create UMAP plot using ggplot2 with labels for each unique cell subtype and with certain experiment
base_vehicle <- ggplot(umap_data %>% filter(experiment == "vehicle"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )   +
  ggtitle("AEM_Vehicle")

# Display the plot
print(base_vehicle)

#again but with different experiment
base_vd <- ggplot(umap_data %>% filter(experiment == "vd"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )  +
  ggtitle("AEM_VD")

# Display the plot
print(base_vd)

#print them together
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP SENmayo experiments.png")),
    width = 20, height = 16, units = "in", res=300)
print((base_vehicle | base_vd) +
        plot_annotation(title = paste0("\n",sample.name,"\n"),
                        theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)))
)
dev.off()

#### Pairwise wilcoxon top 3 cell subtypes

pairwise.wilcox.test(umap_data$SENmayo, umap_data$SingleR.labels,
                     p.adjust.method = "BH")
# Subset the data for "vehicle" and "vd" experiments
vehicle_data <- aem_combined@meta.data %>%
  filter(experiment == "vehicle")

vd_data <- aem_combined@meta.data %>%
  filter(experiment == "vd")

# Get unique labels
unique_labels <- unique(aem_combined@meta.data$SingleR.labels)
# Initialize a data frame to store the results
wilcox_results <- data.frame(
  SingleR_label = character(),
  wilcox_value = numeric(),
  p_value = numeric(),
  experiment_comparison = character(),
  VD_minus_vehicle = numeric(), # New column for the mean difference
  stringsAsFactors = FALSE
)

# Perform Wilcoxon tests for each label
for (label in unique_labels) {
  # Subset data for the current label
  vehicle_label_data <- vehicle_data %>%
    filter(SingleR.labels == label)
  
  vd_label_data <- vd_data %>%
    filter(SingleR.labels == label)
  
  # Check if both datasets have enough observations
  if (nrow(vehicle_label_data) >= 3 && nrow(vd_label_data) >= 3) {
    # Calculate means for each group
    vehicle_mean <- mean(vehicle_label_data$SENmayo)
    vd_mean <- mean(vd_label_data$SENmayo)
    
    # Calculate the difference in means
    mean_difference <- vd_mean - vehicle_mean
    
    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(vehicle_label_data$SENmayo, vd_label_data$SENmayo)
    
    # Store the results along with the calculated mean difference
    wilcox_results <- rbind(
      wilcox_results,
      data.frame(
        SingleR_label = label,
        wilcox_value = wilcox_result$statistic,
        p_value = wilcox_result$p.value,
        experiment_comparison = "vehicle vs vd",
        VD_minus_vehicle = mean_difference,  # Store the mean difference
        stringsAsFactors = FALSE
      )
    )
  } else {
    # Print a message or handle cases where there are not enough observations
    cat("Skipping label", label, "due to insufficient observations.\n")
  }
}

# Order the results by wilcox_value
wilcox_results <- wilcox_results %>% arrange(p_value)

# Get the top 3 results
top_results <- head(wilcox_results, 3)

# Print the top results
print(top_results)
write.csv(x = top_results, 
          file = file.path(output.dir, paste0(sample.name, "_top_3_wilcox_cellsubtypes.csv")))

#   3.3.7 Save ####

saveRDS(object = aem_combined, file = paste0(output.dir,sample.name, "_ES.rds"))

#   3.3.8 clean ####
rm(list=ls())
gc()

#  3.4 lin ####
#   3.4.1 Set values ##########

input.seurat = "/datos/sensence/emilio/VD/output/seurat/prepros/lin_combined.rds"
sample.name = "enrich_lin"
output.dir = "/datos/sensence/emilio/VD/output/seurat/enrichment/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/enrichment/lin/"

#   3.4.2 import datasets ####

### Senescence scores
cellage <- read_csv("input/scores/cellage.csv")

genage <- read_tsv("input/scores/geneage.tsv")

senmayo_mouse <- read_csv("input/scores/senmayo_mouse.csv")

#### Seurat data
aem_combined <- readRDS(input.seurat)

#   3.4.3 Clean data  #####

#### Convert to respective species

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gs_cellage <- getLDS(attributes = c("external_gene_name"),
                     filters = "external_gene_name", values = unique(cellage$`Gene symbol`), mart = human,
                     attributesL = c("external_gene_name"), martL = mouse)
rm(human,mouse, cellage)

#### Prepros integration

aem_combined <- NormalizeData(aem_combined)
aem_combined <- FindVariableFeatures(aem_combined)
aem_combined <- ScaleData(aem_combined)
aem_combined <- RunPCA(aem_combined)

aem_combined <- FindNeighbors(aem_combined, reduction = "pca", dims = 1:30)
aem_combined <- FindClusters(aem_combined, resolution = 1)
aem_combined <- RunUMAP(aem_combined, dims = 1:30, reduction = "pca")

#   3.4.4 Gene sets #####

gs_cellage <- gs_cellage[,2]
gs_geneage <- unique(genage$`Gene Symbol`)
gs_senmayo <- unique(senmayo_mouse$`Gene(murine)`)

gs_cellage <- GeneSet(gs_cellage)
gs_geneage <- GeneSet(gs_geneage)
gs_senmayo <- GeneSet(gs_senmayo)

geneSets <- list(geneSet1=gs_cellage, geneSet2=gs_geneage, geneSet3=gs_senmayo)

geneSets$geneSet1@setName <- "Cellage"
geneSets$geneSet2@setName <- "GenAge"
geneSets$geneSet3@setName <- "SENmayo"
geneSets <- GeneSetCollection(geneSets)

#   3.4.5 Enrichment ####

ES <- escape.matrix(aem_combined, gene.sets = geneSets, method = "AUCell", groups = 1000, BPPARAM = SnowParam(workers = 40))
aem_combined <- Seurat::AddMetaData(aem_combined, ES)

#   3.4.6 Plot #####

### Vln plot
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot SENMayo by subtype and experiment.png")),
    width = 20, height = 16, units = "in", res=300)
vlnplot_aem <- VlnPlot(
  object = aem_combined,
  features = "SENmayo",
  group.by = "SingleR.labels",
  pt.size = 0.5,
  log = FALSE,
  split.by = "experiment",
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  ncol = NULL,
  slot = "data",
  split.plot = FALSE,
  stack = FALSE,
  combine = TRUE,
  fill.by = "feature",
  flip = FALSE,
  add.noise = TRUE
) + ggtitle(paste0(sample.name,"Violin_plot_SENmayo_by_subtype_and_experiment"))
plot(vlnplot_aem)
dev.off()

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$cellage <- aem_combined@meta.data$Cellage
umap_data$geneage <- aem_combined@meta.data$GenAge
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment
label_counts <- table(umap_data$SingleR.labels)
valid_labels <- names(label_counts[label_counts >= 80])
umap_data <- umap_data %>%
  filter(SingleR.labels %in% valid_labels)

#see the overall violin plot 
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot overrall SENmayo.png")),
    width = 20, height = 16, units = "in", res=300)
overal_plot <- ggbetweenstats(
  data  = umap_data,
  x     = experiment,
  y     = SENmayo,
  title = paste0(sample.name, "__overall_SENmayo score values"),
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 7, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
  #pairwise.display = "none"
)
overal_plot + 
  theme(axis.text.x=element_text(size = 17.5))  + theme(plot.title=element_text(size=25))
dev.off()
### Significance 

significance <- data.frame(
  Variable = character(),
  Wilcoxon_Result = character(),
  P_Value = numeric(),
  Median_Vehicle = numeric(),
  Median_VD = numeric(),
  stringsAsFactors = FALSE
)

# List of variables to compare
variables_to_compare <- c("Cellage", "GenAge", "SENmayo")

# Iterate through each variable
for (variable in variables_to_compare) {
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(get(variable) ~ experiment, data = aem_combined@meta.data)
  
  # Check if the Wilcoxon statistic is numeric
  wilcox_statistic <- ifelse(is.numeric(wilcox_result$statistic), as.character(wilcox_result$statistic), NA)
  
  # Calculate median for each group
  median_vehicle <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vehicle", variable], na.rm = TRUE)
  median_vd <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vd", variable], na.rm = TRUE)
  
  # Create a row for the result dataframe
  result_row <- data.frame(
    Variable = variable,
    Wilcoxon_Result = wilcox_statistic,
    P_Value = wilcox_result$p.value,
    Median_Vehicle = median_vehicle,
    Median_VD = median_vd,
    stringsAsFactors = FALSE
  )
  
  # Add the row to the result dataframe
  significance <- bind_rows(significance, result_row)
}

# Print the resulting dataframe
print(significance)
write.csv(x = significance, 
          file = file.path(output.dir, paste0(sample.name, "_significance.csv")))

### UMAP

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment

# Choose the specific labels you want to label
labels_to_label <- unique(umap_data$SingleR.labels)

# Create UMAP plot using ggplot2 with labels for each unique cell subtype and with certain experiment
base_vehicle <- ggplot(umap_data %>% filter(experiment == "vehicle"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )   +
  ggtitle("AEM_Vehicle")

# Display the plot
print(base_vehicle)

#again but with different experiment
base_vd <- ggplot(umap_data %>% filter(experiment == "vd"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )  +
  ggtitle("AEM_VD")

# Display the plot
print(base_vd)

#print them together
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP SENmayo experiments.png")),
    width = 20, height = 16, units = "in", res=300)
print((base_vehicle | base_vd) +
        plot_annotation(title = paste0("\n",sample.name,"\n"),
                        theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)))
)
dev.off()

#### Pairwise wilcoxon top 3 cell subtypes

pairwise.wilcox.test(umap_data$SENmayo, umap_data$SingleR.labels,
                     p.adjust.method = "BH")
# Subset the data for "vehicle" and "vd" experiments
vehicle_data <- aem_combined@meta.data %>%
  filter(experiment == "vehicle")

vd_data <- aem_combined@meta.data %>%
  filter(experiment == "vd")

# Get unique labels
unique_labels <- unique(aem_combined@meta.data$SingleR.labels)
# Initialize a data frame to store the results
wilcox_results <- data.frame(
  SingleR_label = character(),
  wilcox_value = numeric(),
  p_value = numeric(),
  experiment_comparison = character(),
  VD_minus_vehicle = numeric(), # New column for the mean difference
  stringsAsFactors = FALSE
)

# Perform Wilcoxon tests for each label
for (label in unique_labels) {
  # Subset data for the current label
  vehicle_label_data <- vehicle_data %>%
    filter(SingleR.labels == label)
  
  vd_label_data <- vd_data %>%
    filter(SingleR.labels == label)
  
  # Check if both datasets have enough observations
  if (nrow(vehicle_label_data) >= 3 && nrow(vd_label_data) >= 3) {
    # Calculate means for each group
    vehicle_mean <- mean(vehicle_label_data$SENmayo)
    vd_mean <- mean(vd_label_data$SENmayo)
    
    # Calculate the difference in means
    mean_difference <- vd_mean - vehicle_mean
    
    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(vehicle_label_data$SENmayo, vd_label_data$SENmayo)
    
    # Store the results along with the calculated mean difference
    wilcox_results <- rbind(
      wilcox_results,
      data.frame(
        SingleR_label = label,
        wilcox_value = wilcox_result$statistic,
        p_value = wilcox_result$p.value,
        experiment_comparison = "vehicle vs vd",
        VD_minus_vehicle = mean_difference,  # Store the mean difference
        stringsAsFactors = FALSE
      )
    )
  } else {
    # Print a message or handle cases where there are not enough observations
    cat("Skipping label", label, "due to insufficient observations.\n")
  }
}

# Order the results by wilcox_value
wilcox_results <- wilcox_results %>% arrange(p_value)

# Get the top 3 results
top_results <- head(wilcox_results, 3)

# Print the top results
print(top_results)
write.csv(x = top_results, 
          file = file.path(output.dir, paste0(sample.name, "_top_3_wilcox_cellsubtypes.csv")))

#   3.4.7 Save ####

saveRDS(object = aem_combined, file = paste0(output.dir,sample.name, "_ES.rds"))

#   3.4.8 clean ####
rm(list=ls())
gc()

# 4. Gene regulatory networks ####
#  4.1 Libraries ####
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 40)

#  4.2 aem ####
#   4.2.1 values ####
input.seurat = "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_aem_ES.rds"
sample.name = "aem_GRN_pseudobulk"
output.dir = "/datos/sensence/emilio/VD/output/seurat/grn/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/grn/aem/"
#   4.2.2 Load and process for WGCNA #################
seurat_obj <- readRDS(input.seurat)
seurat_obj[["RNA3"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
DefaultAssay(object = seurat_obj) <- "RNA3"
seurat_obj[["RNA"]] <- NULL
RenameAssays(object = seurat_obj, assay.name = "RNA3", new.assay.name = "RNA")
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "pseudobulk" # the name of the hdWGCNA experiment
)

#   4.2.3 Pseudobulk and GRN ##############
datExpr <- ConstructPseudobulk(
  seurat_obj,
  group.by = 'SingleR.labels',
  replicate_col = 'experiment',
  label_col = 'experiment',
  assay = 'RNA3',
  slot = 'counts', # this should always be counts!
  min_reps = 2
)
# compute log2CPM normalization
# You can substitute this with another normalization of your choosing.
cpm <- t(apply(datExpr, 1, function(x){
  y <- (x) / sum(x) * 1000000
  log2(y + 1)
}))

seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = cpm
)

# select the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)
# construct the co-expression network and identify gene modules
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  tom_name=sample.name, 
  overwrite_tom=TRUE,
  mergeCutHeight=0.15
)

PlotDendrogram(seurat_obj, main='pseudobulk dendrogram')
# compute the MEs and kMEs
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

# get MEs from seurat object
MEs <- GetMEs(seurat_obj)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


#   4.2.4 Plots ####

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'SingleR.labels')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient(high='red', low='grey95') + 
  xlab('') + ylab('')


p 

# compute the co-expression network umap 
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.4,
  spread=3,
  supervised=TRUE,
  target_weight=0.3
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) + 
  umap_theme() 

# add the module names to the plot by taking the mean coordinates
centroid_df <- umap_df %>% 
  dplyr::group_by(module) %>%
  dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

p <- p + geom_label(
  data = centroid_df, 
  label=as.character(centroid_df$module), 
  fontface='bold', size=2
)
p

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

#   4.2.5 Save ####
saveRDS(seurat_obj, file = file.path(output.dir, paste0(sample.name,".rds")))
#   4.2.6 clean ####
rm(list=ls())
gc()

#  4.3 hanai ####
#   4.3.1 values ####
input.seurat = "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_hanai_ES.rds"
sample.name = "hanai_GRN_pseudobulk"
output.dir = "/datos/sensence/emilio/VD/output/seurat/grn/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/grn/hanai/"
#   4.3.2 Load and process for WGCNA #################
seurat_obj <- readRDS(input.seurat)
seurat_obj[["RNA3"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
DefaultAssay(object = seurat_obj) <- "RNA3"
seurat_obj[["RNA"]] <- NULL
RenameAssays(object = seurat_obj, assay.name = "RNA3", new.assay.name = "RNA")
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "pseudobulk" # the name of the hdWGCNA experiment
)

#   4.3.3 Pseudobulk and GRN ##############
datExpr <- ConstructPseudobulk(
  seurat_obj,
  group.by = 'SingleR.labels',
  replicate_col = 'experiment',
  label_col = 'experiment',
  assay = 'RNA3',
  slot = 'counts', # this should always be counts!
  min_reps = 2
)
# compute log2CPM normalization
# You can substitute this with another normalization of your choosing.
cpm <- t(apply(datExpr, 1, function(x){
  y <- (x) / sum(x) * 1000000
  log2(y + 1)
}))

seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = cpm
)

# select the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)
# construct the co-expression network and identify gene modules
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  tom_name=sample.name, 
  overwrite_tom=TRUE,
  mergeCutHeight=0.15
)

PlotDendrogram(seurat_obj, main='pseudobulk dendrogram')
# compute the MEs and kMEs
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

# get MEs from seurat object
MEs <- GetMEs(seurat_obj)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


#   4.3.4 Plots ####

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'SingleR.labels')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient(high='red', low='grey95') + 
  xlab('') + ylab('')


p 

# compute the co-expression network umap 
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.4,
  spread=3,
  supervised=TRUE,
  target_weight=0.3
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) + 
  umap_theme() 

# add the module names to the plot by taking the mean coordinates
centroid_df <- umap_df %>% 
  dplyr::group_by(module) %>%
  dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

p <- p + geom_label(
  data = centroid_df, 
  label=as.character(centroid_df$module), 
  fontface='bold', size=2
)
p

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

#   4.3.5 Save ####
saveRDS(seurat_obj, file = file.path(output.dir, paste0(sample.name,".rds")))
#   4.3.6 clean ####
rm(list=ls())
gc()

#  4.4 lin ####
#   4.4.1 values ####
input.seurat = "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_lin_ES.rds"
sample.name = "lin_GRN_pseudobulk"
output.dir = "/datos/sensence/emilio/VD/output/seurat/grn/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/grn/lin/"
#   4.4.2 Load and process for WGCNA #################
seurat_obj <- readRDS(input.seurat)
seurat_obj[["RNA3"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
DefaultAssay(object = seurat_obj) <- "RNA3"
seurat_obj[["RNA"]] <- NULL
RenameAssays(object = seurat_obj, assay.name = "RNA3", new.assay.name = "RNA")
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "pseudobulk" # the name of the hdWGCNA experiment
)

#   4.4.3 Pseudobulk and GRN ##############
datExpr <- ConstructPseudobulk(
  seurat_obj,
  group.by = 'SingleR.labels',
  replicate_col = 'experiment',
  label_col = 'experiment',
  assay = 'RNA3',
  slot = 'counts', # this should always be counts!
  min_reps = 2
)
# compute log2CPM normalization
# You can substitute this with another normalization of your choosing.
cpm <- t(apply(datExpr, 1, function(x){
  y <- (x) / sum(x) * 1000000
  log2(y + 1)
}))

seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = cpm
)

# select the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)
# construct the co-expression network and identify gene modules
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  tom_name=sample.name, 
  overwrite_tom=TRUE,
  mergeCutHeight=0.15
)

PlotDendrogram(seurat_obj, main='pseudobulk dendrogram')
# compute the MEs and kMEs
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

# get MEs from seurat object
MEs <- GetMEs(seurat_obj)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


#   4.4.4 Plots ####

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'SingleR.labels')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient(high='red', low='grey95') + 
  xlab('') + ylab('')


p 

# compute the co-expression network umap 
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.4,
  spread=3,
  supervised=TRUE,
  target_weight=0.3
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) + 
  umap_theme() 

# add the module names to the plot by taking the mean coordinates
centroid_df <- umap_df %>% 
  dplyr::group_by(module) %>%
  dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

p <- p + geom_label(
  data = centroid_df, 
  label=as.character(centroid_df$module), 
  fontface='bold', size=2
)
p

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

#   4.4.5 Save ####
saveRDS(seurat_obj, file = file.path(output.dir, paste0(sample.name,".rds")))
#   4.4.6 clean ####
rm(list=ls())
gc()

# 5. Cellular communication ####
#  5.1 Cellular communication senescence only ####
#   5.1.1 libraries ####
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

#   5.1.2 aem ####
#    5.1.2.1 values ####
input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_aem_ES.rds"
sample.name = "aem_cellchat"
output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat"
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/sconly/aem"

#    5.1.2.2 load ####
## get the seurat, and subset for vehicle and vd
seurat_obj <- readRDS(input.seurat)
seurat_vehicle <- subset(x = seurat_obj, subset = experiment == "vehicle")
seurat_vd <- subset(x = seurat_obj, subset = experiment == "vd")

## subset for only the top 10% of the senescent cells for each group
seurat_vehicle <- subset(x = seurat_vehicle, subset = SENmayo >= quantile(seurat_vehicle@meta.data[["SENmayo"]], probs = c(0.90)))
seurat_vd <- subset(x = seurat_vd, subset = SENmayo >= quantile(seurat_vd@meta.data[["SENmayo"]], probs = c(0.90)))

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

#    5.1.2.3 setup interaction db ####

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

#    5.1.2.4 Cellchat obj prepros ####

#process intervention
cellchat_vd <- subsetData(cellchat_vd) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30) # do parallel
cellchat_vd <- identifyOverExpressedGenes(cellchat_vd)
cellchat_vd <- identifyOverExpressedInteractions(cellchat_vd)
#cellchat_vd <- projectData(cellchat_vd, PPI.mouse) is not supported currently 18-06-24
cellchat_vd <- computeCommunProb(cellchat_vd)
cellchat_vd <- computeCommunProbPathway(cellchat_vd)
df.net_vd <- subsetCommunication(cellchat_vd)
cellchat_vd <- aggregateNet(cellchat_vd)
cellchat_vd <- netAnalysis_computeCentrality(cellchat_vd)
gc()

#process vehicle
cellchat_vehicle <- subsetData(cellchat_vehicle) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30) # do parallel
cellchat_vehicle <- identifyOverExpressedGenes(cellchat_vehicle)
cellchat_vehicle <- identifyOverExpressedInteractions(cellchat_vehicle)
#cellchat_vehicle <- projectData(cellchat_vehicle, PPI.mouse)
cellchat_vehicle <- computeCommunProb(cellchat_vehicle)
cellchat_vehicle <- computeCommunProbPathway(cellchat_vehicle)
df.net_vehicle <- subsetCommunication(cellchat_vehicle)
cellchat_vehicle <- aggregateNet(cellchat_vehicle)
cellchat_vehicle <- netAnalysis_computeCentrality(cellchat_vehicle)
gc()

#merge it boi
object.list <- list(vehicle = cellchat_vehicle, vd = cellchat_vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

rm(seurat_obj, seurat_vd, seurat_vehicle, meta, df.net_vd, df.net_vehicle, CellChatDB.use, CellChatDB)
#    5.1.2.5 Plots ####

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
png(filename = file.path(plots.dir, paste0(sample.name, "_celltype_interactions.png")),
    width = 16, height = 10, units = "cm", res=300)
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
    width = 6.5, height = 7.75, units = "in", res=300)
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

#    5.1.2.6 Save  ####

saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_sconly.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_sconly.rds")))

#    5.1.2.7 clean ####
rm(list=ls())
gc()

#   5.1.3 hanai ####
#    5.1.3.1 values ####
input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_hanai_ES.rds"
sample.name = "hanai_cellchat"
output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat"
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/sconly/hanai"

#    5.1.3.2 load ####
## get the seurat, and subset for vehicle and vd
seurat_obj <- readRDS(input.seurat)
seurat_vehicle <- subset(x = seurat_obj, subset = experiment == "vehicle")
seurat_vd <- subset(x = seurat_obj, subset = experiment == "vd")

## subset for only the top 10% of the senescent cells for each group
seurat_vehicle <- subset(x = seurat_vehicle, subset = SENmayo >= quantile(seurat_vehicle@meta.data[["SENmayo"]], probs = c(0.90)))
seurat_vd <- subset(x = seurat_vd, subset = SENmayo >= quantile(seurat_vd@meta.data[["SENmayo"]], probs = c(0.90)))

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

#    5.1.3.3 setup interaction db ####

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

#    5.1.3.4 Cellchat obj prepros ####

#process intervention
cellchat_vd <- subsetData(cellchat_vd) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30) # do parallel
cellchat_vd <- identifyOverExpressedGenes(cellchat_vd)
cellchat_vd <- identifyOverExpressedInteractions(cellchat_vd)
#cellchat_vd <- projectData(cellchat_vd, PPI.mouse) is not supported currently 18-06-24
cellchat_vd <- computeCommunProb(cellchat_vd)
cellchat_vd <- computeCommunProbPathway(cellchat_vd)
df.net_vd <- subsetCommunication(cellchat_vd)
cellchat_vd <- aggregateNet(cellchat_vd)
cellchat_vd <- netAnalysis_computeCentrality(cellchat_vd)
gc()

#process vehicle
cellchat_vehicle <- subsetData(cellchat_vehicle) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30) # do parallel
cellchat_vehicle <- identifyOverExpressedGenes(cellchat_vehicle)
cellchat_vehicle <- identifyOverExpressedInteractions(cellchat_vehicle)
#cellchat_vehicle <- projectData(cellchat_vehicle, PPI.mouse)
cellchat_vehicle <- computeCommunProb(cellchat_vehicle)
cellchat_vehicle <- computeCommunProbPathway(cellchat_vehicle)
df.net_vehicle <- subsetCommunication(cellchat_vehicle)
cellchat_vehicle <- aggregateNet(cellchat_vehicle)
cellchat_vehicle <- netAnalysis_computeCentrality(cellchat_vehicle)
gc()

#merge it boi
object.list <- list(vehicle = cellchat_vehicle, vd = cellchat_vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

rm(seurat_obj, seurat_vd, seurat_vehicle, meta, df.net_vd, df.net_vehicle, CellChatDB.use, CellChatDB)
#    5.1.3.5 Plots ####

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
png(filename = file.path(plots.dir, paste0(sample.name, "_celltype_interactions.png")),
    width = 16, height = 10, units = "cm", res=300)
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
    width = 6.5, height = 7.75, units = "in", res=300)
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

#    5.1.3.6 Save  ####

saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_sconly.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_sconly.rds")))

#    5.1.3.7 clean ####
rm(list=ls())
gc()

#   5.1.4 lin ####
#    5.1.4.1 values ####
input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_lin_ES.rds"
sample.name = "lin_cellchat"
output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat"
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/sconly/lin"


#    5.1.4.2 load ####
## get the seurat, and subset for vehicle and vd
seurat_obj <- readRDS(input.seurat)
seurat_vehicle <- subset(x = seurat_obj, subset = experiment == "vehicle")
seurat_vd <- subset(x = seurat_obj, subset = experiment == "vd")

## subset for only the top 10% of the senescent cells for each group
seurat_vehicle <- subset(x = seurat_vehicle, subset = SENmayo >= quantile(seurat_vehicle@meta.data[["SENmayo"]], probs = c(0.90)))
seurat_vd <- subset(x = seurat_vd, subset = SENmayo >= quantile(seurat_vd@meta.data[["SENmayo"]], probs = c(0.90)))

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

#    5.1.4.3 setup interaction db ####

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

#    5.1.4.4 Cellchat obj prepros ####

#process intervention
cellchat_vd <- subsetData(cellchat_vd) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30) # do parallel
cellchat_vd <- identifyOverExpressedGenes(cellchat_vd)
cellchat_vd <- identifyOverExpressedInteractions(cellchat_vd)
#cellchat_vd <- projectData(cellchat_vd, PPI.mouse) is not supported currently 18-06-24
cellchat_vd <- computeCommunProb(cellchat_vd)
cellchat_vd <- computeCommunProbPathway(cellchat_vd)
df.net_vd <- subsetCommunication(cellchat_vd)
cellchat_vd <- aggregateNet(cellchat_vd)
cellchat_vd <- netAnalysis_computeCentrality(cellchat_vd)
gc()

#process vehicle
cellchat_vehicle <- subsetData(cellchat_vehicle) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30) # do parallel
cellchat_vehicle <- identifyOverExpressedGenes(cellchat_vehicle)
cellchat_vehicle <- identifyOverExpressedInteractions(cellchat_vehicle)
#cellchat_vehicle <- projectData(cellchat_vehicle, PPI.mouse)
cellchat_vehicle <- computeCommunProb(cellchat_vehicle)
cellchat_vehicle <- computeCommunProbPathway(cellchat_vehicle)
df.net_vehicle <- subsetCommunication(cellchat_vehicle)
cellchat_vehicle <- aggregateNet(cellchat_vehicle)
cellchat_vehicle <- netAnalysis_computeCentrality(cellchat_vehicle)
gc()

#merge it boi
object.list <- list(vehicle = cellchat_vehicle, vd = cellchat_vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

rm(seurat_obj, seurat_vd, seurat_vehicle, meta, df.net_vd, df.net_vehicle, CellChatDB.use, CellChatDB)
#    5.1.4.5 Plots ####

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
png(filename = file.path(plots.dir, paste0(sample.name, "_celltype_interactions.png")),
    width = 16, height = 10, units = "cm", res=300)
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
    width = 6.5, height = 7.75, units = "in", res=300)
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

#    5.1.4.6 Save  ####

saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_sconly.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_sconly.rds")))

#    5.1.4.7 clean ####
rm(list=ls())
gc()

#  5.2 Cellular communication senescence all ####

#   5.2.1 libraries ####
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

#   5.2.2 aem ####
#    5.2.2.1 values ####

input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_aem_ES.rds"
sample.name = "aem_cellchat"
output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/scall/aem"

#    5.2.2.2 input prep ####
## get the object
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

#    5.2.2.3 setup interaction db ####

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

#    5.2.2.4 cellchat obj prepros ####
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
#    5.2.2.5 plots ####
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


#    5.2.2.6 save ####
saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_sc&all.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_sc&all.rds")))

#    5.2.2.7 clean ####
rm(list=ls())
gc()

#   5.2.3 hanai ####
#    5.2.3.1 values ####

input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/aem_hanai_ES.rds"
sample.name = "hanai_cellchat"
output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/scall/hanai"

#    5.2.3.2 input prep ####
## get the object
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

#    5.2.3.3 setup interaction db ####

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

#    5.2.3.4 cellchat obj prepros ####
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
#    5.2.3.5 plots ####
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


#    5.2.3.6 save ####
saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_sc&all.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_sc&all.rds")))

#    5.2.3.7 clean ####
rm(list=ls())
gc()

#   5.2.4 lin ####
#    5.2.4.1 values ####

input.seurat= "/datos/sensence/emilio/VD/output/seurat/enrichment/aem_lin_ES.rds"
sample.name = "lin_cellchat"
output.dir = "/datos/sensence/emilio/VD/output/seurat/cellchat/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/cellchat/scall/lin"

#    5.2.4.2 input prep ####
## get the object
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

#    5.2.4.3 setup interaction db ####

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellchat_vd@DB <- CellChatDB.use
cellchat_vehicle@DB <- CellChatDB.use

#    5.2.4.4 cellchat obj prepros ####
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
#    5.2.4.5 plots ####
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


#    5.2.4.6 save ####
saveRDS(object = cellchat_vd, file = file.path(output.dir, paste0(sample.name, "_vd_sc&all.rds")))
saveRDS(object = cellchat_vehicle, file = file.path(output.dir, paste0(sample.name, "_vehicle_sc&all.rds")))



#    5.2.4.7 clean ####
rm(list=ls())
gc()

