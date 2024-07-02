############ Values #################

input.seurat = "/datos/sensence/emilio/VD/output/seurat/enrichment/enrich_hanai_ES.rds"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"

sample.name = "hanai_GRN_pseudobulk"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat/grn/"
#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/grn/hanai"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"

############ Libraries ##############
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 40)

############ Load and process for WGCNA #################
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

############## Pseudobulk and GRN ##############
datExpr <- ConstructPseudobulk(
  seurat_obj,
  group.by = 'SingleR.labels',
  replicate_col = 'experiment',
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


############### Plots ###############

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


saveRDS(seurat_obj, file = file.path(output.dir, paste0(sample.name,".rds")))

