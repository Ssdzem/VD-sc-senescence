###### Variables ####

input.seurat = "/datos/sensence/emilio/VD/output/seurat/hanai_combined.rds"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"

sample.name = "hanai_combined_GRN"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat"
#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/pub/chondro_shenanigans"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"

# cell.type = "Chondrocytes"
#check Seurat object meta data i dont know why this doesnt work

###### Set up ######
suppressPackageStartupMessages({ 
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(igraph)
})

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 50)
seurat_myobj <- readRDS(input.seurat)
seurat_myobj[["RNA3"]] <- as(object = seurat_myobj[["RNA"]], Class = "Assay")
DefaultAssay(object = seurat_myobj) <- "RNA3"
seurat_myobj[["RNA"]] <- NULL
RenameAssays(object = seurat_myobj, assay.name = "RNA3", new.assay.name = "RNA")

###### Vibe check ######
p <- DimPlot(seurat_myobj, group.by='SingleR.labels', label=F) +
  umap_theme() + ggtitle(sample.name)
p

###### Setup for WGCNA ######
seurat_myobj <- SetupForWGCNA(
  seurat_myobj,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "sensence" 
)

###### Create metacells (ominous) ######
seurat_myobj <- MetacellsByGroups(
  seurat_obj = seurat_myobj,
  group.by = c("SingleR.labels"), 
  reduction = 'umap', 
  k = 25, 
  max_shared = 10, 
  ident.group = 'SingleR.labels' 
)
seurat_myobj <- NormalizeMetacells(seurat_myobj)

####### process metacells you never know ########
metacell_myobj <- GetMetacellObject(seurat_myobj)

seurat_myobj <- NormalizeMetacells(seurat_myobj)
seurat_myobj <- FindVariableFeatures(object = seurat_myobj)
seurat_myobj <- ScaleMetacells(seurat_myobj, features=VariableFeatures(seurat_myobj))
seurat_myobj <- RunPCAMetacells(seurat_myobj, features=VariableFeatures(seurat_myobj))
seurat_myobj <- RunUMAPMetacells(seurat_myobj, dims=1:15)

p1 <- DimPlotMetacells(seurat_myobj, group.by='SingleR.labels') + umap_theme() + ggtitle("Cell Type")
p1 

####### Coexpression network analysis #########

seurat_myobj <- SetDatExpr(
  seurat_myobj,
  group_name = "Chondrocytes", # the name of the group of interest in the group.by column
  group.by='SingleR.labels', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA3', # using RNA assay
  slot = 'data' # using normalized data
)

####### Selecting soft power threshold ######

# Test different soft powers:
seurat_myobj <- TestSoftPowers(
  seurat_myobj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_myobj)

# assemble with patchwork
png(filename = file.path(plots.dir, paste0(sample.name, "Chondrocytes", "_soft_powers.png")),
    width = 20, height = 16, units = "in", res=300)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_myobj)
head(power_table)
p1 <- wrap_plots(plot_list, ncol=2)
########################## Construct co expression network ########################################

# construct co-expression network:
seurat_myobj <- ConstructNetwork(
  seurat_myobj, soft_power=9,
  setDatExpr=FALSE,
  tom_name = "Chondrocytes", # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE
)
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_Dendrogram.png")),
    width = 20, height = 16, units = "in", res=300)
PlotDendrogram(seurat_myobj, main='Chondrocytes hdWGCNA Dendrogram')
dev.off()

TOM <- GetTOM(seurat_myobj)
####### Module Eigengenes and Connectivity #######

# need to run ScaleData first or else harmony throws an error:
seurat_myobj <- ScaleData(seurat_myobj, features=VariableFeatures(seurat_myobj))

# compute all MEs in the full single-cell dataset
seurat_myobj <- ModuleEigengenes(
  seurat_myobj,
  verbose = T,
)

# harmonized module eigengenes:
MEs <- GetMEs(seurat_myobj)

# module eigengenes:
MEs <- GetMEs(seurat_myobj, harmonized=FALSE)
####### Compute module connectivity #######

# compute eigengene-based connectivity (kME):
seurat_myobj <- ModuleConnectivity(
  seurat_myobj,
  group.by = 'SingleR.labels', group_name = "Chondrocytes"
)
# rename the modules
seurat_myobj <- ResetModuleNames(
  seurat_myobj,
  new_name = paste0(sample.name, "_", "Chondrocytes") 
)
# plot genes ranked by kME for each module
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_module_kME.png")),
    width = 20, height = 16, units = "in", res=300)
PlotKMEs(seurat_myobj, ncol=5)
dev.off()


####### Getting the module assignment table ######
# get the module assignment table:
modules <- GetModules(seurat_myobj)

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
myhub_df <- GetHubGenes(seurat_myobj, n_hubs = 10)

head(myhub_df)

saveRDS(seurat_myobj, file = file.path(output.dir, paste(sample.name, "Chondrocytes", '_hdWGCNA_object.rds', sep = "_")))

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_myobj <- ModuleExprScore(
  seurat_myobj,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_myobj <- ModuleExprScore(
  seurat_myobj,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_myobj, file = file.path(output.dir, paste(sample.name, "Chondrocytes", '_hdWGCNA_object.rds', sep = "_")))

####### Module feature plots ###########

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_myobj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_featureplot_kME.png")),
    width = 20, height = 16, units = "in", res=300)
wrap_plots(plot_list, ncol=6)
dev.off()

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_myobj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_featureplot_hubgenes.png")),
    width = 20, height = 16, units = "in", res=300)
wrap_plots(plot_list, ncol=6)
dev.off()

####### Module Correlations ##########

# plot module correlagram
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_correlogram.png")),
    width = 20, height = 16, units = "in", res=300)
ModuleCorrelogram(seurat_myobj)
dev.off()

####### Seurat plotting functions ########
# get hMEs from seurat object
MEs <- GetMEs(seurat_myobj, harmonized=F)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_myobj@meta.data <- cbind(seurat_myobj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_myobj, features=mods, group.by = 'SingleR.labels')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_dotplot_hME.png")),
    width = 20, height = 16, units = "in", res=300)
p
dev.off()

####### Individual module plots ###########
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_individual_network_hME.png")),
    width = 20, height = 16, units = "in", res=300)
p2 <- ModuleNetworkPlot(seurat_myobj,  return_graph=TRUE, mods= "hanai_combined_GRN_Chondrocytes13", n_inner = 5, n_outer = 10)
print(p2)
dev.off()
####### Combined hub gene network plots #########

# hubgene network
png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_hubgene_network.png")),
    width = 20, height = 16, units = "in", res=300)
HubGeneNetworkPlot(
  seurat_myobj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
g <- HubGeneNetworkPlot(seurat_myobj,  return_graph=TRUE)
g
dev.off()
####### Applying UMAP to co-expression networks #######

png(filename = file.path(plots.dir, paste0(sample.name,"Chondrocytes", "_GRN_UMAP_.png")),
    width = 20, height = 16, units = "in", res=300)
seurat_myobj <- RunModuleUMAP(
  seurat_myobj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
# get the hub gene UMAP table from the seurat object
umap_mydf <- GetModuleUMAP(seurat_myobj)

# plot with ggplot
ggplot(umap_mydf, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_mydf$color, # color each point by WGCNA module
    size=umap_mydf$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_myobj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)

g <- ModuleUMAPPlot(seurat_myobj,  return_graph=TRUE)
dev.off()

saveRDS(seurat_myobj, file = file.path(plots.dir, paste(sample.name, "Chondrocytes", '_hdWGCNA_object.rds', sep = "_")))
