############################## seurat mode ########################################
library(Seurat)
library(SeuratData)
library(ggplot2)
ifnb <- LoadData("ifnb")

# Normalize the data
ifnb <- NormalizeData(ifnb)

# Find DE features between CD16 Mono and CD1 Mono
Idents(ifnb) <- "seurat_annotations"
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = "CD4 Naive T_CTRL")
# view results
head(monocyte.de.markers)
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD4 Memory T_CTRL", verbose = FALSE)
head(mono.de, n = 10)

# load the inferred sample IDs of each cell
ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)

# only keep the cells with high-confidence sample ID
info <- info[grep(pattern = "SNG", x = info$BEST), ]

# remove cells with duplicated IDs in both ctrl and stim groups
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]

# now add the sample IDs to ifnb 
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
ifnb <- AddMetaData(ifnb, metadata = info)

# remove cells without donor IDs
ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
ifnb <- subset(ifnb, subset = donor_id != "unknown")

#pseudobulk
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))
pseudo_ifnb$celltype.stim <- paste(ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
Idents(pseudo_ifnb) <- "celltype.stim"

bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "CTRL_SNG-101_CD14 Mono", 
                            ident.2 = "STIM_SNG-1016_CD14 Mono",
                            test.use = "DESeq2")
head(bulk.mono.de, n = 15)



##################################### hdWGCNA mode ###############################

# my idea is that i need a wgcna and metacells preprocess to start the pseudobulk, lets try it

# single-cell analysis package
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
enableWGCNAThreads(nThreads = 50)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('Zhou_2020.rds')
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "INH", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
datExpr <- ConstructPseudobulk(
  seurat_obj,
  group.by = 'cell_type',
  replicate_col = 'Sample',
  label_col = 'Sample',
  assay = 'RNA',
  slot = 'counts', # this should always be counts!
  min_reps = 10
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
# we only want the data from the astrocytes 
cur_group <- 'ASC'

# subset the matrix for just this cell type
cur_cpm <- cpm[grepl(cur_group, rownames(cpm)),]

seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = cur_cpm
)

# select the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

# construct the co-expression network and identify gene modules
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  tom_name='pseudobulk', 
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

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient(high='red', low='grey95') + 
  xlab('') + ylab('')


p 