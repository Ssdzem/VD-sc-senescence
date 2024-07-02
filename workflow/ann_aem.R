###########################################################################
################## Assign Cell type identity ##############################
###########################################################################
#Get the DGE with rows as genes and column as cells 
barcodes_annotation <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_barcodes_anno.csv")
dge <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_dge.csv")
gene_list <- read_csv("/datos/sensence/emilio/VD/input/ref/aem/Adult-Prostate_gene.csv")
#clean data
ref <- merge(gene_list, dge, by = "...1")
ref <- ref[,2:8289]
genes <- ref[,1]
rownames(ref) <- genes
ref <- ref[,2:8288]
ref <- as.matrix(ref)
#his goes in the principal code (seurat.analysis.VD)
#test <- pbmc.seurat@assays[["RNA"]]@counts
#test <- as.matrix(test)
##run SingleR using the expression matrix of the seurat object 
#test_singleR <- SingleR( test = test, ref = ref, labels = barcodes_annotation$Idents.pbmc.)
