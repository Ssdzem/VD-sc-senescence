######## Experimento con SIngle R ########

#pbmc.seurat <- Read10X("/datos/sensence/emilio/VD/input/Abu_el_maaty_2021/ref/Adult_Prostate")
Adult_Prostate_barcodes_anno <- read_csv("input/Abu_el_maaty_2021/ref/Adult_Prostate/Adult-Prostate_barcodes_anno.csv")
Adult_Prostate_dge <- read_csv("input/Abu_el_maaty_2021/ref/Adult_Prostate/Adult-Prostate_dge.csv")
Adult_Prostate_gene <- read_csv("input/Abu_el_maaty_2021/ref/Adult_Prostate/Adult-Prostate_gene.csv")
ref <- merge(Adult_Prostate_gene, Adult_Prostate_dge, by = "...1")
ref <- ref[,2:8289]
genes <- ref[,1]
rownames(ref) <- genes
ref <- ref[,2:8288]
ref <- as.matrix(ref)
test <- as.matrix(test)
test_singleR <- SingleR( test = test, ref = ref, labels = Adult_Prostate_barcodes_anno$Idents.pbmc.)
pbmc.seurat[["SingleR.labels"]] <- test_singleR$labels
library(SingleR)
