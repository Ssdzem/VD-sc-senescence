library(Matrix)
library(Seurat)

control_path = "/datos/sensence/emilio/VD/input/sc_data/lin/vehicle/GSM5266942_C5_matrix.tsv"
VD_path = "/datos/sensence/emilio/VD/input/sc_data/lin/vd/GSM5266944_VD_matrix.tsv"

data_control <- read.delim(control_path, header = TRUE, sep = '\t', row.names = "X")
counts_matrix = Matrix(as.matrix(data_control) , sparse=TRUE)
options(Seurat.object.assay.version = "v3")
pbmc.seurat <- CreateSeuratObject(counts = counts_matrix, min.cells = 3, min.features = 200, project = "vehicle")



data_VD <- read.delim(VD_path, header = TRUE, sep = '\t', row.names = "X")
counts_matrix = Matrix(as.matrix(data_VD) , sparse=TRUE)
options(Seurat.object.assay.version = "v3")
pbmc.seurat <- CreateSeuratObject(counts = counts_matrix, min.cells = 3, min.features = 200)

