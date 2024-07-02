library(Seurat)
library(SeuratData)
library(SeuratDisk)

pbmc.bone <- Connect(filename = "/datos/sensence/emilio/VD/input/ref/hanai/SCA_v.0.5.0_preprint.loom", mode = "r")
pbmc.bone <- as.Seurat(pbmc.bone)
selected_cell_types <- c("Osteoclasts", "Osteoblasts","Chondrocytes", "Endothelium", "LepR+ BM-MSCs", "Fibroblasts", "Osteocytes", "Pericytes","Periosteum", "Adipocytes", "Immune cells", "Vascular smooth muscle", "Tenocytes")
pbmc.boneM <- subset(pbmc.bone, subset = CellType %in% selected_cell_types)
