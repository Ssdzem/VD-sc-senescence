####### instructions ##############

# adjust values (line 15-20) and enable lines 27-31 for hanai dataset

####### set up libraries and variables #########

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)

######## Set values ##########

input.seurat.vd = "/datos/sensence/emilio/VD/output/seurat/prepros/hanai_vd.rds"
input.seurat.vehicle = "/datos/sensence/emilio/VD/output/seurat/prepros/hanai_vehicle.rds"
sample.name = "hanai_combined"
output.dir = "/datos/sensence/emilio/VD/output/seurat/prepros/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/prepros/hanai/"
set.seed(12345)

######## import datasets ########

vehicle <- readRDS(input.seurat.vehicle)
vd <- readRDS(input.seurat.vd)

# Here we will start the purging of the immune cells since they are from bone marrow and adds much noise 
vehicle <- subset(vehicle, SingleR.labels != 'Immune cells')
saveRDS(vehicle, "/datos/sensence/emilio/VD/output/seurat/hanai_control.rds")
vd <- subset (vd, SingleR.labels != 'Immune cells')
saveRDS(vd, "/datos/sensence/emilio/VD/output/seurat/hanai_vd.rds")

####### Set idents and metadata #######

Idents(vehicle) <- 'vehicle'
df <- data.frame(Idents(vehicle))
vehicle <- AddMetaData(object = vehicle, metadata = df, col.name = 'experiment')
Idents(vehicle)

Idents(vd) <- 'vd'
df <- data.frame(Idents(vd))
vd <- AddMetaData(object = vd, metadata = df, col.name = 'experiment')
Idents(vd)


########### Merge/integration/prepros ###########

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

######### Plotting individually #########
combined <- RunUMAP(combined, dims = 1:30, reduction = "harmony")

DimPlot(combined, reduction = "umap", group.by = c("experiment", "SingleR.labels"))
umap <- DimPlot(combined, reduction = "umap", split.by = "experiment", group.by = ("SingleR.labels"))
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP_annotated_integrated.png")),
    width = 20, height = 16, units = "in", res=300)
umap + labs(title= paste0(sample.name, "_UMAP_annnotated"))
dev.off()

######## Save it boi #####

saveRDS(object = combined, file = file.path(output.dir, paste0(sample.name, ".rds")))
