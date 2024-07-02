library(tidyverse)
library(vroom)
library(Seurat)

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
######## vibe check  #######

search_value <- "AGCCTCTGCTAAGC-1"

# Use the subset function to filter the rows where the value is located
matching_rows <- subset(cell_ann, barcodes == search_value)

# Print the matching row(s)
if (nrow(matching_rows) > 0) {
  print(matching_rows)
} else {
  cat("The value", search_value, "is not found in the 'barcodes' column.\n")
}

pbmc.test <- subset(pbmc.skin, subset = cell.type == "IFE SB1")
