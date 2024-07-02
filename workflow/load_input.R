################## instructions ################

#This code needs the metadata input files to get the rest of the input data, check it out
#for directory management, it will set on the working direcotry on R, so set your directory where you will manage all files of this project

################# values ######################

scores.path = paste0(getwd(),"/metadata/scores.csv")
ref.path = paste0(getwd(),"/metadata/ref.csv")
sc_rna.path = paste0(getwd(),"/metadata/input.csv")
scores.dir = paste0(getwd(),"/input/scores/")
ref.dir = paste0(getwd(),"/input/ref/")
sc_data.dir = paste0(getwd(),"/input/sc_data/")
input.dir =paste0(getwd(),"/input/")



################ libraries ####################

library(downloader)
library(readr)
library(readxl)
library(googledrive)

############### load data ####################
#0. make the input directory
dir.create(input.dir, recursive = T)

#1. load the metadata files to fill the input directory
scores <- read_csv(scores.path)
ref <- read_csv(ref.path)
input <- read_csv(sc_rna.path)

############## Senescent genes sets #########

#0. make a directory to save it for further analysis
dir.create(scores.dir, recursive = T)

# 1. Download and extract the TSV file from the ZIP (geneage)
download(scores$`download link`[1], destfile = paste0(scores.dir,"geneage.tsv") , mode = "wb")

# 2. Download cellage (from zip)
zip_file <- tempfile(fileext = ".zip")
download(scores$`download link`[2], destfile = zip_file, mode = "wb")
temp_dir <- tempdir()
unzip(zip_file, exdir = temp_dir)
cellage_tsv_path <- file.path(temp_dir, "cellage3.tsv")
cellage <- read_tsv(cellage_tsv_path)
write_csv(cellage, paste0(scores.dir, "cellage.csv"))


# 3. Download and read the XLSX file (senmayo)
download(scores$`download link`[3], destfile = paste0(scores.dir,"senmayo.xlsx"), mode = "wb")
senmayo_mouse <- read_excel(paste0(scores.dir,"senmayo.xlsx"), sheet = "mouse")
write_csv(senmayo_mouse, paste0(scores.dir, "senmayo_mouse.csv"))

#############  Ref (annotation) of the input data ############

#0. make the dirs to operate normally

dir.create(ref.dir, recursive = T)
#for aem annotation ref data
dir.create(paste0(ref.dir,"aem/"))
#for hanai annotation ref data
dir.create(paste0(ref.dir,"hanai/"))
#for lin annotation ref data
dir.create(paste0(ref.dir,"lin/"))

#1. aem input
aemref_url <- ref$`Download link`[1]
zip_file <- tempfile(fileext = ".zip")
download(aemref_url, destfile = zip_file, mode = "wb")
unzip(zip_file, exdir = paste0(ref.dir,"aem/"))

#2. hanai input 
drive_deauth()
drive_user()
public_file <-  drive_get(as_id(ref$`Download link`[5]))
drive_download(public_file, overwrite = TRUE, path = paste0(ref.dir,"hanai/","SCA_v.0.5.0_preprint.loom"))

#3. lin input 
linref_barcodes_url <- ref$`Download link`[2]
download(linref_barcodes_url, destfile = paste0(ref.dir,"lin/barcodes.tsv.gz"), mode = "wb")
linref_features_url <- ref$`Download link`[3]
download(linref_features_url, destfile = paste0(ref.dir,"lin/features.tsv.gz"), mode = "wb")
linref_matrix_url <- ref$`Download link`[4]
download(linref_matrix_url, destfile = paste0(ref.dir,"lin/matrix.mtx.gz"), mode = "wb")

############ Input data (scRNAseq data) ##############

#0. make the dirs to operate normally
dir.create(sc_data.dir, recursive = T)
#for aem annotation ref data
dir.create(paste0(sc_data.dir,"aem/"))
dir.create(paste0(sc_data.dir,"aem/vd/"))
dir.create(paste0(sc_data.dir,"aem/vehicle/"))
#for hanai annotation ref data
dir.create(paste0(sc_data.dir,"hanai/"))
dir.create(paste0(sc_data.dir,"hanai/vd/"))
dir.create(paste0(sc_data.dir,"hanai/vehicle/"))
#for lin annotation ref data
dir.create(paste0(sc_data.dir,"lin/"))
dir.create(paste0(sc_data.dir,"lin/vd/"))
dir.create(paste0(sc_data.dir,"lin/vehicle/"))

#1. aem vd input
download(input$matrix[1], destfile = paste0(sc_data.dir,"aem/vd/matrix.mtx.gz"), mode = "wb")
download(input$features[1], destfile = paste0(sc_data.dir,"aem/vd/features.tsv.gz"), mode = "wb")
download(input$barcodes[1], destfile = paste0(sc_data.dir,"aem/vd/barcodes.tsv.gz"), mode = "wb")

#2. aem vehicle input 
download(input$matrix[2], destfile = paste0(sc_data.dir,"aem/vehicle/matrix.mtx.gz"), mode = "wb")
download(input$features[2], destfile = paste0(sc_data.dir,"aem/vehicle/features.tsv.gz"), mode = "wb")
download(input$barcodes[2], destfile = paste0(sc_data.dir,"aem/vehicle/barcodes.tsv.gz"), mode = "wb")

#3. hanai vd input 
download(input$matrix[3], destfile = paste0(sc_data.dir,"hanai/vd/matrix.mtx.gz"), mode = "wb")
download(input$features[3], destfile = paste0(sc_data.dir,"hanai/vd/features.tsv.gz"), mode = "wb")
download(input$barcodes[3], destfile = paste0(sc_data.dir,"hanai/vd/barcodes.tsv.gz"), mode = "wb")

#4. hanai vehicle input 
download(input$matrix[4], destfile = paste0(sc_data.dir,"hanai/vehicle/matrix.mtx.gz"), mode = "wb")
download(input$features[4], destfile = paste0(sc_data.dir,"hanai/vehicle/features.tsv.gz"), mode = "wb")
download(input$barcodes[4], destfile = paste0(sc_data.dir,"hanai/vehicle/barcodes.tsv.gz"), mode = "wb")

#5. lin vd input 
download(input$matrix[5], destfile = paste0(sc_data.dir,"lin/vd/GSM5266944_VD_matrix.tsv"), mode = "wb")

#6. lin vehicle input 
download(input$matrix[6], destfile = paste0(sc_data.dir,"lin/vehicle/GSM5266942_C5_matrix.tsv"), mode = "wb")


######## output directories #####
# Define the base directory
base.dir <- getwd()

# Create the output directory
output.dir <- file.path(base.dir, "output")
dir.create(output.dir, recursive = TRUE)

# Define the subdirectories structure relative to the output directory
plots.subdirs <- c(
  "plots/cellchat/scall",
  "plots/cellchat/sconly",
  "plots/cellchat/scall/hanai",
  "plots/cellchat/scall/aem",
  "plots/cellchat/scall/lin",
  "plots/cellchat/sconly/hanai",
  "plots/cellchat/sconly/aem",
  "plots/cellchat/sconly/lin",
  "plots/enrichment/aem",
  "plots/enrichment/hanai",
  "plots/enrichment/lin",
  "plots/grn/aem",
  "plots/grn/hanai",
  "plots/grn/lin",
  "plots/prepros/aem/vd",
  "plots/prepros/aem/vehicle",
  "plots/prepros/hanai/vd",
  "plots/prepros/hanai/vehicle",
  "plots/prepros/lin/vd",
  "plots/prepros/lin/vehicle",
  "seurat/cellchat",
  "seurat/enrichment",
  "seurat/grn",
  "seurat/prepros"
)

# Create all directories inside the output directory
for (subdir in plots.subdirs) {
  full.path <- file.path(paste0(output.dir, "/",subdir))
  dir.create(full.path, recursive = TRUE, showWarnings = FALSE)
}

# Check the directory structure
system(paste("tree", base.dir))
