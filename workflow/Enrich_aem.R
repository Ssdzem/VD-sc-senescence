######### instructions ########
#adjust values of lines 21-25, input is from code "integrate_aem.R"

######## Libraries #########
library(Seurat)
library(GSEABase)
library(escape)
library(AUCell)
library(tidyverse)
library(dittoSeq)
library(gridExtra)
library(patchwork)
library(biomaRt)
library(ggrepel)
library(ggstatsplot)
library(rstatix)
library(BiocParallel)
set.seed(12345)
######## Set values ##########

input.seurat = "/datos/sensence/emilio/VD/output/seurat/prepros/hanai_combined.rds"
sample.name = "enrich_hanai"
output.dir = "/datos/sensence/emilio/VD/output/seurat/enrichment/"
plots.dir = "/datos/sensence/emilio/VD/output/plots/enrichment/hanai/"

######## import datasets #############

### Senescence scores
cellage <- read_csv("input/scores/cellage.csv")

genage <- read_tsv("input/scores/geneage.tsv")

senmayo_mouse <- read_csv("input/scores/senmayo_mouse.csv")

#### Seurat data
aem_combined <- readRDS(input.seurat)

######## Clean data  #######

#### Convert to respective species

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gs_cellage <- getLDS(attributes = c("external_gene_name"),
               filters = "external_gene_name", values = unique(cellage$`Gene symbol`), mart = human,
               attributesL = c("external_gene_name"), martL = mouse)
rm(human,mouse, cellage)

#### Prepros integration

aem_combined <- NormalizeData(aem_combined)
aem_combined <- FindVariableFeatures(aem_combined)
aem_combined <- ScaleData(aem_combined)
aem_combined <- RunPCA(aem_combined)

aem_combined <- FindNeighbors(aem_combined, reduction = "pca", dims = 1:30)
aem_combined <- FindClusters(aem_combined, resolution = 1)
aem_combined <- RunUMAP(aem_combined, dims = 1:30, reduction = "pca")

######## Gene sets ############

gs_cellage <- gs_cellage[,2]
gs_geneage <- unique(genage$`Gene Symbol`)
gs_senmayo <- unique(senmayo_mouse$`Gene(murine)`)

gs_cellage <- GeneSet(gs_cellage)
gs_geneage <- GeneSet(gs_geneage)
gs_senmayo <- GeneSet(gs_senmayo)

geneSets <- list(geneSet1=gs_cellage, geneSet2=gs_geneage, geneSet3=gs_senmayo)

geneSets$geneSet1@setName <- "Cellage"
geneSets$geneSet2@setName <- "GenAge"
geneSets$geneSet3@setName <- "SENmayo"
geneSets <- GeneSetCollection(geneSets)

########## Enrichment ############

ES <- escape.matrix(aem_combined, gene.sets = geneSets, method = "AUCell", groups = 1000, BPPARAM = SnowParam(workers = 40))
aem_combined <- Seurat::AddMetaData(aem_combined, ES)

########### Plot ################

### Vln plot
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot SENMayo by subtype and experiment.png")),
    width = 20, height = 16, units = "in", res=300)
vlnplot_aem <- VlnPlot(
  object = aem_combined,
  features = "SENmayo",
  group.by = "SingleR.labels",
  pt.size = 0.5,
  log = FALSE,
  split.by = "experiment",
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  ncol = NULL,
  slot = "data",
  split.plot = FALSE,
  stack = FALSE,
  combine = TRUE,
  fill.by = "feature",
  flip = FALSE,
  add.noise = TRUE
) + ggtitle(paste0(sample.name,"Violin_plot_SENmayo_by_subtype_and_experiment"))
plot(vlnplot_aem)
dev.off()

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$cellage <- aem_combined@meta.data$Cellage
umap_data$geneage <- aem_combined@meta.data$GenAge
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment
label_counts <- table(umap_data$SingleR.labels)
valid_labels <- names(label_counts[label_counts >= 80])
umap_data <- umap_data %>%
  filter(SingleR.labels %in% valid_labels)

#see the overall violin plot 
png(filename = file.path(plots.dir, paste0(sample.name, "_violin plot overrall SENmayo.png")),
    width = 20, height = 16, units = "in", res=300)
overal_plot <- ggbetweenstats(
  data  = umap_data,
  x     = experiment,
  y     = SENmayo,
  title = paste0(sample.name, "__overall_SENmayo score values"),
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 7, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
  #pairwise.display = "none"
)
overal_plot + 
  theme(axis.text.x=element_text(size = 17.5))  + theme(plot.title=element_text(size=25))
dev.off()
### Significance 

significance <- data.frame(
  Variable = character(),
  Wilcoxon_Result = character(),
  P_Value = numeric(),
  Median_Vehicle = numeric(),
  Median_VD = numeric(),
  stringsAsFactors = FALSE
)

# List of variables to compare
variables_to_compare <- c("Cellage", "GenAge", "SENmayo")

# Iterate through each variable
for (variable in variables_to_compare) {
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(get(variable) ~ experiment, data = aem_combined@meta.data)
  
  # Check if the Wilcoxon statistic is numeric
  wilcox_statistic <- ifelse(is.numeric(wilcox_result$statistic), as.character(wilcox_result$statistic), NA)
  
  # Calculate median for each group
  median_vehicle <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vehicle", variable], na.rm = TRUE)
  median_vd <- median(aem_combined@meta.data[aem_combined@meta.data$experiment == "vd", variable], na.rm = TRUE)
  
  # Create a row for the result dataframe
  result_row <- data.frame(
    Variable = variable,
    Wilcoxon_Result = wilcox_statistic,
    P_Value = wilcox_result$p.value,
    Median_Vehicle = median_vehicle,
    Median_VD = median_vd,
    stringsAsFactors = FALSE
  )
  
  # Add the row to the result dataframe
  significance <- bind_rows(significance, result_row)
}

# Print the resulting dataframe
print(significance)
write.csv(x = significance, 
          file = file.path(output.dir, paste0(sample.name, "_significance.csv")))

### UMAP

#grab the data for plotting
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment

# Choose the specific labels you want to label
labels_to_label <- unique(umap_data$SingleR.labels)

# Create UMAP plot using ggplot2 with labels for each unique cell subtype and with certain experiment
base_vehicle <- ggplot(umap_data %>% filter(experiment == "vehicle"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )   +
  ggtitle("AEM_Vehicle")

# Display the plot
print(base_vehicle)

#again but with different experiment
base_vd <- ggplot(umap_data %>% filter(experiment == "vd"), aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels, color = NULL,group= NULL),
    size = 3, color = "black", box.padding = 0.5, 
  )  +
  ggtitle("AEM_VD")

# Display the plot
print(base_vd)

#print them together
png(filename = file.path(plots.dir, paste0(sample.name, "_UMAP SENmayo experiments.png")),
    width = 20, height = 16, units = "in", res=300)
print((base_vehicle | base_vd) +
        plot_annotation(title = paste0("\n",sample.name,"\n"),
                        theme = theme(plot.title = element_text(face = "bold", size = 22, hjust=0.5)))
)
dev.off()

#### Pairwise wilcoxon top 3 cell subtypes

pairwise.wilcox.test(umap_data$SENmayo, umap_data$SingleR.labels,
                     p.adjust.method = "BH")
# Subset the data for "vehicle" and "vd" experiments
vehicle_data <- aem_combined@meta.data %>%
  filter(experiment == "vehicle")

vd_data <- aem_combined@meta.data %>%
  filter(experiment == "vd")

# Get unique labels
unique_labels <- unique(aem_combined@meta.data$SingleR.labels)
# Initialize a data frame to store the results
wilcox_results <- data.frame(
  SingleR_label = character(),
  wilcox_value = numeric(),
  p_value = numeric(),
  experiment_comparison = character(),
  VD_minus_vehicle = numeric(), # New column for the mean difference
  stringsAsFactors = FALSE
)

# Perform Wilcoxon tests for each label
for (label in unique_labels) {
  # Subset data for the current label
  vehicle_label_data <- vehicle_data %>%
    filter(SingleR.labels == label)
  
  vd_label_data <- vd_data %>%
    filter(SingleR.labels == label)
  
  # Check if both datasets have enough observations
  if (nrow(vehicle_label_data) >= 3 && nrow(vd_label_data) >= 3) {
    # Calculate means for each group
    vehicle_mean <- mean(vehicle_label_data$SENmayo)
    vd_mean <- mean(vd_label_data$SENmayo)
    
    # Calculate the difference in means
    mean_difference <- vd_mean - vehicle_mean
    
    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(vehicle_label_data$SENmayo, vd_label_data$SENmayo)
    
    # Store the results along with the calculated mean difference
    wilcox_results <- rbind(
      wilcox_results,
      data.frame(
        SingleR_label = label,
        wilcox_value = wilcox_result$statistic,
        p_value = wilcox_result$p.value,
        experiment_comparison = "vehicle vs vd",
        VD_minus_vehicle = mean_difference,  # Store the mean difference
        stringsAsFactors = FALSE
      )
    )
  } else {
    # Print a message or handle cases where there are not enough observations
    cat("Skipping label", label, "due to insufficient observations.\n")
  }
}

# Order the results by wilcox_value
wilcox_results <- wilcox_results %>% arrange(p_value)

# Get the top 3 results
top_results <- head(wilcox_results, 3)

# Print the top results
print(top_results)
write.csv(x = top_results, 
          file = file.path(output.dir, paste0(sample.name, "_top_3_wilcox_cellsubtypes.csv")))

############# Save ###########

saveRDS(object = aem_combined, file = paste0(output.dir,sample.name, "_ES.rds"))
