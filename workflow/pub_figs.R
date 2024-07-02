# 1.PREPROS ##########
#  1.1 values ####
plots.dir = paste0(getwd(),"/output/plots/pub/fig1")
input.dir = paste0(getwd(), "/output/seurat/prepros/")
#  1.2 libraries ####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)

#  1.3 load ####
aem_combined <- readRDS(paste0(input.dir, "aem_combined.rds"))
hanai_combined <- readRDS(paste0(input.dir, "hanai_combined.rds"))
lin_combined <- readRDS(paste0(input.dir, "lin_combined.rds"))

#  1.4 Figure 1 ####
#   1.4.1 Figure 1A ####
aem_umap <- DimPlot(aem_combined, reduction = "umap", split.by = "experiment", group.by = "SingleR.labels")  +
  ggtitle("Prostate VD vs Vehicle anno-UMAP") + 
  theme(plot.title = element_text(hjust = 0.5))
A <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "A", hjust = -0.1, vjust = -0.1, size = 10, fontface = "bold") +
  theme_void()
aem_umap <- aem_umap + inset_element(A, left = 0, bottom = 1, right = 0.1, top = 1.1)

#   1.4.2 FIgure 2C ####
hanai_umap <- DimPlot(hanai_combined, reduction = "umap", split.by = "experiment", group.by = ("SingleR.labels")) +
  ggtitle("Bone VD vs Vehicle anno-UMAP") + 
  theme(plot.title = element_text(hjust = 0.5))
C <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "C", hjust = -0.1, vjust = -0.1, size = 10, fontface = "bold") +
  theme_void()
hanai_umap <- hanai_umap + inset_element(C, left = 0, bottom = 1, right = 0.1, top = 1.1)

#   1.4.3 Figure 2B ####
lin_umap <- DimPlot(lin_combined, reduction = "umap", split.by = "experiment", group.by = ("SingleR.labels")) +
  ggtitle("Skin VD vs Vehicle anno-UMAP") + 
  theme(plot.title = element_text(hjust = 0.5))
B <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "B", hjust = -0.1, vjust = -0.1, size = 10, fontface = "bold") +
  theme_void()
lin_umap <- lin_umap + inset_element(B, left = 0, bottom = 1, right = 0.1, top = 1.1)

#   1.4.4 Overrall ####
umap_pub <- aem_umap + hanai_umap | lin_umap + plot_spacer()
umap_pub <- umap_pub + 
  plot_annotation(
    title = "Figure 1. UMAP plots of VD vs Vehicle for each tissue",
    theme = theme(
      plot.title = element_text(hjust = 0, size = 25, face = "bold")  # Adjusted for size and alignment
    )
  )

png(filename = file.path(plots.dir, paste0("figure_1.png")),
    width = 18, height = 16, units = "in", res=300)
umap_pub
dev.off()




#   1.4.5 clean ####
rm(list=ls())
gc()

# 2.ENRICHMENT ####
#  2.1 Values #############

plots.dir = paste0(getwd(),"/output/plots/pub/fig2")
input.dir = paste0(getwd(), "/output/seurat/enrichment/")

#  2.2 Libraries #################
library(Seurat)
library(ggstatsplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(tidyverse)

set.seed(12345)

# ---Enrichment
#  2.3 load ####

aem <- readRDS(paste0(input.dir, "enrich_aem_ES.rds"))
hanai <- readRDS(paste0(input.dir, "enrich_hanai_ES.rds"))
lin <- readRDS(paste0(input.dir, "enrich_lin_ES.rds"))

#  2.4 FIgure 2 ####
#   2.4.1 Figure 2A ####
#prostate dataframe

umap_aem <- as.data.frame(Embeddings(aem, "umap"))
umap_aem$SENmayo <- aem@meta.data$SENmayo
umap_aem$SingleR.labels <- aem@meta.data$SingleR.labels
umap_aem$experiment <- aem@meta.data$experiment
label_counts <- table(umap_aem$SingleR.labels)
valid_labels <- names(label_counts[label_counts >= 80])
umap_aem <- umap_aem %>%
  filter(SingleR.labels %in% valid_labels)

#Prostate violin-boxplot

aem_ES <- ggbetweenstats(
  data = umap_aem,
  x = experiment,
  y = SENmayo,
  title = "Prostate overall SENmayo score violin-boxplots",
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  messages = FALSE,
  plot = FALSE
)

category_counts <- umap_aem %>% 
  group_by(experiment) %>%
  summarise(count = n(), .groups = 'drop')

experiment_labels <- category_counts %>% 
  mutate(label = paste(experiment, "\n(n =", count, ")")) %>%
  pull(label)

aem_ES$layers <- aem_ES$layers[-1]

final_plot <- ggplot(umap_aem, aes(x = experiment, y = SENmayo)) +
  geom_point(aes(color = SingleR.labels), position = position_jitter(width = 0.1, height = 0), alpha = 0.6) + # Add custom points first
  scale_color_viridis_d(name = "Cell Type")

for (layer in aem_ES$layers) {
  final_plot <- final_plot + layer
}

aem_ES <- final_plot +
  labs(
    title = "Prostate overall SENmayo score violin-boxplots",
    subtitle = aem_ES$labels$subtitle,
    caption = aem_ES$labels$caption,
    x = NULL,
    y = aem_ES$labels$y
  ) +
  scale_x_discrete(labels = experiment_labels) +
  theme_minimal() +
  theme(legend.position = "right")

print(aem_ES)

#   2.4.2 Figure 2B ####
#Bone dataframe

umap_hanai <- as.data.frame(Embeddings(hanai, "umap"))
umap_hanai$SENmayo <- hanai@meta.data$SENmayo
umap_hanai$SingleR.labels <- hanai@meta.data$SingleR.labels
umap_hanai$experiment <- hanai@meta.data$experiment
umap_hanai <- umap_hanai %>%
  dplyr::filter(!SingleR.labels ==  "Immune cells")

#Bone violin-boxplot

hanai_ES <- ggbetweenstats(
  data = umap_hanai,
  x = experiment,
  y = SENmayo,
  title = "Bone overall SENmayo score violin-boxplots",
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  messages = FALSE,
  plot = FALSE
)

category_counts <- umap_hanai %>% 
  group_by(experiment) %>%
  summarise(count = n(), .groups = 'drop')

experiment_labels <- category_counts %>% 
  mutate(label = paste(experiment, "\n(n =", count, ")")) %>%
  pull(label)

hanai_ES$layers <- hanai_ES$layers[-1]

final_plot <- ggplot(umap_hanai, aes(x = experiment, y = SENmayo)) +
  geom_point(aes(color = SingleR.labels), position = position_jitter(width = 0.1, height = 0), alpha = 0.6) + 
  scale_color_viridis_d(name = "Cell Type")

for (layer in hanai_ES$layers) {
  final_plot <- final_plot + layer
}

hanai_ES <- final_plot +
  labs(
    title = "Bone overall SENmayo score violin-boxplots",
    subtitle = hanai_ES$labels$subtitle,
    caption = hanai_ES$labels$caption,
    x = NULL, 
    y = hanai_ES$labels$y
  ) +
  scale_x_discrete(labels = experiment_labels) + 
  theme_minimal() +
  theme(legend.position = "right")

print(hanai_ES)

#   2.4.3 Figure 2C ####
#Skin dataframe

umap_lin <- as.data.frame(Embeddings(lin, "umap"))
umap_lin$SENmayo <- lin@meta.data$SENmayo
umap_lin$SingleR.labels <- lin@meta.data$SingleR.labels
umap_lin$experiment <- lin@meta.data$experiment
#label_counts <- table(umap_lin$SingleR.labels)
#valid_labels <- names(label_counts[label_counts >= 80])
#umap_lin <- umap_lin %>%
#  filter(SingleR.labels %in% valid_labels)

#Skin violinboxplot
lin_ES <- ggbetweenstats(
  data = umap_lin,
  x = experiment,
  y = SENmayo,
  title = "Prostate overall SENmayo score violin-boxplots",
  type = "nonparametric",
  k = 3,
  ggsignif.args = list(textsize = 10),
  messages = FALSE,
  plot = FALSE
)

category_counts <- umap_lin %>% 
  group_by(experiment) %>%
  summarise(count = n(), .groups = 'drop')

experiment_labels <- category_counts %>% 
  mutate(label = paste(experiment, "\n(n =", count, ")")) %>%
  pull(label)

lin_ES$layers <- lin_ES$layers[-1]

final_plot <- ggplot(umap_lin, aes(x = experiment, y = SENmayo)) +
  geom_point(aes(color = SingleR.labels), position = position_jitter(width = 0.1, height = 0), alpha = 0.6) +
  scale_color_viridis_d(name = "Cell Type")

for (layer in lin_ES$layers) {
  final_plot <- final_plot + layer
}

lin_ES <- final_plot +
  labs(
    title = "Skin overall SENmayo score violin-boxplots",
    subtitle = lin_ES$labels$subtitle,
    caption = lin_ES$labels$caption,
    x = NULL, 
    y = lin_ES$labels$y
  ) +
  scale_x_discrete(labels = experiment_labels) + 
  theme_minimal() +
  theme(legend.position = "right")

print(lin_ES)

#   2.4.4 overall ####

#Reference each plot with letters

A <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "A", hjust = -0.1, vjust = -0.1, size = 10, fontface = "bold") +
  theme_void()
aem_ES <- aem_ES + inset_element(A, left = 2, bottom = 1, right = 0, top = 1.1)

B <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "B", hjust = -0.1, vjust = -0.1, size = 10, fontface = "bold") +
  theme_void()
hanai_ES <- hanai_ES + inset_element(B, left = 2, bottom = 1, right = 0, top = 1.1)

C <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "C", hjust = -0.1, vjust = -0.1, size = 10, fontface = "bold") +
  theme_void()
lin_ES <- lin_ES + inset_element(C, left = 2, bottom = 1, right = 0, top = 1.1)

umap_pub <- aem_ES | hanai_ES | lin_ES

umap_pub <- umap_pub + 
  plot_annotation(
    title = "Figure 2. Enrichment violin-boxplots of VD vs Vehicle for each tissue",
    theme = theme(
      plot.title = element_text(hjust = 0, size = 25, face = "bold")  # Adjusted for size and alignment
    )
  )


png(filename = file.path(plots.dir, paste0("figure_2.png")),
    width = 24, height = 8, units = "in", res=300)
umap_pub
dev.off()

#  2.5 Supplementary table 1 wilcox senescent gene sets #####


# Define the IDs and the result dataframe
table_id <- data.frame(id = c("aem", "lin", "hanai"))
significance <- data.frame(
  Variable = character(),
  Wilcoxon_Result = character(),
  P_Value = numeric(),
  Median_Vehicle = numeric(),
  Median_VD = numeric(),
  ID = character(),  # Add ID to track which dataset each row belongs to
  stringsAsFactors = FALSE
)

# List of variables to compare
variables_to_compare <- c("Cellage", "GenAge", "SENmayo")

# Iterate through each ID and variable
# Iterate through each ID and variable
for (id_name in table_id$id) {
  # Dynamically access the combined data object for each ID
  data_object <- get(id_name)
  
  for (variable in variables_to_compare) {
    # Check if the variable exists in the data frame
    if (!variable %in% names(data_object@meta.data)) {
      cat("Variable", variable, "not found in", id_name, "\n")
      next  # Skip this iteration
    }
    
    # Perform Wilcoxon test on the current variable using the formula interface
    # Here `reformulate()` dynamically creates the formula for the test
    wilcox_result <- wilcox.test(reformulate("experiment", response = variable), data = data_object@meta.data)
    
    # Convert Wilcoxon statistic if numeric
    wilcox_statistic <- ifelse(is.numeric(wilcox_result$statistic), as.character(wilcox_result$statistic), NA)
    
    # Calculate median for each group
    median_vehicle <- median(data_object@meta.data[data_object@meta.data$experiment == "vehicle", variable], na.rm = TRUE)
    median_vd <- median(data_object@meta.data[data_object@meta.data$experiment == "vd", variable], na.rm = TRUE)
    
    # Create a row for the result dataframe
    result_row <- data.frame(
      Variable = variable,
      Wilcoxon_Result = wilcox_statistic,
      P_Value = wilcox_result$p.value,
      Median_Vehicle = median_vehicle,
      Median_VD = median_vd,
      ID = id_name,  # Record the ID
      stringsAsFactors = FALSE
    )
    
    # Add the row to the result dataframe
    significance <- bind_rows(significance, result_row)
  }
}
write_csv(significance, paste0(plots.dir, "/", "supplementary_table_1.csv"))


#  2.6 Supplementary table 2 wilcox cells  ####

table_id <- data.frame(id = c("aem","lin","hanai"))
### here begin the loop for each value of table_id (id will be the variable name in the loop)

for (id_name in table_id$id) {
  id <- get(id_name)
vehicle_data <- id@meta.data %>%
  filter(experiment == "vehicle")
vd_data <- id@meta.data %>%
  filter(experiment == "vd")

# Get unique labels
unique_labels <- unique(id@meta.data$SingleR.labels)
# Initialize a data frame to store the results
wilcox_results <- data.frame(
  dataset = character(),
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
        dataset = id_name,
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
assign(paste0(id_name, "_results"), wilcox_results, envir = .GlobalEnv)
write_csv(wilcox_results, paste0(plots.dir, "/", id_name, "_celltypes_ES_results.csv"))
}
 
sup_tab_2 <- rbind(lin_results, hanai_results)
sup_tab_2 <- rbind(sup_tab_2, aem_results)
write_csv(sup_tab_2, paste0(plots.dir, "/", "supplementary_table_2.csv"))

#  2.7 clean ####
rm(list=ls())
gc()

# 3. GENE REGULATORY NETWORKS ####
#  3.1 Values ####

plots.dir = paste0(getwd(),"/output/plots/pub/fig3")
input.dir = paste0(getwd(), "/output/seurat/grn/")

#  3.2 libraries ####

library(Seurat)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
library(biomaRt)
library(data.table)
library(UCell)
library(ggplot2)
library(ComplexUpset)
library(enrichR)
library(GeneOverlap)
library(patchwork)

theme_set(theme_cowplot())
set.seed(12345)
 
#  3.3 functions ####
#for legend plots (volcano plot with legend)
create_volcano_plot_legends <- function(group_name, DMEs_all) {
  # Subset data for the specific group
  subset_data <- subset(DMEs_all, group == group_name)
  
  # Map module colors based on the 'module' column values
  module_colors <- setNames(unique(subset_data$module), unique(subset_data$module))
  
  # Create volcano plot
  volcano_plot <- ggplot(subset_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = module_colors, name = "module") +  # Disable legend for individual plots
    labs(title = paste("Volcano Plot for", group_name, "Cell Subtype"),
         x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    geom_rect(data = subset(subset_data, p_val_adj < 0.05), 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
              fill = "grey", alpha = 0.2, color = NA) +
    theme(plot.margin = margin(10, 10, 10, 10, "mm")) + 
    ggplot2::theme(legend.position = "top")
  
  return(volcano_plot)
}
#for volcano plots (volcano plot without legends)
create_volcano_plot <- function(group_name, DMEs_all) {
  # Subset data for the specific group
  subset_data <- subset(DMEs_all, group == group_name)
  
  # Map module colors based on the 'module' column values
  module_colors <- setNames(unique(subset_data$module), unique(subset_data$module))
  
  # Create volcano plot
  volcano_plot <- ggplot(subset_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = module_colors, name = "Cell Subtype", guide = "none") +  # Disable legend for individual plots
    labs(title = paste(group_name),
         x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    geom_rect(data = subset(subset_data, p_val_adj < 0.05), 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
              fill = "grey", alpha = 0.2, color = NA) +
    theme(plot.margin = margin(10, 10, 10, 10, "mm"))
  
  return(volcano_plot)
}
#for getting the module genes 
preupset <- function(data, column_name) {
  unique_modules <- data %>%
    dplyr::select({{column_name}}) %>%
    distinct() %>%
    pull() %>% 
    as.character()
  
  module_genes <- lapply(unique_modules, function(single_module) {
    data %>%
      filter(module == single_module) %>%
      pull(gene)
  })
  
  module_list <- setNames(module_genes, unique_modules)
  
  return(module_list)
}

#  3.4 load ####

aem <- readRDS(paste0(input.dir, "aem_GRN_pseudobulk.rds"))
hanai <- readRDS(paste0(input.dir, "hanai_GRN_pseudobulk.rds"))
lin <- readRDS(paste0(input.dir, "lin_GRN_pseudobulk.rds"))
cellage <- read_csv("input/scores/cellage.csv")
genage <- read_tsv("input/scores/geneage.tsv")
senmayo_mouse <- read_csv("input/scores/senmayo_mouse.csv")

#  3.5 Process loads ####

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gs_cellage <- getLDS(attributes = c("external_gene_name"),
                     filters = "external_gene_name", values = unique(cellage$`Gene symbol`), mart = human,
                     attributesL = c("external_gene_name"), martL = mouse)
rm(human,mouse, cellage)

gs_cellage <- gs_cellage[,2]
gs_geneage <- unique(genage$`Gene Symbol`)
gs_senmayo <- unique(senmayo_mouse$`Gene(murine)`)
rm(cellage, genage, senmayo_mouse)

#  3.6 Figure 3 ####
#   3.6.1 Figure 3A ####
group1 <- aem@meta.data %>% subset(experiment == "vehicle") %>% rownames
group2 <- aem@meta.data %>% subset(experiment == "vd") %>% rownames

group.by = 'SIngleR.labels'
cell_count_by_group <- table(aem@meta.data$SingleR.labels)
valid_groups <- names(cell_count_by_group[cell_count_by_group >= 70])
aem <- aem[, aem@meta.data$SingleR.labels %in% valid_groups]
DMEs_all <- FindAllDMEs(
  aem,
  group.by = 'SingleR.labels',
  wgcna_name = 'pseudobulk'
)
DMEs_all_aem <- DMEs_all
#volcano plot AEM (figure 2A)
# Get unique group names
unique_groups <- unique(DMEs_all$group)

# Create a list to store individual plots
plot_list <- list()

# Get the legend
legend_plot <- create_volcano_plot_legends(unique_groups[[1]], DMEs_all)
legend_plot <- cowplot::get_plot_component(legend_plot, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(legend_plot)

# Iterate over unique group names and create plots
for (group_name in unique_groups) {
  plot_list[[group_name]] <- create_volcano_plot(group_name, DMEs_all)
}

# Add the legend plot to the list
plot_list[["Legend"]] <- legend_plot
# Combine plots and arrange them in a grid layout
fig_3a <- plot_grid(plotlist = plot_list, ncol = 3)

#   3.6.2 Figure 3B ####
# HANAI DME's
group1 <- hanai@meta.data %>% subset(experiment == "vehicle") %>% rownames
group2 <- hanai@meta.data %>% subset(experiment == "vd") %>% rownames

group.by = 'SIngleR.labels'
cell_count_by_group <- table(hanai@meta.data$SingleR.labels)
valid_groups <- names(cell_count_by_group[cell_count_by_group >= 70])
hanai <- hanai[, hanai@meta.data$SingleR.labels %in% valid_groups]
DMEs_all <- FindAllDMEs(
  hanai,
  group.by = 'SingleR.labels',
  wgcna_name = 'pseudobulk'
)
DMEs_all_hanai <- DMEs_all
#Volcano plot HANAI

# Get unique group names
unique_groups <- unique(DMEs_all$group)

# Create a list to store individual plots
plot_list <- list()

# Get the legend
legend_plot <- create_volcano_plot_legends(unique_groups[[1]], DMEs_all)
legend_plot <- cowplot::get_plot_component(legend_plot, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(legend_plot)

# Iterate over unique group names and create plots
for (group_name in unique_groups) {
  plot_list[[group_name]] <- create_volcano_plot(group_name, DMEs_all)
}

# Add the legend plot to the list
plot_list[["Legend"]] <- legend_plot
# Combine plots and arrange them in a grid layout
fig_3b <- plot_grid(plotlist = plot_list, ncol = 3)

#   3.6.3 Figure 3C ####
# LIN DME's
group1 <- lin@meta.data %>% subset(experiment == "vehicle") %>% rownames
group2 <- lin@meta.data %>% subset(experiment == "vd") %>% rownames

group.by = 'SIngleR.labels'
cell_count_by_group <- table(lin@meta.data$SingleR.labels)
valid_groups <- names(cell_count_by_group[cell_count_by_group >= 70])
lin <- lin[, lin@meta.data$SingleR.labels %in% valid_groups]
DMEs_all <- FindAllDMEs(
  lin,
  group.by = 'SingleR.labels',
  wgcna_name = 'pseudobulk'
)
DMEs_all_lin <- DMEs_all
#Volcano plot lin

# Get unique group names
unique_groups <- unique(DMEs_all$group)

# Create a list to store individual plots
plot_list <- list()

# Get the legend
legend_plot <- create_volcano_plot_legends(unique_groups[[1]], DMEs_all)
legend_plot <- cowplot::get_plot_component(legend_plot, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(legend_plot)

# Iterate over unique group names and create plots
for (group_name in unique_groups) {
  plot_list[[group_name]] <- create_volcano_plot(group_name, DMEs_all)
}

# Add the legend plot to the list
plot_list[["Legend"]] <- legend_plot
# Combine plots and arrange them in a grid layout
fig_3c <- plot_grid(plotlist = plot_list, ncol = 3)

#   3.6.4 Figure 3D ####
#    3.6.4.1 Upset plot aem ####
module_umap_df <- as_tibble(aem@misc[["pseudobulk"]][["module_umap"]])
all_genes <- unique(c(module_umap_df$gene, gs_cellage, gs_geneage, gs_senmayo))
#get the module genes
module_list <- preupset(module_umap_df, "module")

# Create a tibble with the modules as columns and genes as rows
gzuschrist <- tibble(
  Genes = all_genes,
  kME = module_umap_df$kME[match(all_genes, module_umap_df$gene)]
)

# Add columns for each module
gzuschrist <- gzuschrist %>%
  bind_cols(as_tibble(lapply(module_list, function(module_genes) as.integer(all_genes %in% module_genes))))

# Add additional columns as needed (e.g., cellage, geneage, senmayo)
gzuschrist <- gzuschrist %>%
  mutate(
    cellage = as.integer(all_genes %in% gs_cellage),
    geneage = as.integer(all_genes %in% gs_geneage),
    senmayo = as.integer(all_genes %in% gs_senmayo)
  )

#get the gene sets for intersections in the upset plot
gene_sets = colnames(gzuschrist)[3:length(gzuschrist)]
module_umap_df$module <- as.character(module_umap_df$module)

# Create a tibble with the sets and save if they are gene modules or senescent gene sets
gs_groups <- tibble(
  set = unique(c(module_umap_df$module, "cellage", "geneage", "senmayo")),
  gene_sets_groups = ifelse(set %in% unique(module_umap_df$module), "module", "senescence")
)

#get the raw upset plot to get the intra senescence intersections and remove them
hm <- ComplexUpset::upset(
  gzuschrist,
  gene_sets, 
  name='gene_sets', 
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=20,
  sort_sets=FALSE, 
  min_degree = 2, 
  mode = "inclusive_intersection",
  annotations = list(
    'kME'=list(
      aes=aes(x=intersection, y=kME),
      geom=geom_boxplot(na.rm=TRUE)
    )),
  stripes=upset_stripes(
    mapping=aes(color = gene_sets_groups),
    colors=c(
      'module'='#ADD8E6',
      'senescence'='#ff6961'
    ),
    data=gs_groups,
  )
)

# Create a tibble with intersections
intersections_tbl <- tibble(intersection = levels(hm[[6]][["data"]][["intersection"]]))
groups <- as_tibble(hm[[6]][["layers"]][[1]][["data"]]) 
groups <- groups %>%
  mutate(group = as.integer(group))


# Reverse the order of the intersection levels
intersections_tbl <- intersections_tbl %>%
  mutate(intersection = rev(intersection))

# Separate the intersections into two groups
intersections_tbl <- intersections_tbl %>%
  separate(intersection, into = c("group1", "group2"), sep = "-", convert = TRUE)

# Join to get the gene_sets_groups
intersections_tbl <- intersections_tbl %>%
  left_join(groups, by = c("group1" = "group")) %>%
  dplyr::rename(gene_sets_group1 = gene_sets_groups) %>%
  left_join(groups, by = c("group2" = "group")) %>%
  dplyr::rename(gene_sets_group2 = gene_sets_groups)

# Filter out intersections where both groups are "senescence"
intersections_tbl <- intersections_tbl %>%
  filter(!(gene_sets_group1 == "senescence" & gene_sets_group2 == "senescence"))

# Retrieve the top 10 intersections from the reversed order
intersections_tbl <- intersections_tbl %>%
  slice_head(n = 10)

groups_y <- intersections_tbl$group_name.y
groups_x <- intersections_tbl$group_name.x

conditions <- mapply(function(y, x) {
  paste0("(", y, " == 1 & ", x, " == 1)")
}, y = groups_y, x = groups_x)

# Combine all conditions into one string with a single separator
combined_conditions <- paste(conditions, collapse = " | ")

# Apply the filter dynamically using eval and parse
gzuschrists <- gzuschrist %>%
  filter(eval(parse(text = combined_conditions)))

# get the official upset plot with top 10 intersections
p <-ComplexUpset::upset(
  gzuschrists,
  gene_sets, 
  name='gene_sets',
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=10,
  sort_sets=FALSE, 
  min_degree = 2, 
  mode = "inclusive_intersection",
  annotations = list(
    'kME'=list(
      aes=aes(x=intersection, y=kME),
      geom=geom_boxplot(na.rm=TRUE)
    )),
  stripes=upset_stripes(
    mapping=aes(color = gene_sets_groups),
    colors=c(
      'module'='#ADD8E6',
      'senescence'='#ff6961'
    ),
    data=gs_groups,
  )
)

#get the legend and update the plot to have no legends:

legend_p <- cowplot::get_plot_component(p[[6]], 'guide-box-right', return_all = TRUE)
legend_p <- ggdraw(legend_p)
p[[6]] <- p[[6]] + theme(legend.position = "none")

#    3.6.4.2 GO/intersect aem ####
#for now ill go with overlaping the intersection genes with the enrichment of the WHOLE module, but see code for enrichment of intersection genes if otherwise decided
#get the intersections as a list for GO 
intersections <- intersections_tbl %>%
  dplyr::select(group_name.y, group_name.x) %>%
  apply(1, as.list) %>%
  lapply(unlist)

# enrichr databases to test
dbs <- c('GO_Biological_Process_2023')

# perform enrichment tests
aem <- RunEnrichr(
  aem,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(aem)

#### Then we need a dataframe of the genes and their kME of each intersection

final_result <- list()
for (intersection in intersections) {
  # Extract the column names from the intersection
  column_names <- intersection
  
  # Create a dynamic filtering condition
  filter_condition <- paste0("gzuschrist$", paste(column_names, collapse = " == 1 & gzuschrist$"))
  
  # Evaluate the filtering condition
  filtered_rows <- which(eval(parse(text = filter_condition)))
  
  # Extract relevant information from the filtered rows
  result_data <- gzuschrist[filtered_rows, c("Genes", "kME")]
  
  # Create a new column with the intersection name
  result_data$Intersection <- paste(intersection, collapse = "-")
  
  # Append the result to the list
  final_result[[length(final_result) + 1]] <- result_data
}
final_result <- bind_rows(final_result)

#### now we'll make a dataframe that sees the overlaps between biological proccesses and intersections

#get both dataframes that we need
enrich_intersection <- as_tibble(enrich_df) %>% 
  dplyr::select(Term, Combined.Score, Genes, Adjusted.P.value)

final_result <- final_result %>%
  mutate(Genes = toupper(Genes))

# Create a list to store the results
result_list <- list()

# Loop through each intersection in final_result
for (intersection in unique(final_result$Intersection)) {
  
  # Filter genes and kME values for the current intersection
  current_intersection <- final_result %>% 
    filter(Intersection == intersection) 
  
  # Extract unique genes for the current intersection
  unique_genes <- unique(current_intersection$Genes)
  
  # Separate rows for each gene in enrich_intersection
  enrich_intersection_separated <- enrich_intersection %>%
    separate_rows(Genes, sep = ";")
  
  # Filter based on genes in the current intersection
  current_result <- enrich_intersection_separated %>% 
    filter(Genes %in% toupper(unique_genes))
  
  # Add the result to the list
  result_list[[intersection]] <- current_result
}

#### now well get the heatmap for visualizing it 

combined_results <- bind_rows(result_list, .id = "Intersection")

# Step 1: Filter terms with adjusted p-values ≤ 0.05 and get their top 3 combined.scores for each intersection
filtered_results <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(
    Combined.Score = max(Combined.Score),
    Adjusted.P.value = Adjusted.P.value[which.max(Combined.Score)]
  ) %>%
  filter(Adjusted.P.value <= 0.05) %>%
  group_by(Intersection) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  mutate(P.value.Indicator = "0.05 or less") %>% 
  ungroup()

# Step 2: Identify intersections where all terms have adjusted p-values > 0.05 and retain their highest combined score terms
remaining_intersections <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(
    Combined.Score = max(Combined.Score),
    Adjusted.P.value = Adjusted.P.value[which.max(Combined.Score)]
  ) %>%
  group_by(Intersection) %>%
  filter(all(Adjusted.P.value > 0.05)) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  mutate(P.value.Indicator = "Greater than 0.05") %>% 
  ungroup()
# Step 3: Add a highlight column and combine the results

combined_results <- bind_rows(filtered_results, remaining_intersections) %>%
  arrange(Intersection, desc(Combined.Score))

# Ensure the Intersections in combined_results follow the order in final_result
ordered_intersections <- unique(final_result$Intersection)
combined_results$Intersection <- factor(combined_results$Intersection, levels = ordered_intersections)

# Create a heatmap with swapped axes and log scale for Combined.Score with legend
h <- ggplot(combined_results, aes(x = Intersection, y = Term)) +
  geom_tile(aes(fill = ifelse(P.value.Indicator == "Greater than 0.05", NA, log(Combined.Score))), color = "white") +
  scale_fill_viridis_c(na.value = "grey", name = "log(Combined.Score)") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Biological Process") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        legend.position = "right")

#get the legend and the y axis whole
legend_h <- cowplot::get_plot_component(h, 'guide-box-right', return_all = TRUE)
y_axis_grob <- cowplot::get_plot_component(h, "axis-l")
y_axis_title <- cowplot::get_plot_component(h + theme(axis.title.y = element_text()), "ylab")
y_axis_full_grob <- cowplot::ggdraw(y_axis_title) + cowplot::ggdraw(y_axis_grob)

#create the heatmap without the legend
h <- ggplot(combined_results, aes(x = Intersection, y = Term)) +
  geom_tile(aes(fill = ifelse(P.value.Indicator == "Greater than 0.05", NA, log(Combined.Score))), color = "white") +
  scale_fill_viridis_c(na.value = "grey", name = "log(Combined.Score)") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

#    3.6.4.3 barplot DME's aem ####

median_avg_log2FC <- DMEs_all_aem %>%
  filter(p_val_adj <= 0.05) %>% 
  group_by(module) %>%
  summarise(median_avg_log2FC = median(avg_log2FC, na.rm = TRUE),
            p_adj_value = median(p_val_adj))

barplot_data <- median_avg_log2FC
y_labs <- ggplot_build(p[[6]])$layout$panel_params[[1]]$y$get_labels()
y_labs <- unname(y_labs)
# Filter barplot_data to include only the relevant modules
filtered_barplot_data <- barplot_data[barplot_data$module %in% y_labs, ]
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = rev(y_labs))
senmayo_row <- data.frame(module = "senmayo", median_avg_log2FC = 0, p_adj_value = 1)
cellage_row <- data.frame(module = "cellage", median_avg_log2FC = 0, p_adj_value = 1)

# Bind rows to the filtered_barplot_data dataframe
filtered_barplot_data <- bind_rows(senmayo_row, cellage_row, filtered_barplot_data)

# Ensure the module column is a factor with the correct levels
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = y_labs)
filtered_barplot_data <- filtered_barplot_data[complete.cases(filtered_barplot_data$module), ]

# Create the barplot with filtered data
barplot <- ggplot(filtered_barplot_data, aes(x = median_avg_log2FC, y = module)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(data = filtered_barplot_data %>% filter(p_adj_value < 0.05),
            aes(x = median_avg_log2FC, y = module, label = "*"),
            position = position_stack(vjust = 0.5), size = 8, color = "red") +  # Add asterisks for p_val > 0.05
  labs(title = NULL,
       x = "Median Avg_Log2FC_all",
       y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank())

#    3.6.4.4 Combined plot aem ####

fig_3d <- y_axis_full_grob  / legend_h / legend_p / barplot | h / p[[2]] /p[[4]] / p[[6]]



#   3.6.5 Figure 3E ####
#    3.6.5.1 Upset plot hanai ####
module_umap_df <- as_tibble(hanai@misc[["pseudobulk"]][["module_umap"]])
all_genes <- unique(c(module_umap_df$gene, gs_cellage, gs_geneage, gs_senmayo))
#get the module genes
module_list <- preupset(module_umap_df, "module")

# Create a tibble with the modules as columns and genes as rows
gzuschrist <- tibble(
  Genes = all_genes,
  kME = module_umap_df$kME[match(all_genes, module_umap_df$gene)]
)

# Add columns for each module
gzuschrist <- gzuschrist %>%
  bind_cols(as_tibble(lapply(module_list, function(module_genes) as.integer(all_genes %in% module_genes))))

# Add additional columns as needed (e.g., cellage, geneage, senmayo)
gzuschrist <- gzuschrist %>%
  mutate(
    cellage = as.integer(all_genes %in% gs_cellage),
    geneage = as.integer(all_genes %in% gs_geneage),
    senmayo = as.integer(all_genes %in% gs_senmayo)
  )

#get the gene sets for intersections in the upset plot
gene_sets = colnames(gzuschrist)[3:length(gzuschrist)]
module_umap_df$module <- as.character(module_umap_df$module)

# Create a tibble with the sets and save if they are gene modules or senescent gene sets
gs_groups <- tibble(
  set = unique(c(module_umap_df$module, "cellage", "geneage", "senmayo")),
  gene_sets_groups = ifelse(set %in% unique(module_umap_df$module), "module", "senescence")
)

#get the raw upset plot to get the intra senescence intersections and remove them
hm <- ComplexUpset::upset(
  gzuschrist,
  gene_sets, 
  name='gene_sets', 
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=20,
  sort_sets=FALSE, 
  min_degree = 2, 
  mode = "inclusive_intersection",
  annotations = list(
    'kME'=list(
      aes=aes(x=intersection, y=kME),
      geom=geom_boxplot(na.rm=TRUE)
    )),
  stripes=upset_stripes(
    mapping=aes(color = gene_sets_groups),
    colors=c(
      'module'='#ADD8E6',
      'senescence'='#ff6961'
    ),
    data=gs_groups,
  )
)

# Create a tibble with intersections
intersections_tbl <- tibble(intersection = levels(hm[[6]][["data"]][["intersection"]]))
groups <- as_tibble(hm[[6]][["layers"]][[1]][["data"]]) 
groups <- groups %>%
  mutate(group = as.integer(group))


# Reverse the order of the intersection levels
intersections_tbl <- intersections_tbl %>%
  mutate(intersection = rev(intersection))

# Separate the intersections into two groups
intersections_tbl <- intersections_tbl %>%
  separate(intersection, into = c("group1", "group2"), sep = "-", convert = TRUE)

# Join to get the gene_sets_groups
intersections_tbl <- intersections_tbl %>%
  left_join(groups, by = c("group1" = "group")) %>%
  dplyr::rename(gene_sets_group1 = gene_sets_groups) %>%
  left_join(groups, by = c("group2" = "group")) %>%
  dplyr::rename(gene_sets_group2 = gene_sets_groups)

# Filter out intersections where both groups are "senescence"
intersections_tbl <- intersections_tbl %>%
  filter(!(gene_sets_group1 == "senescence" & gene_sets_group2 == "senescence"))

# Retrieve the top 10 intersections from the reversed order
intersections_tbl <- intersections_tbl %>%
  slice_head(n = 10)

groups_y <- intersections_tbl$group_name.y
groups_x <- intersections_tbl$group_name.x

conditions <- mapply(function(y, x) {
  paste0("(", y, " == 1 & ", x, " == 1)")
}, y = groups_y, x = groups_x)

# Combine all conditions into one string with a single separator
combined_conditions <- paste(conditions, collapse = " | ")

# Apply the filter dynamically using eval and parse
gzuschrists <- gzuschrist %>%
  filter(eval(parse(text = combined_conditions)))

# get the official upset plot with top 10 intersections
p <-ComplexUpset::upset(
  gzuschrists,
  gene_sets, 
  name='gene_sets',
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=10,
  sort_sets=FALSE, 
  min_degree = 2, 
  mode = "inclusive_intersection",
  annotations = list(
    'kME'=list(
      aes=aes(x=intersection, y=kME),
      geom=geom_boxplot(na.rm=TRUE)
    )),
  stripes=upset_stripes(
    mapping=aes(color = gene_sets_groups),
    colors=c(
      'module'='#ADD8E6',
      'senescence'='#ff6961'
    ),
    data=gs_groups,
  )
)

#get the legend and update the plot to have no legends:

legend_p <- cowplot::get_plot_component(p[[6]], 'guide-box-right', return_all = TRUE)
legend_p <- ggdraw(legend_p)
p[[6]] <- p[[6]] + theme(legend.position = "none")

#    3.6.5.2 GO/intersect hanai ####
#for now ill go with overlaping the intersection genes with the enrichment of the WHOLE module, but see code for enrichment of intersection genes if otherwise decided
#get the intersections as a list for GO 
intersections <- intersections_tbl %>%
  dplyr::select(group_name.y, group_name.x) %>%
  apply(1, as.list) %>%
  lapply(unlist)

# enrichr databases to test
dbs <- c('GO_Biological_Process_2023')

# perform enrichment tests
hanai <- RunEnrichr(
  hanai,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(hanai)

#### Then we need a dataframe of the genes and their kME of each intersection

final_result <- list()
for (intersection in intersections) {
  # Extract the column names from the intersection
  column_names <- intersection
  
  # Create a dynamic filtering condition
  filter_condition <- paste0("gzuschrist$", paste(column_names, collapse = " == 1 & gzuschrist$"), " == 1")
  
  # Evaluate the filtering condition
  filtered_rows <- which(eval(parse(text = filter_condition)))
  
  # Extract relevant information from the filtered rows
  result_data <- gzuschrist[filtered_rows, c("Genes", "kME")]
  
  # Create a new column with the intersection name
  result_data$Intersection <- paste(intersection, collapse = "-")
  
  # Append the result to the list
  final_result[[length(final_result) + 1]] <- result_data
}
## this to add manually the geneage-cellage-turquoise intersection
filter_condition <- "gzuschrist$geneage == 1 & gzuschrist$turquoise == 1 & gzuschrist$cellage == 1"
filtered_rows <- which(eval(parse(text = filter_condition)))
result_data <- gzuschrist[filtered_rows, c("Genes", "kME")]
result_data$Intersection <- paste("cellage-geneage-turquoise", collapse = "-")
#replace the 6th value since is a fake intersection
final_result[[6]] <- result_data

final_result <- bind_rows(final_result)

#### now we'll make a dataframe that sees the overlaps between biological proccesses and intersections

#get both dataframes that we need
enrich_intersection <- as_tibble(enrich_df) %>% 
  dplyr::select(Term, Combined.Score, Genes, Adjusted.P.value)

final_result <- final_result %>%
  mutate(Genes = toupper(Genes))

# Create a list to store the results
result_list <- list()

# Loop through each intersection in final_result
for (intersection in unique(final_result$Intersection)) {
  
  # Filter genes and kME values for the current intersection
  current_intersection <- final_result %>% 
    filter(Intersection == intersection) 
  
  # Extract unique genes for the current intersection
  unique_genes <- unique(current_intersection$Genes)
  
  # Separate rows for each gene in enrich_intersection
  enrich_intersection_separated <- enrich_intersection %>%
    separate_rows(Genes, sep = ";")
  
  # Filter based on genes in the current intersection
  current_result <- enrich_intersection_separated %>% 
    filter(Genes %in% toupper(unique_genes))
  
  # Add the result to the list
  result_list[[intersection]] <- current_result
}

#### now well get the heatmap for visualizing it 

combined_results <- bind_rows(result_list, .id = "Intersection")

# Step 1: Filter terms with adjusted p-values ≤ 0.05 and get their top 3 combined.scores for each intersection
filtered_results <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(
    Combined.Score = max(Combined.Score),
    Adjusted.P.value = Adjusted.P.value[which.max(Combined.Score)]
  ) %>%
  filter(Adjusted.P.value <= 0.05) %>%
  group_by(Intersection) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  mutate(P.value.Indicator = "0.05 or less") %>% 
  ungroup()

# Step 2: Identify intersections where all terms have adjusted p-values > 0.05 and retain their highest combined score terms
remaining_intersections <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(
    Combined.Score = max(Combined.Score),
    Adjusted.P.value = Adjusted.P.value[which.max(Combined.Score)]
  ) %>%
  group_by(Intersection) %>%
  filter(all(Adjusted.P.value > 0.05)) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  mutate(P.value.Indicator = "Greater than 0.05") %>% 
  ungroup()
# Step 3: Add a highlight column and combine the results

combined_results <- bind_rows(filtered_results, remaining_intersections) %>%
  arrange(Intersection, desc(Combined.Score))

# Ensure the Intersections in combined_results follow the order in final_result
ordered_intersections <- unique(final_result$Intersection)
combined_results$Intersection <- factor(combined_results$Intersection, levels = ordered_intersections)

# Create a heatmap with swapped axes and log scale for Combined.Score with legend
h <- ggplot(combined_results, aes(x = Intersection, y = Term)) +
  geom_tile(aes(fill = ifelse(P.value.Indicator == "Greater than 0.05", NA, log(Combined.Score))), color = "white") +
  scale_fill_viridis_c(na.value = "grey", name = "log(Combined.Score)") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Biological Process") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        legend.position = "right")

#get the legend and the y axis whole
legend_h <- cowplot::get_plot_component(h, 'guide-box-right', return_all = TRUE)
y_axis_grob <- cowplot::get_plot_component(h, "axis-l")
y_axis_title <- cowplot::get_plot_component(h + theme(axis.title.y = element_text()), "ylab")
y_axis_full_grob <- cowplot::ggdraw(y_axis_title) + cowplot::ggdraw(y_axis_grob)

#create the heatmap without the legend
h <- ggplot(combined_results, aes(x = Intersection, y = Term)) +
  geom_tile(aes(fill = ifelse(P.value.Indicator == "Greater than 0.05", NA, log(Combined.Score))), color = "white") +
  scale_fill_viridis_c(na.value = "grey", name = "log(Combined.Score)") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

#    3.6.5.3 barplot DME's hanai ####

median_avg_log2FC <- DMEs_all_hanai %>%
  filter(p_val_adj <= 0.05) %>% 
  group_by(module) %>%
  summarise(median_avg_log2FC = median(avg_log2FC, na.rm = TRUE),
            p_adj_value = median(p_val_adj))

barplot_data <- median_avg_log2FC
y_labs <- ggplot_build(p[[6]])$layout$panel_params[[1]]$y$get_labels()
y_labs <- unname(y_labs)
# Filter barplot_data to include only the relevant modules
filtered_barplot_data <- barplot_data[barplot_data$module %in% y_labs, ]
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = rev(y_labs))
senmayo_row <- data.frame(module = "senmayo", median_avg_log2FC = 0, p_adj_value = 1)
geneage_row <- data.frame(module = "geneage", median_avg_log2FC = 0, p_adj_value = 1)
cellage_row <- data.frame(module = "cellage", median_avg_log2FC = 0, p_adj_value = 1)

# Bind rows to the filtered_barplot_data dataframe
filtered_barplot_data <- bind_rows(senmayo_row,geneage_row, cellage_row, filtered_barplot_data)

# Ensure the module column is a factor with the correct levels
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = y_labs)
filtered_barplot_data <- filtered_barplot_data[complete.cases(filtered_barplot_data$module), ]

# Create the barplot with filtered data
barplot <- ggplot(filtered_barplot_data, aes(x = median_avg_log2FC, y = module)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(data = filtered_barplot_data %>% filter(p_adj_value < 0.05),
            aes(x = median_avg_log2FC, y = module, label = "*"),
            position = position_stack(vjust = 0.5), size = 8, color = "red") +  # Add asterisks for p_val > 0.05
  labs(title = NULL,
       x = "Median Avg_Log2FC_all",
       y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank())


#    3.6.5.4 Combined plot hanai ####

fig_3e <- y_axis_full_grob  / legend_h / legend_p / barplot | h / p[[2]] /p[[4]] / p[[6]]

#   3.6.6 Figure 3F ####
#    3.6.6.1 Upset plot lin ####
module_umap_df <- as_tibble(lin@misc[["pseudobulk"]][["module_umap"]])
all_genes <- unique(c(module_umap_df$gene, gs_cellage, gs_geneage, gs_senmayo))
#get the module genes
module_list <- preupset(module_umap_df, "module")

# Create a tibble with the modules as columns and genes as rows
gzuschrist <- tibble(
  Genes = all_genes,
  kME = module_umap_df$kME[match(all_genes, module_umap_df$gene)]
)

# Add columns for each module
gzuschrist <- gzuschrist %>%
  bind_cols(as_tibble(lapply(module_list, function(module_genes) as.integer(all_genes %in% module_genes))))

# Add additional columns as needed (e.g., cellage, geneage, senmayo)
gzuschrist <- gzuschrist %>%
  mutate(
    cellage = as.integer(all_genes %in% gs_cellage),
    geneage = as.integer(all_genes %in% gs_geneage),
    senmayo = as.integer(all_genes %in% gs_senmayo)
  )

#get the gene sets for intersections in the upset plot
gene_sets = colnames(gzuschrist)[3:length(gzuschrist)]
module_umap_df$module <- as.character(module_umap_df$module)

# Create a tibble with the sets and save if they are gene modules or senescent gene sets
gs_groups <- tibble(
  set = unique(c(module_umap_df$module, "cellage", "geneage", "senmayo")),
  gene_sets_groups = ifelse(set %in% unique(module_umap_df$module), "module", "senescence")
)

#get the raw upset plot to get the intra senescence intersections and remove them
hm <- ComplexUpset::upset(
  gzuschrist,
  gene_sets, 
  name='gene_sets', 
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=20,
  sort_sets=FALSE, 
  min_degree = 2, 
  mode = "inclusive_intersection",
  annotations = list(
    'kME'=list(
      aes=aes(x=intersection, y=kME),
      geom=geom_boxplot(na.rm=TRUE)
    )),
  stripes=upset_stripes(
    mapping=aes(color = gene_sets_groups),
    colors=c(
      'module'='#ADD8E6',
      'senescence'='#ff6961'
    ),
    data=gs_groups,
  )
)

# Create a tibble with intersections
intersections_tbl <- tibble(intersection = levels(hm[[6]][["data"]][["intersection"]]))
groups <- as_tibble(hm[[6]][["layers"]][[1]][["data"]]) 
groups <- groups %>%
  mutate(group = as.integer(group))


# Reverse the order of the intersection levels
intersections_tbl <- intersections_tbl %>%
  mutate(intersection = rev(intersection))

# Separate the intersections into two groups
intersections_tbl <- intersections_tbl %>%
  separate(intersection, into = c("group1", "group2"), sep = "-", convert = TRUE)

# Join to get the gene_sets_groups
intersections_tbl <- intersections_tbl %>%
  left_join(groups, by = c("group1" = "group")) %>%
  dplyr::rename(gene_sets_group1 = gene_sets_groups) %>%
  left_join(groups, by = c("group2" = "group")) %>%
  dplyr::rename(gene_sets_group2 = gene_sets_groups)

# Filter out intersections where both groups are "senescence"
intersections_tbl <- intersections_tbl %>%
  filter(!(gene_sets_group1 == "senescence" & gene_sets_group2 == "senescence"))

# Retrieve the top 10 intersections from the reversed order
intersections_tbl <- intersections_tbl %>%
  slice_head(n = 10)

groups_y <- intersections_tbl$group_name.y
groups_x <- intersections_tbl$group_name.x

conditions <- mapply(function(y, x) {
  paste0("(", y, " == 1 & ", x, " == 1)")
}, y = groups_y, x = groups_x)

# Combine all conditions into one string with a single separator
combined_conditions <- paste(conditions, collapse = " | ")

# Apply the filter dynamically using eval and parse
gzuschrists <- gzuschrist %>%
  filter(eval(parse(text = combined_conditions)))

# get the official upset plot with top 10 intersections
p <-ComplexUpset::upset(
  gzuschrists,
  gene_sets, 
  name='gene_sets',
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=10,
  sort_sets=FALSE, 
  min_degree = 2, 
  mode = "inclusive_intersection",
  annotations = list(
    'kME'=list(
      aes=aes(x=intersection, y=kME),
      geom=geom_boxplot(na.rm=TRUE)
    )),
  stripes=upset_stripes(
    mapping=aes(color = gene_sets_groups),
    colors=c(
      'module'='#ADD8E6',
      'senescence'='#ff6961'
    ),
    data=gs_groups,
  )
)

#get the legend and update the plot to have no legends:

legend_p <- cowplot::get_plot_component(p[[6]], 'guide-box-right', return_all = TRUE)
legend_p <- ggdraw(legend_p)
p[[6]] <- p[[6]] + theme(legend.position = "none")

#    3.6.6.2 GO/intersect lin ####
#for now ill go with overlaping the intersection genes with the enrichment of the WHOLE module, but see code for enrichment of intersection genes if otherwise decided
#get the intersections as a list for GO 
intersections <- intersections_tbl %>%
  dplyr::select(group_name.y, group_name.x) %>%
  apply(1, as.list) %>%
  lapply(unlist)

# enrichr databases to test
dbs <- c('GO_Biological_Process_2023')

# perform enrichment tests
lin <- RunEnrichr(
  lin,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(lin)

#### Then we need a dataframe of the genes and their kME of each intersection

final_result <- list()
for (intersection in intersections) {
  # Extract the column names from the intersection
  column_names <- intersection
  
  # Create a dynamic filtering condition
  filter_condition <- paste0("gzuschrist$", paste(column_names, collapse = " == 1 & gzuschrist$"))
  
  # Evaluate the filtering condition
  filtered_rows <- which(eval(parse(text = filter_condition)))
  
  # Extract relevant information from the filtered rows
  result_data <- gzuschrist[filtered_rows, c("Genes", "kME")]
  
  # Create a new column with the intersection name
  result_data$Intersection <- paste(intersection, collapse = "-")
  
  # Append the result to the list
  final_result[[length(final_result) + 1]] <- result_data
}
final_result <- bind_rows(final_result)

#### now we'll make a dataframe that sees the overlaps between biological proccesses and intersections

#get both dataframes that we need
enrich_intersection <- as_tibble(enrich_df) %>% 
  dplyr::select(Term, Combined.Score, Genes, Adjusted.P.value)

final_result <- final_result %>%
  mutate(Genes = toupper(Genes))

# Create a list to store the results
result_list <- list()

# Loop through each intersection in final_result
for (intersection in unique(final_result$Intersection)) {
  
  # Filter genes and kME values for the current intersection
  current_intersection <- final_result %>% 
    filter(Intersection == intersection) 
  
  # Extract unique genes for the current intersection
  unique_genes <- unique(current_intersection$Genes)
  
  # Separate rows for each gene in enrich_intersection
  enrich_intersection_separated <- enrich_intersection %>%
    separate_rows(Genes, sep = ";")
  
  # Filter based on genes in the current intersection
  current_result <- enrich_intersection_separated %>% 
    filter(Genes %in% toupper(unique_genes))
  
  # Add the result to the list
  result_list[[intersection]] <- current_result
}

#### now well get the heatmap for visualizing it 

combined_results <- bind_rows(result_list, .id = "Intersection")

# Step 1: Filter terms with adjusted p-values ≤ 0.05 and get their top 3 combined.scores for each intersection
filtered_results <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(
    Combined.Score = max(Combined.Score),
    Adjusted.P.value = Adjusted.P.value[which.max(Combined.Score)]
  ) %>%
  filter(Adjusted.P.value <= 0.05) %>%
  group_by(Intersection) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  mutate(P.value.Indicator = "0.05 or less") %>% 
  ungroup()

# Step 2: Identify intersections where all terms have adjusted p-values > 0.05 and retain their highest combined score terms
remaining_intersections <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(
    Combined.Score = max(Combined.Score),
    Adjusted.P.value = Adjusted.P.value[which.max(Combined.Score)]
  ) %>%
  group_by(Intersection) %>%
  filter(all(Adjusted.P.value > 0.05)) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  mutate(P.value.Indicator = "Greater than 0.05") %>% 
  ungroup()
# Step 3: Add a highlight column and combine the results

combined_results <- bind_rows(filtered_results, remaining_intersections) %>%
  arrange(Intersection, desc(Combined.Score))

# Ensure the Intersections in combined_results follow the order in final_result
ordered_intersections <- unique(final_result$Intersection)
combined_results$Intersection <- factor(combined_results$Intersection, levels = ordered_intersections)

# Create a heatmap with swapped axes and log scale for Combined.Score with legend
h <- ggplot(combined_results, aes(x = Intersection, y = Term)) +
  geom_tile(aes(fill = ifelse(P.value.Indicator == "Greater than 0.05", NA, log(Combined.Score))), color = "white") +
  scale_fill_viridis_c(na.value = "grey", name = "log(Combined.Score)") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Biological Process") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        legend.position = "right")

#get the legend and the y axis whole
legend_h <- cowplot::get_plot_component(h, 'guide-box-right', return_all = TRUE)
y_axis_grob <- cowplot::get_plot_component(h, "axis-l")
y_axis_title <- cowplot::get_plot_component(h + theme(axis.title.y = element_text()), "ylab")
y_axis_full_grob <- cowplot::ggdraw(y_axis_title) + cowplot::ggdraw(y_axis_grob)

#create the heatmap without the legend
h <- ggplot(combined_results, aes(x = Intersection, y = Term)) +
  geom_tile(aes(fill = ifelse(P.value.Indicator == "Greater than 0.05", NA, log(Combined.Score))), color = "white") +
  scale_fill_viridis_c(na.value = "grey", name = "log(Combined.Score)") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

#    3.6.6.3 barplot DME's lin ####

median_avg_log2FC <- DMEs_all_lin %>%
  filter(p_val_adj <= 0.05) %>% 
  group_by(module) %>%
  summarise(median_avg_log2FC = median(avg_log2FC, na.rm = TRUE),
            p_adj_value = median(p_val_adj))

barplot_data <- median_avg_log2FC
y_labs <- ggplot_build(p[[6]])$layout$panel_params[[1]]$y$get_labels()
y_labs <- unname(y_labs)
# Filter barplot_data to include only the relevant modules
filtered_barplot_data <- barplot_data[barplot_data$module %in% y_labs, ]
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = rev(y_labs))
senmayo_row <- data.frame(module = "senmayo", median_avg_log2FC = 0, p_adj_value = 1)
geneage_row <- data.frame(module = "geneage", median_avg_log2FC = 0, p_adj_value = 1)
cellage_row <- data.frame(module = "cellage", median_avg_log2FC = 0, p_adj_value = 1)

# Bind rows to the filtered_barplot_data dataframe
filtered_barplot_data <- bind_rows(senmayo_row, geneage_row, cellage_row, filtered_barplot_data)

# Ensure the module column is a factor with the correct levels
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = y_labs)
filtered_barplot_data <- filtered_barplot_data[complete.cases(filtered_barplot_data$module), ]

# Create the barplot with filtered data
barplot <- ggplot(filtered_barplot_data, aes(x = median_avg_log2FC, y = module)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(data = filtered_barplot_data %>% filter(p_adj_value < 0.05),
            aes(x = median_avg_log2FC, y = module, label = "*"),
            position = position_stack(vjust = 0.5), size = 8, color = "red") +  # Add asterisks for p_val > 0.05
  labs(title = NULL,
       x = "Median Avg_Log2FC_all",
       y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank())


#    3.6.6.4 Combined plot lin ####

fig_3f <- y_axis_full_grob  / legend_h / legend_p / barplot | h / p[[2]] /p[[4]] / p[[6]]

#   3.6.7 overrall ####

figure3 <- fig_3a / fig_3d | fig_3b / fig_3e | fig_3c / fig_3f

png(filename = paste0(plots.dir, "/figure_3.png"), width = 100, height = 75, units = "cm", res=600)
figure3
dev.off()
ggsave(file=paste0(plots.dir, "/figure_3.svg"), plot=figure3, width=200, height=150, units = "cm", limitsize = FALSE)
#   3.6.8 clean ####
rm(list=ls())
gc()

# 4. CELLULAR COMMUNICATION ####
#  4.1 libraries ####

library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

#  4.2 sconly ####
#   4.2.1 Figure 4 ####
#   4.2.2 values ####

plots.dir = paste0(getwd(),"/output/plots/pub/fig4")
input.dir = paste0(getwd(), "/output/seurat/cellchat/")

#   4.2.3 load ####

aem.vd <- readRDS(paste0(input.dir, "aem_cellchat_vd_sconly.rds"))
aem.vehicle <-readRDS(paste0(input.dir, "aem_cellchat_vehicle_sconly.rds"))
hanai.vd <- readRDS(paste0(input.dir, "hanai_cellchat_vd_sconly.rds"))
hanai.vehicle <- readRDS(paste0(input.dir, "hanai_cellchat_vehicle_sconly.rds"))
lin.vd <- readRDS(paste0(input.dir, "lin_cellchat_vd_sconly.rds"))
lin.vehicle <- readRDS(paste0(input.dir, "lin_cellchat_vehicle_sconly.rds"))

#   4.2.4 Figure 4A (aem) ####
object.list <- list(vehicle = aem.vehicle, vd = aem.vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
png(filename = file.path(plots.dir, "figure_4a.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4a <- gg1 + gg2
fig_4a
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4a.svg"), plot=fig_4a, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.2.5 FIgure 4D ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, "figure_4d.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4d <- gg1 + gg2
fig_4d
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4d.svg"), plot=fig_4d, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.2.6 Figure 4G ####
fig_4g <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
png(filename = file.path(plots.dir, "figure_4g.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4g
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4g.svg"), plot=fig_4g, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.2.7 Figure 4B (hanai) ####
object.list <- list(vehicle = hanai.vehicle, vd = hanai.vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
png(filename = file.path(plots.dir, "figure_4b.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4b <- gg1 + gg2
fig_4b
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4b.svg"), plot=fig_4b, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.2.8 FIgure 4E ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, "figure_4e.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4e <- gg1 + gg2
fig_4e
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4e.svg"), plot=fig_4e, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.2.9 Figure 4H ####
fig_4h <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
png(filename = file.path(plots.dir, "figure_4h.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4h
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4h.svg"), plot=fig_4h, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.2.10 Figure 4C (lin) ####
object.list <- list(vehicle = lin.vehicle, vd = lin.vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
png(filename = file.path(plots.dir, "figure_4c.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4c <- gg1 + gg2
fig_4c
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4c.svg"), plot=fig_4c, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.2.11 Figure 4F ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, "figure_4f.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4f <- gg1 + gg2
fig_4f
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4f.svg"), plot=fig_4f, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.2.12 Figure 4H ####
fig_4i <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
png(filename = file.path(plots.dir, "figure_4i.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4i
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4i.svg"), plot=fig_4i, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.2.13 overall ####

figure4 <- fig_4a / fig_4d / fig_4g | fig_4b / fig_4e / fig_4h | fig_3c / fig_3f / fig_3i

png(filename = paste0(plots.dir, "/figure_4.png"), width = 100, height = 75, units = "cm", res=600)
figure4
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4.svg"), plot=figure4, width=200, height=150, units = "cm", limitsize = FALSE)
#   4.2.14 clean ####
rm(list=ls())
gc()

#  4.3 scall ####
#   4.3.1 Figure 5 ####

#   4.3.2 values ####
plots.dir = paste0(getwd(),"/output/plots/pub/fig5")
input.dir = paste0(getwd(), "/output/seurat/cellchat/")

#   4.3.3 Load ####
aem.vd <- readRDS(paste0(input.dir, "aem_cellchat_vd_scall.rds"))
aem.vehicle <-readRDS(paste0(input.dir, "aem_cellchat_vehicle_scall.rds"))
hanai.vd <- readRDS(paste0(input.dir, "hanai_cellchat_vd_scall.rds"))
hanai.vehicle <- readRDS(paste0(input.dir, "hanai_cellchat_vehicle_scall.rds"))
lin.vd <- readRDS(paste0(input.dir, "lin_cellchat_vd_scall.rds"))
lin.vehicle <- readRDS(paste0(input.dir, "lin_cellchat_vehicle_scall.rds"))
#   4.3.4 Figure 5A (aem) ####
object.list <- list(vehicle = aem.vehicle, vd = aem.vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
png(filename = file.path(plots.dir, "figure_4a.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4a <- gg1 + gg2
fig_4a
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4a.svg"), plot=fig_4a, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.3.5 FIgure 5D ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, "figure_4d.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4d <- gg1 + gg2
fig_4d
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4d.svg"), plot=fig_4d, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.3.6 Figure 5G ####
fig_4g <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
png(filename = file.path(plots.dir, "figure_4g.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4g
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4g.svg"), plot=fig_4g, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.3.7 Figure 5B (hanai) ####
object.list <- list(vehicle = hanai.vehicle, vd = hanai.vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
png(filename = file.path(plots.dir, "figure_4b.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4b <- gg1 + gg2
fig_4b
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4b.svg"), plot=fig_4b, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.3.8 FIgure 5E ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, "figure_4e.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4e <- gg1 + gg2
fig_4e
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4e.svg"), plot=fig_4e, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.3.9 Figure 5H ####
fig_4h <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
png(filename = file.path(plots.dir, "figure_4h.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4h
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4h.svg"), plot=fig_4h, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.3.10 Figure 5C (lin) ####
object.list <- list(vehicle = lin.vehicle, vd = lin.vd) #VERY IMPORTANT SINCE VEHICLE IS GROUP 1 AND VD IS GROUP 2
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
png(filename = file.path(plots.dir, "figure_4c.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4c <- gg1 + gg2
fig_4c
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4c.svg"), plot=fig_4c, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.3.11 Figure 5F ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
gg1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
gg2 <- rankSimilarity(cellchat, type = "functional")
png(filename = file.path(plots.dir, "figure_4f.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4f <- gg1 + gg2
fig_4f
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4f.svg"), plot=fig_4f, width=1280, height=720, units = "px", limitsize = FALSE)
#   4.3.12 Figure 5H ####
fig_4i <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
png(filename = file.path(plots.dir, "figure_4i.png"),
    width = 16, height = 10, units = "cm", res=300)
fig_4i
dev.off()
ggsave(file=paste0(plots.dir, "/figure_4i.svg"), plot=fig_4i, width=1280, height=720, units = "px", limitsize = FALSE)

#   4.3.13 overall ####
figure5 <- fig_5a / fig_5d / fig_5g | fig_5b / fig_5e / fig_5h | fig_5c / fig_5f / fig_3i

png(filename = paste0(plots.dir, "/figure_5.png"), width = 100, height = 75, units = "cm", res=600)
figure5
dev.off()
ggsave(file=paste0(plots.dir, "/figure_5.svg"), plot=figure5, width=200, height=150, units = "cm", limitsize = FALSE)


#   4.3.14 clean ####
rm(list=ls())
gc()
