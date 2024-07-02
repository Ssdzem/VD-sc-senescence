##### enrichment

input.seurat = "/datos/sensence/emilio/VD/output/seurat/hanai_GRN_pseudobulk.rds"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"

sample.name = "hanai_combined_GRN"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat"
#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/hanai/GRN"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"


# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

#for your crazy upsetplot
library(ggplot2)
library(ComplexUpset)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS(input.seurat)


# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

#Print the results for each module
EnrichrBarPlot(
  seurat_obj,
  outdir = plots.dir, # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

######### Intersections in GO terms ##############

enrich_intersection <- as_tibble(enrich_df) %>% 
  dplyr::select(Term, Combined.Score, Genes)

# Convert all genes in enrich_intersection to uppercase
final_result <- final_result %>%
  mutate(Genes = toupper(Genes))
# Create a list to store the results
result_list <- list()

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

# Display the results
result_list

# Combine the individual dataframes from result_list into one dataframe
combined_results <- bind_rows(result_list, .id = "Intersection")

# Rank the processes within each intersection based on Combined.Score
combined_results <- combined_results %>%
  group_by(Intersection) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 5) %>%
  ungroup()

# Create a heatmap with swapped axes and log scale for Combined.Score
heatmap_plot <- ggplot(combined_results, aes(x = Intersection, y = Term, fill = log(Combined.Score + 1))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +  # Use a color scale of your choice (here, using viridis)
  theme_minimal() +
  labs(title = "Top 5 Biological Processes for Each Intersection",
       x = "Intersection",
       y = "Biological Process") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Show the plot
print(heatmap_plot)





# Combine the individual dataframes from result_list into one dataframe
combined_results <- bind_rows(result_list, .id = "Intersection")

# Rank the processes within each intersection based on Combined.Score
combined_results <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(Combined.Score = max(Combined.Score)) %>%
  group_by(Intersection) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 5) %>%
  ungroup()

# Create a heatmap with swapped axes and log scale for Combined.Score
heatmap_plot <- ggplot(combined_results, aes(x = Intersection, y = Term, fill = log(Combined.Score + 1))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +  # Use a color scale of your choice (here, using viridis)
  theme_minimal() +
  labs(title = "Top 5 Biological Processes for Each Intersection",
       x = "Intersection",
       y = "Biological Process") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Show the plot
png(filename = file.path(plots.dir, paste0(sample.name,"IFE SB2", "_GO_heatmap.png")),
    width = 20, height = 16, units = "in", res=300)
print(heatmap_plot)
dev.off()
