###### Instructions ##########

# This code has to get 2 manual inputs, the intersections desired to be preserved
# (line 432) and the variables in the beginning (line 8)

###### Variables ####

input.seurat = "/datos/sensence/emilio/VD/output/seurat/aem_GRN_pseudobulk .rds"
#input.seurat = "results/seurat_raw/500_PBMC_3p_LT_Chromium/average/500_PBMC_3p_LT_Chromium_clusters_seurat_qc_raw_average.Rds"

sample.name = "lin_GRN_pb_dif"
#sample.name = "500_PBMC_3p_LT_Chromium_X_clusters_average"

output.dir = "/datos/sensence/emilio/VD/output/seurat"
#output.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average"

########## consider separating output plots by stellarscope method
plots.dir = "/datos/sensence/emilio/VD/output/plots/lin/GRN"
#plots.dir = "results/seurat_analysis/500_PBMC_3p_LT_Chromium_X/average/plots"

# cell.type = "IFE SB2"
#check Seurat object meta data i dont know why this doesnt work

###### Set up ######

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
# re-load
seurat_myobj <- readRDS(input.seurat)

######### DME analysis comparing two groups #########

group1 <- seurat_myobj@meta.data %>% subset(experiment == "vehicle") %>% rownames
group2 <- seurat_myobj@meta.data %>% subset(experiment == "vd") %>% rownames

DMEs <- FindDMEs(
  seurat_myobj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='pseudobulk'
)

head(DMEs)
write.csv(x = DMEs, 
          file = file.path(output.dir, paste0(sample.name, ".csv")))
#volcano plot from DME of all cells (pseudobulk kinda)
png(filename = file.path(plots.dir, paste0(sample.name, "_DMEvolcano.png")),
    width = 20, height = 16, units = "in", res=300)
modules <- unique(seurat_myobj@misc[["pseudobulk"]][["module_umap"]][, c("module", "color")])
volcano_df <- merge(DMEs, modules, by = "module")
modules$ID <- seq_len(nrow(modules))
rownames(modules) <- modules$ID
modules <- modules[, -which(names(modules) == "ID")]
modules$module <- gsub("_", "-", modules$module)
merged_data <- merge(DMEs, modules, by = "module")
modules <- modules[order(modules$module), ]

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
  geom_point(size = 3) +
  scale_color_manual(values = modules$color) +  # Use sorted colors from modules dataframe
  labs(title = "Volcano plot of module expression", x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)")

# Highlight significant points
volcano_plot <- volcano_plot +
  geom_rect(data = subset(merged_data, p_val_adj < 0.05), 
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
            fill = "grey", alpha = 0.2, color = NA)

volcano_plot
print(volcano_plot)
dev.off()

#Volcano plot of all cell subtypes vs modules

group.by = 'SIngleR.labels'
cell_count_by_group <- table(seurat_myobj@meta.data$SingleR.labels)
valid_groups <- names(cell_count_by_group[cell_count_by_group >= 70])
seurat_filtered <- seurat_myobj[, seurat_myobj@meta.data$SingleR.labels %in% valid_groups]
DMEs_all <- FindAllDMEs(
  seurat_filtered,
  group.by = 'SingleR.labels',
  wgcna_name = 'pseudobulk'
)
write.csv(x = DMEs_all, 
          file = file.path(output.dir, paste0(sample.name, ".csv")))


###### with pseudobulk ############

# Function to create volcano plot
create_volcano_plot_legends <- function(group_name, DMEs_all) {
  # Subset data for the specific group
  subset_data <- subset(DMEs_all, group == group_name)
  
  # Map module colors based on the 'module' column values
  module_colors <- setNames(unique(subset_data$module), unique(subset_data$module))
  
  # Create volcano plot
  volcano_plot <- ggplot(subset_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = module_colors, name = "Cell Subtype") +  # Disable legend for individual plots
    labs(title = paste("Volcano Plot for", group_name, "Cell Subtype"),
         x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    geom_rect(data = subset(subset_data, p_val_adj < 0.05), 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
              fill = "grey", alpha = 0.2, color = NA) +
    theme(plot.margin = margin(10, 10, 10, 10, "mm"))
  
  return(volcano_plot)
}
legend_plot <- cowplot::get_legend(create_volcano_plot_legends('IFE SB2', DMEs_all))

create_volcano_plot <- function(group_name, DMEs_all) {
  # Subset data for the specific group
  subset_data <- subset(DMEs_all, group == group_name)
  
  # Map module colors based on the 'module' column values
  module_colors <- setNames(unique(subset_data$module), unique(subset_data$module))
  
  # Create volcano plot
  volcano_plot <- ggplot(subset_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = module_colors, name = "Cell Subtype", guide = "none") +  # Disable legend for individual plots
    labs(title = paste("Volcano Plot for", group_name, "Cell Subtype"),
         x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    geom_rect(data = subset(subset_data, p_val_adj < 0.05), 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
              fill = "grey", alpha = 0.2, color = NA) +
    theme(plot.margin = margin(10, 10, 10, 10, "mm"))
  
  return(volcano_plot)
}

# Get unique group names
unique_groups <- unique(DMEs_all$group)

# Create a list to store individual plots
plot_list <- list()

# Iterate over unique group names and create plots
for (group_name in unique_groups) {
  plot_list[[group_name]] <- create_volcano_plot(group_name, DMEs_all)
}

# Add the legend plot to the list
plot_list[["Legend"]] <- legend_plot

# Combine plots and arrange them in a grid layout
combined_plot <- plot_grid(plotlist = plot_list, ncol = 3)  # You can adjust 'ncol' based on your preference

# Print or save the combined plot
print(combined_plot)




####################################

png(filename = file.path(plots.dir, paste0(sample.name, "_DMEallvolcano.png")),
    width = 20, height = 16, units = "in", res=300)
print(combined_plot)  # To display the arranged plots
dev.off()

####### Overlap with senescence gene sets ########

cellage <- read_delim("scores/cellage3.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

genage <- read_delim("scores/genage_mouse.tsv", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

senmayo_mouse <- read_csv("scores/senmayo_mouse.csv")

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gs_cellage <- getLDS(attributes = c("external_gene_name"),
                     filters = "external_gene_name", values = unique(cellage$`Gene symbol`), mart = human,
                     attributesL = c("external_gene_name"), martL = mouse)
rm(human,mouse, cellage)

gs_cellage <- gs_cellage[,2]
gs_geneage <- unique(genage$`Gene Symbol`)
gs_senmayo <- unique(senmayo_mouse$`Gene(murine)`)

####### overlap senescence genes code ######

# Your original data frame
module_umap_df <- as.data.frame(seurat_myobj@misc[["pseudobulk"]][["module_umap"]])
#Get the list for the genes in each module
genes_by_module <- split(module_umap_df$gene, module_umap_df$module)
# Remove modules with no genes
genes_by_module <- genes_by_module[sapply(genes_by_module, length) > 0]


#run the functions boi
overlap_counts2 <- lapply(genes_by_module, function(module_genes) {
  overlap_count2 <- list(
    gs_cellage_percentage = (length(intersect(gs_cellage, module_genes))/length(gs_cellage))*100,
    gs_geneage_percentage = (length(intersect(gs_geneage, module_genes))/length(gs_geneage))*100,
    gs_senmayo_percentage = (length(intersect(gs_senmayo, module_genes))/length(gs_senmayo))*100
  )
  return(overlap_count2)
})

# Convert the result list to a data frame
overlap_counts2_df <- as.data.frame(do.call(rbind, overlap_counts2))
overlap_counts2_df <- overlap_counts2_df %>%
  mutate(module = row.names(overlap_counts2_df))
overlap_counts2_df <- cbind(overlap_counts2_df['module'], overlap_counts2_df[, setdiff(names(overlap_counts2_df), 'module')])
overlap_counts2_df <- overlap_counts2_df %>%
  unnest(cols = starts_with("gs_"))

# see which genes overlap in each module 

overlap <- lapply(genes_by_module, function(module_genes) {
  overlap <- list(
    gs_cellage_percentage = intersect(gs_cellage, module_genes),
    gs_geneage_percentage = intersect(gs_geneage, module_genes),
    gs_senmayo_percentage = intersect(gs_senmayo, module_genes)
  )
  return(overlap)
})

# Create a function to get kME values for a given gene
get_kME <- function(gene) {
  kME_value <- module_umap_df$kME[module_umap_df$gene == gene]
  if (length(kME_value) > 0) {
    return(kME_value)
  } else {
    return(NA)
  }
}

# Apply the function to each gene in the overlap list
overlap_with_kME <- lapply(overlap, function(sublist) {
  sublist_with_kME <- lapply(sublist, function(sub) {
    if (is.character(sub)) {
      # If it's a character vector (gene names), get kME values
      kME_values <- sapply(sub, get_kME)
      return(kME_values)
    } else {
      # If it's not a character vector, keep it as is
      return(sub)
    }
  })
  return(sublist_with_kME)
})

create_gene_kME_df <- function(overlap_list) {
  # Initialize an empty list to store individual data frames
  df_list <- list()
  
  # Iterate through each sublist in the provided list
  for (sublist_name in names(overlap_list)) {
    subset_list <- overlap_list[[sublist_name]]
    
    # Determine the lengths of gene sets within the sublist
    cellage_len <- length(subset_list$gs_cellage_percentage)
    geneage_len <- length(subset_list$gs_geneage_percentage)
    senmayo_len <- length(subset_list$gs_senmayo_percentage)
    
    # Combine gene names and kME values
    all_genes <- c(
      names(subset_list$gs_cellage_percentage),
      names(subset_list$gs_geneage_percentage),
      names(subset_list$gs_senmayo_percentage)
    )
    
    all_kME <- c(
      gs_cellage_percentage = as.numeric(subset_list$gs_cellage_percentage),
      gs_geneage_percentage = as.numeric(subset_list$gs_geneage_percentage),
      gs_senmayo_percentage = as.numeric(subset_list$gs_senmayo_percentage)
    )
    
    # Create a data frame for the current sublist
    df <- data.frame(
      genes = all_genes,
      kME = all_kME,
      module = sublist_name,
      gene_set = rep(c("gs_cellage", "gs_geneage", "gs_senmayo"),
                     times = c(cellage_len, geneage_len, senmayo_len))
    )
    
    # Append the data frame to the list
    df_list[[sublist_name]] <- df
  }
  
  # Combine all data frames into a single data frame
  combined_df <- do.call(rbind, df_list)
  
  return(combined_df)
}

overlap_with_kME_df <- create_gene_kME_df(overlap_with_kME)
overlap_with_kME_df %>%
  rownames_to_column(var = "row_id") %>%
  mutate(row_id = row_number()) %>%
  dplyr::select(-row_id)

write.csv(x = overlap_with_kME_df, 
          file = file.path(output.dir, paste0(sample.name, ".csv")))
####### Plot the UMAP network for only the overlap genes #########

seurat_myobj <- RunModuleUMAP(
  seurat_myobj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1, # min distance between points in UMAP space
  genes_use = overlap_with_kME_df$genes
)
umap_df <- GetModuleUMAP(seurat_myobj)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
png(filename = file.path(plots.dir, paste0(sample.name,"_UMAPoverlaptop10.png")),
    width = 20, height = 16, units = "in", res=300)
ModuleUMAPPlot(
  seurat_myobj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=20 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()

################## Upset plot finally boya #################
module_umap_df <- as_tibble(module_umap_df)
all_genes <- unique(c(module_umap_df$gene, gs_cellage, gs_geneage, gs_senmayo))

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
write.csv(x = gzuschrist, 
          file = file.path(output.dir, paste0(sample.name, "_upset.csv")))

gene_sets = colnames(gzuschrist)[3:length(gzuschrist)]
#Start the graphic
module_umap_df$module <- as.character(module_umap_df$module)

# Create a tibble with the sets and their corresponding values
gs_groups <- tibble(
  set = unique(c(module_umap_df$module, "cellage", "geneage", "senmayo")),
  gene_sets_groups = ifelse(set %in% unique(module_umap_df$module), "module", "senescence")
)
png(filename = file.path(plots.dir, paste0(sample.name, "_upsetplotchafa.png")),
    width = 20, height = 16, units = "in", res=300)
ComplexUpset::upset(
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
dev.off()
#create a subset of gzuschrist for only getting the top 10 without intra senescence intersections and its list because will need it in GO enrichment
gzuschrists <- gzuschrist %>%
  filter(  (cellage == 1 & saddlebrown == 1) |
           (cellage == 1 & green == 1) |
           (cellage == 1 & brown == 1) |
           (cellage == 1 & turquoise == 1) |
             (cellage == 1 & purple == 1) |
           (cellage == 1 & yellow == 1) |
           (cellage == 1 & red == 1) |
           (cellage == 1 & greenyellow == 1) |
           (cellage == 1 & tan == 1) |
           (cellage == 1 & darkmagenta == 1) 
           )

intersections <- list(
  c("cellage", "saddlebrown"),
  c("cellage", "green"),
  c("cellage", "brown"),
  c("cellage", "turquoise"),
  c("cellage", "purple"),
  c("cellage", "yellow"),
  c("cellage", "red"),
  c("cellage", "greenyellow"),
  c("cellage", "tan"),
  c("cellage", "darkmagenta")
)
#############################################

png(filename = file.path(plots.dir, paste0(sample.name, "_upsetplot.png")),
    width = 20, height = 16, units = "in", res=300)
ComplexUpset::upset(
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
dev.off()

#rm(preupset, get_kME, create_gene_kME_df, create_volcano_plot, create_volcano_plot_legends, valid_groups, unique_groups, gs_cellage, gs_geneage, gs_senmayo, group1, group2, group.by, group_name, cell_count_by_group, all_genes, volcano_plot, volcano_df, umap_df, seurat_filtered, senmayo_mouse, plot_list, overlap, overlap_counts2, overlap_counts2_df, overlap_with_kME, overlap_with_kME_df, modules, module_umap_df, module_list, merged_data, legend_plot, genes_by_module, genage, DMEs, DMEs_all, combined_plot)
########### GO enrichment with intersections ################

### Here we get the enrichment for each module of the GRN

# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
seurat_myobj <- RunEnrichr(
  seurat_myobj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_myobj)

#Print the results for each module
EnrichrBarPlot(
  seurat_myobj,
  outdir = plots.dir, # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

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
  dplyr::select(Term, Combined.Score, Genes)

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

# Rank the processes within each intersection based on Combined.Score
combined_results <- combined_results %>%
  group_by(Intersection, Term) %>%
  summarise(Combined.Score = max(Combined.Score)) %>%
  group_by(Intersection) %>%
  arrange(desc(Combined.Score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  ungroup()

# Create a heatmap with swapped axes and log scale for Combined.Score
heatmap_plot <- ggplot(combined_results, aes(x = Intersection, y = Term, fill = log(Combined.Score))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +  # Use a color scale of your choice (here, using viridis)
  theme_minimal() +
  labs(title = "Top 3 Biological Processes for Each Intersection",
       x = "Intersection",
       y = "Biological Process") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

png(filename = file.path(plots.dir, paste0(sample.name, "_heatmap_GO.png")),
    width = 20, height = 16, units = "in", res=300)
plot(heatmap_plot)
dev.off()





######### mother of all plots ####################

#upset p
p <- ComplexUpset::upset(
  gzuschrists,
  gene_sets, 
  name='gene_sets',
  set_sizes = F,
  width_ratio=0.1,
  sort_intersections=FALSE,
  n_intersections=8,
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

# Code for heatmap
top_terms <- combined_results %>%
  group_by(Intersection) %>%
  top_n(3, Combined.Score) %>%
  ungroup()
top_terms$Intersection <- gsub("-", "_", top_terms$Intersection)

desired_order <- sapply(intersections, function(x) paste(x, collapse = "_"))
desired_order <- desired_order[1:8]
top_terms$Intersection <- factor(top_terms$Intersection, levels = desired_order)
top_terms <- na.omit(top_terms)
# Plot with ordered x-axis
heatmap_plot <- ggplot(top_terms, aes(x = Intersection, y = Term, fill = log(Combined.Score))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +  
  theme_minimal() +
  labs(title = NULL,
       x = NULL,  # Remove x-axis title
       y = "Biological Process") +
 theme(axis.text.x = element_blank(),
       axis.text.y = element_text(size = 6),
       legend.position = "right") + # Remove x-axis labels
  scale_y_discrete(position = "right")

# Code for barplot
barplot_data <- DMEs %>%
  mutate(module = factor(module, levels = rev(unique(DMEs$module))))
y_labs <- ggplot_build(p[[3]])$layout$panel_params[[1]]$y$get_labels()
y_labs <- unname(y_labs)
# Filter barplot_data to include only the relevant modules
filtered_barplot_data <- barplot_data[barplot_data$module %in% y_labs, ]
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = rev(y_labs))
senmayo_row <- data.frame(module = "senmayo", avg_log2FC = 0, p_val = 1)
cellage_row <- data.frame(module = "cellage", avg_log2FC = 0, p_val = 1)

# Bind rows to the filtered_barplot_data dataframe
filtered_barplot_data <- bind_rows(senmayo_row, cellage_row, filtered_barplot_data)

# Ensure the module column is a factor with the correct levels
filtered_barplot_data$module <- factor(filtered_barplot_data$module, levels = y_labs)
filtered_barplot_data <- filtered_barplot_data[complete.cases(filtered_barplot_data$module), ]

# Create the barplot with filtered data
barplot <- ggplot(filtered_barplot_data, aes(x = avg_log2FC, y = module)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(data = filtered_barplot_data %>% filter(p_val < 0.05),
            aes(x = avg_log2FC, y = module, label = "*"),
            position = position_stack(vjust = 0.5), size = 8, color = "red") +  # Add asterisks for p_val > 0.05
  labs(title = NULL,
       x = "Average Log2FC",
       y = NULL) +
  theme_minimal() +
 theme(axis.text.y = element_blank())
# Combine plots

png(filename = file.path(plots.dir, paste0(sample.name, "_big_plot.png")),
    width = 13, height = 8, units = "in", res=300)
mother_of_plots <- plot_spacer() / plot_spacer() / plot_spacer() / barplot | heatmap_plot / p[[1]] /p[[2]] / p[[3]]
mother_of_plots
dev.off()
