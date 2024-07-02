###### Variables ####

input.seurat = "/datos/sensence/emilio/VD/output/seurat/lin_GRN_pseudobulk .rds"
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
png(filename = file.path(plots.dir, paste0(sample.name,"IFE SB2", "_DMEvolcano.png")),
    width = 20, height = 16, units = "in", res=300)
modules <- unique(seurat_myobj@misc[["pseudobulk"]][["module_umap"]][, c("module", "color")])
volcano_df <- merge(DMEs, modules, by = "module")
modules$ID <- seq_len(nrow(modules))
rownames(modules) <- modules$ID
modules <- modules[, -which(names(modules) == "ID")]
modules$module <- gsub("_", "-", modules$module)
merged_data <- merge(DMEs, modules, by = "module")
volcano_plot <- ggplot(merged_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
  geom_point(size = 3) +
  scale_color_manual(values = unique(modules$color)) +
  labs(title = "volcano plot of module expression", x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)")
volcano_plot <- volcano_plot +
  geom_rect(data = subset(merged_data, p_val_adj < 0.05), 
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
            fill = "grey", alpha = 0.2, color = NA)

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
# Function to create volcano plots for each group
create_volcano_plot <- function(group_name) {
  # Subset data for the specific group
  subset_data <- subset(DMEs_all, group == group_name)
  
  # Replace complex 'module' names with simplified ones
  subset_data$module <- as.factor(paste0("SB", seq_along(unique(subset_data$module))))
  module_labels <- setNames(levels(subset_data$module), unique(subset_data$module))
  
  # Map module colors based on the simplified module names
  module_colors <- c("SB1" = "turquoise",
                     "SB2" = "green",
                     "SB3" = "brown",
                     "SB4" = "yellow",
                     "SB5" = "blue",
                     "SB6" = "red")
  
  # Create volcano plot
  volcano_plot <- ggplot(subset_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = module)) +
    geom_point(size = 3, shape = 16) +  # Use shape 16 for solid dots
    scale_color_manual(values = module_colors, name = "Cell Subtype") +
    labs(title = paste("Volcano Plot for", group_name, "Cell Subtype"),
         x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    # Add shading for points with p_val_adj < 0.05
    geom_rect(data = subset(subset_data, p_val_adj < 0.05), 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)),
              fill = "grey", alpha = 0.2, color = NA) +
    # Add a margin around the plot
    theme(plot.margin = margin(10, 10, 10, 10, "mm"))
  
  return(volcano_plot)
}


###### with pseudobulk ############


library(ggplot2)
library(cowplot)

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

# Create and arrange volcano plots using facet_wrap
unique_groups <- unique(DMEs_all$group)
plots <- lapply(unique_groups, create_volcano_plot)

# Arrange plots in three columns using facet_wrap
facet_wrap_cols <- 3
plots_arranged <- cowplot::plot_grid(plotlist = plots, ncol = facet_wrap_cols)

# Save or display the arranged plots
ggsave("volcano_plots_facet_wrap.png", plots_arranged, width = 15, height = 8)  # To save the arranged plots as an image
# or
png(filename = file.path(plots.dir, paste0(sample.name,"IFE SB2", "_DMEallvolcano.png")),
    width = 20, height = 16, units = "in", res=300)
print(plots_arranged)  # To display the arranged plots
dev.off()

####### Overlap with senescence gene sets ########

cellage <- read_delim("scores/cellage3.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

genage <- read_delim("scores/genage_mouse.tsv", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

senmayo_mouse <- read_csv("scores/senmayo_mouse.csv")

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2021.archive.ensembl.org") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://may2021.archive.ensembl.org")

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
          file = file.path(output.dir, paste0(sample.name, "_DMEs IFE SB2.csv")))
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
png(filename = file.path(plots.dir, paste0(sample.name,"IFE SB2", "_UMAPoverlaptop10.png")),
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

png(filename = file.path(plots.dir, paste0(sample.name,"IFE SB2", "_upset_plot.png")),
    width = 20, height = 16, units = "in", res=300)
ComplexUpset::upset(
  gzuschrist,
  gene_sets, 
  name='gene_sets', 
  width_ratio=0.1,
  sort_intersections=FALSE,
  #intersections=list(
  #  c('cellage','SB21'),
  #  c('cellage','SB24'),
  #  c('cellage','SB25'),
  #  c('cellage','SB23'),
  #  c('geneage','SB21'),
  #  c('cellage','SB22'),
  #  c('geneage','cellage','SB21'),
  #  c('senmayo','SB21'),
  #  c('cellage','SB26'),
  #  c('senmayo','SB22')
  #),
  sort_sets=FALSE, 
  min_degree = 2, 
  n_intersections=15,
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

#############################################