####### GGplot lab #######
aem_combined <- readRDS("/datos/sensence/emilio/VD/output/seurat/aem_combined_ES.rds")
umap_data <- readRDS("/datos/sensence/emilio/VD/output/seurat/aem_umap_data")
# Extract UMAP coordinates and metadata
umap_data <- as.data.frame(Embeddings(aem_combined, "umap"))
umap_data$SENmayo <- aem_combined@meta.data$SENmayo
#umap_data$SENmayo <- cut(aem_combined@meta.data$SENmayo, breaks = 100)
umap_data$SingleR.labels <- aem_combined@meta.data$SingleR.labels
umap_data$experiment <- aem_combined@meta.data$experiment

saveRDS(object = umap_data, file = file.path("/datos/sensence/emilio/VD/output/seurat/aem_umap_data"))
######### yessir #########

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

umap_ggplot <- ggplot(umap_data, aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +  # Adding alpha for transparency
  scale_fill_manual(values = my_palette) 

# Display the plot
print(umap_ggplot)


base <- ggplot(umap_data, aes(x = umap_1, y = umap_2, fill = SENmayo)) + 
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) 
base + scale_fill_continuous(breaks = 100, labels = T)

num_steps <- 10  # You can adjust this based on your preference

# Choose the specific labels you want to label
labels_to_label <- unique(umap_data$SingleR.labels)
# Create UMAP plot using ggplot2 with labels for each unique cell subtype
base <- ggplot(umap_data, aes(x = umap_1, y = umap_2, fill = SENmayo)) +
  geom_point(shape = 21, size = 1.5, stroke = 0, alpha = 1) +
  facet_wrap(~experiment)+
  scale_fill_viridis_c(name = "SENmayo") + 
  geom_text_repel(
    data = umap_data %>% distinct(SingleR.labels, .keep_all = TRUE) %>% filter(SingleR.labels %in% labels_to_label),
    aes(label = SingleR.labels),
    size = 3, color = "black", box.padding = 0.5
  )

# Display the plot
print(base)
