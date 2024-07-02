################## fk remove upset chafa (automatize it) #################
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
  rename(gene_sets_group1 = gene_sets_groups) %>%
  left_join(groups, by = c("group2" = "group")) %>%
  rename(gene_sets_group2 = gene_sets_groups)

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

################# overlap with enrich df vs re-enrich #############
# enrichr databases to test
dbs <- c('KEGG_2019_Mouse','GO_Biological_Process_2023')


exp <- final_result %>% 
  filter(Intersection == "cellage-purple") %>% pull(Genes)

exp_enrich <- enrichr(exp, dbs)
exp_enrich1 <- cbind(exp_enrich)
exp


################# get the long ass y labels of the heatmap as an object #####

plot_build <- ggplot_build(h)
y_labels <- plot_build$layout$panel_params[[1]]$y$labels
y <- ggplot_build(h)$layout$panel_params[[1]]$y$get_labels()

y_axis_grob <- cowplot::get_plot_component(h, "axis-l")
y_axis_title <- cowplot::get_plot_component(h + theme(axis.title.y = element_text()), "ylab")
cowplot::ggdraw(y_axis_grob)
cowplot::ggdraw(y_axis_title)
y_axis_full_grob <- cowplot::ggdraw(y_axis_title) + cowplot::ggdraw(y_axis_grob)
cowplot::ggdraw(y_axis_full_grob)

