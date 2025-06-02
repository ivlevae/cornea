##################### DE function by cell type by condition ############## 

extract_genes <- function(gene_list, pattern) {
  regex <- sprintf("^%s\\d{1,2}([A-OR-Z]?\\d?)?$", pattern)
  pattern_genes <- grep(regex, gene_list, value = TRUE)
  return(pattern_genes)
}
all_MMP <- extract_genes(rownames(cornea), 'MMP')
all_MMP


cornea$condition_detailed <- as.factor(cornea$condition_detailed)


condition <- 'Keratoconus'
DE_function <- function(condition) {
  
  # Subset to two conditions
  cornea_2cond <- subset(cornea, condition_detailed %in% c('Healthy', condition))
  cornea_2cond$condition_detailed <- droplevels(cornea_2cond$condition_detailed)
  cornea_2cond@active.ident <- as.factor(cornea_2cond$condition_detailed)
  
  # Create table of condition × cluster
  table_counts <- table(cornea_2cond$condition_detailed, cornea_2cond$leiden_annot_V3_main_clusters)
  
  # Find clusters to skip if any condition has ≤ 5 cells
  clusters_to_skip <- colnames(table_counts)[apply(table_counts, 2, function(x) any(x <= 5))]
  clusters_to_check <- setdiff(unique(cornea_2cond$leiden_annot_V3_main_clusters), clusters_to_skip)
  results_mmp <- list()
  
  # Run statistical analysis for each cell type
  for (annot in clusters_to_check) {
    message("Processing cluster: ", annot)
    
    
    subset_celltype <- subset(cornea_2cond, leiden_annot_V3_main_clusters == annot)
    subset_celltype@active.ident <- as.factor(subset_celltype$condition_detailed)
    # Run differential expression for MMP genes
    markers <- FindMarkers(subset_celltype,
                           ident.1 = condition,
                           ident.2 = "Healthy",
                           features = all_MMP ,#c('MMP15'),
                          # test.use = 'MAST',
                           min.pct = 0,
                           logfc.threshold = 0)
    markers$celltype <- annot
    markers$gene <- rownames(markers)
    results_mmp[[annot]] <- markers
  }
  
  results_mmp_df <- dplyr::bind_rows(results_mmp)
  return(results_mmp_df)
}



results_kc$p_val_adj_subset <- p.adjust(results_kc$p_val, method = "BH")


results_kc <- DE_function(condition = 'Keratoconus')
View(results_kc)
results_kc_filtered  <- results_kc %>%
  filter(p_val_adj < 0.05) %>%
  arrange(celltype, avg_log2FC)

write_csv(results_kc_filtered, 'results_kc.csv')

View(results_kc_filtered)

results_lscd <- DE_function(condition = 'Limbal Dysplasia')
result_lscd_filtered  <- results_lscd %>%
  filter(p_val_adj < 0.05) %>%
  arrange(celltype, avg_log2FC)

write_csv(result_lscd_filtered, 'results_lscd.csv')

results_ctc <- DE_function(condition = 'Cataract')
results_ctc_filtered  <- results_ctc %>%
  filter(p_val_adj < 0.05) %>%
  arrange(celltype, avg_log2FC)

write_csv(results_ctc_filtered, 'results_ctc.csv')


# expr <- FetchData(cornea, vars = c("MMP15", 'condition_detailed'))
# t <- 
#   expr %>%
#   group_by(condition_detailed) %>%
#   summarize(n = sum(MMP15 > 0))
# View(t)


################## Heatmap #############################################

condition_order <- c("Healthy",
                     "Limbal Dysplasia",
                     "Keratoconus",
                     "Cataract")

cornea$condition_detailed <- factor(cornea$condition_detailed, levels = condition_order)

# avg_table_check <- AverageExpression(cornea,  assays = 'RNA', features =  c('MMP15'), group.by = c(  'leiden_annot_V3_main_clusters', 'condition_detailed'))
# avg_table_check <- as.data.frame(avg_table_check)
# avg_table_check_t <- as.data.frame(t(avg_table_check))
# 
# head(avg_table_check)
# 
# split_colnames <- strsplit(colnames(avg_table_check), "_")
# cell_group <- sapply(split_colnames, function(x) x[1])
# cell_group
# condition <- sapply(split_colnames, function(x) x[2]) #
# condition
# 
# avg_table_check_t$condition <- condition
# 
# avg_table_check_t$celltype <- cell_group
# 
# 
# head(avg_table_check_t)
# heatmap_data <- avg_table_check_t %>%
#   pivot_wider(names_from = celltype, values_from = V1)
# 
# heatmap_data <- as.data.frame(heatmap_data)
# head(heatmap_data)
# 
# rownames(heatmap_data) <- heatmap_data$condition
# heatmap_data$condition <- NULL
# 
# heatmap_data_mat <- as.matrix(as.data.frame(heatmap_data))
# display_numbers <- matrix(sprintf("%.3f", heatmap_data_mat), nrow = nrow(heatmap_data_mat), ncol = ncol(heatmap_data_mat))
# 
# library(RColorBrewer)
# 
# # Define a color palette (e.g., from RColorBrewer)
# color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# pheatmap(heatmap_data, 
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE, 
#          display_numbers = display_numbers,
#          fontsize_number = 14,
#          cellwidth = 50,         # Adjust cell width (default is 10)
#          cellheight = 50,   
#          main = "Expression Heatmap",
#          #breaks = seq(0.006, 0.072, length.out = 100),
#          color = color_palette
# )



gene_expr_df <- FetchData(cornea, vars = c("MMP15", "condition_detailed", "leiden_annot_V3_main_clusters"), slot = "data")
gene_expr_df$MMP15_nonlog <- expm1(gene_expr_df$MMP15)

expression_summary <- gene_expr_df %>%
  group_by(condition_detailed, leiden_annot_V3_main_clusters) %>%
  summarise(
    mean_expr       = mean(MMP15_nonlog, na.rm = TRUE),
    median_expr     = log1p(median(MMP15_nonlog, na.rm = TRUE)),
    n_cells         = n(),
    n_expr_gt0      = sum(MMP15_nonlog > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pct_expr_gt0 = n_expr_gt0 / n_cells * 100)
# View results
View(expression_summary)

expression_summary <- expression_summary %>%
  mutate(
    imputed_expr = case_when(
      n_cells < 25 & mean_expr == 0 ~ 0,
      n_cells < 25 & mean_expr != 0 ~ NA_real_,
      TRUE ~ mean_expr
    )
  )

View(expression_summary)


heatmap_data_2 <- expression_summary %>%
  pivot_wider(names_from = leiden_annot_V3_main_clusters, values_from = imputed_expr, id_cols  =condition_detailed )

heatmap_data_2 <- as.data.frame(heatmap_data_2)
head(heatmap_data_2)

rownames(heatmap_data_2) <- heatmap_data_2$condition_detailed
heatmap_data_2
heatmap_data_2$condition_detailed <- NULL
cell_types <- colnames(heatmap_data_2)

heatmap_data_2_mat <- as.matrix(as.data.frame(t(heatmap_data_2))) ### add t if neeeded
heatmap_data_2_mat <- heatmap_data_2_mat[ nrow(heatmap_data_2_mat):1,]
head(heatmap_data_2_mat)
heatmap_data_2_mat <- heatmap_data_2_mat[c( "Endothelium", "Limbal fibroblasts", 'Myofibroblasts' , "Keratocytes",   "Limbal Epithelium", "TAC", "Corneal Epithelium"),
                                         c('Healthy', 'Keratoconus', 'Cataract', 'Limbal Dysplasia')]
display_numbers <- matrix(sprintf("%.3f", heatmap_data_2_mat), nrow = nrow(heatmap_data_2_mat), ncol = ncol(heatmap_data_2_mat))

library(RColorBrewer)

# Define a color palette (e.g., from RColorBrewer)
color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
pheatmap(heatmap_data_2_mat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = display_numbers,
         fontsize_number = 14,
         cellwidth = 50,         # Adjust cell width (default is 10)
         cellheight = 50,   
         main = "Expression Heatmap",
         #breaks = seq(0.006, 0.072, length.out = 100),
         color = color_palette
)