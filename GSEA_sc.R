library(msigdbr)
library(escape)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library("radarchart")
library(ggplot2) 

gene.sets1 <- getGeneSets(library = "C5")

gene.sets_df <- data.frame(
  gs_name = names(gene.sets1),
  stringsAsFactors = FALSE)
colnames(gene.sets_df) <- "gs_name"

gs_df <- readLines('/home/bnvlab2/cornea/filtered_cornea_pathways.txt')

View(gs_df)

gs_df <- c("GOBP_POSITIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION", "GOBP_REGULATION_OF_STEM_CELL_DIFFERENTIATION") 



filtered_gs <- gene.sets_df %>%
  filter( gs_name %in% gs_df)

filtered_gene_sets <- gene.sets1[filtered_gs$gs_name]
length(filtered_gene_sets)


ES <- enrichIt(obj = prefinalObj,
               gene.sets = filtered_gene_sets,
               cores = 16,
               groups = 1000)

??enrichIt
# # Save the list object to the specified file
# file_path <- "/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/ES.RData" 
# save(ES, file = file_path)
# 
# ES <- load("/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/ES.RData")

ES_copy <- ES
ES@active.ident  <- ES$annot_clusters

ES2 <- data.frame(prefinalObj[[]], Idents(prefinalObj))
#colnames(ES2)[ncol(ES2)] <- "cluster_annot"
prefinalObj <- AddMetaData(prefinalObj, ES)
View(prefinalObj@meta.data)

typeof(ES2)

desired_order <- c(
  "Endothelium",
  "Limbal fibroblasts",
  "Myofibroblasts",
  "Keratocytes",
  'Limbal Basal',
  "Limbal Suprabasal",
  "TAC",
  "Corneal Basal",
  "Corneal Wing",
  "Corneal Superficial"
)

# Reorder the annot_clusters factor based on the desired order
ES2$detailed_annot <- factor(ES2$detailed_annot, levels = desired_order)
ridgeEnrichment(ES2, gene.set = colnames(ES2[188]), group = 'detailed_annot', add.rug = TRUE) 

write.csv(ES2, "/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/Pathways/ES.csv", row.names = FALSE)

colnames(ES2[188])

start_index <- 183
end_index <- start_index + length(filtered_gene_sets) - 1

# Loop through each pathway number
for (pathway_num in start_index:end_index)  {
  
  # Perform ridge enrichment analysis
  p <- ridgeEnrichment(ES2, gene.set = colnames(ES2)[pathway_num], group = 'detailed_annot', add.rug = TRUE)
  
  # Define the file name based on the column name
  file_name <- paste0( '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/Pathways/', colnames(ES2)[pathway_num], '.png')
  
  # Save the plot as a PNG file
  ggsave(filename = file_name, plot = p, width = 10, height = 7)  # Adjust width and height as needed
}

SaveH5Seurat(prefinalObj,  'prefinalObj_iter4.h5Seurat', overwrite = TRUE)



FeaturePlot(prefinalObj, features = c(colnames(ES2[195])), raster=FALSE, min.cutoff = 17000)


pathway_categories <- read.csv(file = '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/Pathways/pathway_categories.csv', quote = '"')
View(pathway_categories)


pathway_categories <- pathway_categories %>%
  mutate(class_name = ifelse(class_name == "keratocyte", "keratin production", class_name))

View(ES2[181:end_index+1])
View(pathway_categories)
#write.csv(pathway_categories, file = "pathway_categories", row.names = FALSE)
ES2_0 <- ES2[180:end_index+1]
dim(ES2_0)[2]
library(tidyr)
df_long <- ES2_0 %>%
  pivot_longer(cols = 1:55, names_to = "gs_name", values_to = "NES")
View(df_long)

df_join <- df_long %>%
  left_join(pathway_categories, by = "gs_name")

View(df_join)



min_max_values <- df_join %>%
  group_by(cluster_annot) %>%
  summarize(
    min_value = min(NES, na.rm = TRUE),
    max_value = max(NES, na.rm = TRUE),
    .groups = 'drop'
  )

# Step 2: Join the min and max values back to df_join
df_join_with_min_max <- df_join %>%
  left_join(min_max_values, by = c( "cluster_annot"))


View(df_join_with_min_max)
# Step 3: Normalize the NES values
df_join_normalized <- df_join_with_min_max %>%
  mutate(
    NES_normalized = (NES - min_value) / (max_value - min_value)
  )

# Print the result
print(df_join_normalized)

#### Mean for each class_name in cluster_annot

df_join_normalized_aver <- df_join_normalized %>%
  group_by( cluster_annot, class_name) %>%
  summarize(mean_exp = mean(NES_normalized))

print(df_join_normalized_aver)


#### Min-max for each class_name in cluster_annot  
min_max_values1 <- df_join_normalized_aver %>%
  group_by(cluster_annot) %>%
  summarize(
    min_exp = min(mean_exp, na.rm = TRUE),
    max_exp = max(mean_exp, na.rm = TRUE),
    .groups = 'drop'
  )
View(min_max_values1)
# Step 2: Join the min and max values back to the original data frame
df_with_min_max1 <- df_join_normalized_aver %>%
  left_join(min_max_values1, by = c("cluster_annot"))

View(df_with_min_max1)
# Step 3: Normalize the mean_exp values
df_normalized1 <- df_with_min_max1 %>%
  mutate(
    mean_exp_normalized = (mean_exp - min_exp) / (max_exp - min_exp)
  )

# Print the result
print(df_join_normalized_aver)

df_normalized12 <- df_normalized1[c('cluster_annot', 'class_name', 'mean_exp_normalized')]

df_final <- df_normalized12 %>%
  pivot_wider(
    id_cols = cluster_annot,  # This specifies the columns that will stay as identifiers
    names_from = class_name,  # This specifies the column whose unique values will become new column names
    values_from = mean_exp_normalized  # This specifies the column whose values will fill the new columns
  )
View(df_final)
chartJSRadar(df_final, maxScale = 1, showToolTipLabel=TRUE)



View(df_wide)

chartJSRadar(scores, maxScale = 10, showToolTipLabel=TRUE)




