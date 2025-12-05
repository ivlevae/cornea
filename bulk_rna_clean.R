library(edgeR) 
library(limma)

#################################################################################

file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_replication_gene_counts.txt.gz"
data_GSE77938_rep <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

colnames(data_GSE77938_rep) <- paste0(colnames(data_GSE77938_rep), "_1")
metadata_rep <- data.frame(
  row.names = colnames(data_GSE77938_rep),  
  Condition = c( rep("KC", 17), rep("KR", 17)) 
)

metadata_rep$Condition <- as.factor(metadata_rep$Condition)
metadata_rep$Condition <- relevel(metadata_rep$Condition, ref = "KR")
metadata_rep <- metadata_rep
metadata_rep$Batch = 'Batch_1'

head(metadata_rep)

#################################################################################

file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_discovery_gene_counts.txt.gz"
data_GSE77938 <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

colnames(data_GSE77938) <- paste0(colnames(data_GSE77938), "_2")

metadata <- data.frame(
  row.names = colnames(data_GSE77938),  
  Condition = c(rep("KR", 8), rep("KC", 8)) 
)
metadata$Condition <- as.factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "KR")
metadata$Batch <- 'Batch_2'

####################################################################################

file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE151631_Raw.counts.csv.gz"
data_GSE151631 <- read.csv(file_name, header = TRUE)

meta_GSE151631_raw <- read.csv("/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE151631_meta.csv")

metadata_GSE151631 <- data.frame(
  row.names = meta_GSE151631_raw$Sample_ID,  
  Condition = meta_GSE151631_raw$Group  
)
metadata_GSE151631$Condition <- as.factor(metadata_GSE151631$Condition)
metadata_GSE151631$Condition <- relevel(metadata_GSE151631$Condition, ref = "KR")
metadata_GSE151631$Batch <- 'Batch_3'

################## 
head(data_GSE151631)  ### has human Gene id only
head(data_GSE77938_rep)  ### hase ENSEMBL Gene id only


####################
library(tidyverse)

#### There is some duplicates by gene name 

duplicates <- data_GSE151631 |>
  group_by(Gene) |>
  filter(n() > 1) |>
  ungroup()

#### merge and replace with mean rounding to the floor

data_GSE1516311 <- data_GSE151631 |>
  group_by(Gene) |>
  summarise(
    across(where(is.numeric), ~ as.integer(floor(mean(.x, na.rm = TRUE)))), 
    across(where(~ !is.numeric(.)), ~ first(.x)),                           
    .groups = "drop"
  )




### Load ensemlid convertation database 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  
  mart = ensembl
)


################# Replace external_gene_name with ensembl_gene_id if its missing 

gene_mapping_fill <- gene_mapping %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "", ensembl_gene_id, external_gene_name))



# gene_mapping_clean <- gene_mapping_fill %>%
#   group_by(external_gene_name) %>%
#   filter(row_number() == 1) %>%
#   ungroup()


duplicates_in_mapping <- gene_mapping_fill[duplicated(gene_mapping_fill$external_gene_name), ]

################################### Conbine GSE77938 studies ###########


common <- intersect(rownames(data_GSE77938), rownames(data_GSE77938_rep))
length(common)
dim(data_GSE77938)
combined_matrix_2_3 <- cbind(data_GSE77938, data_GSE77938_rep)
head(combined_matrix_2_3)
data_GSE77938 <- cbind(ensembl_gene_id = rownames(combined_matrix_2_3), combined_matrix_2_3)
head(data_GSE77938)



data_GSE77938 <- merge(data_GSE77938, gene_mapping_fill, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

duplicates_in_mapping_2 <- data_GSE77938[duplicated(data_GSE77938$external_gene_name), ]
dim(duplicates_in_mapping_2)

sum(is.na(data_GSE77938$external_gene_name))
data_GSE77938_filtered <- data_GSE77938[!is.na(data_GSE77938$external_gene_name), ]


duplicates_in_mapping_3 <- data_GSE77938_filtered %>%
  group_by(external_gene_name) %>%
  filter(n() > 1) %>%
  ungroup()

head(data_GSE77938_filtered)

data_GSE77938_filtered2 <- data_GSE77938_filtered |>
  group_by(external_gene_name) |>
  summarise(
    across(where(is.numeric), ~ as.integer(floor(mean(.x, na.rm = TRUE)))), 
    across(where(~ !is.numeric(.)), ~ first(.x)),                           
    .groups = "drop"
  )

dim(data_GSE77938_filtered)
dim(data_GSE77938_filtered2)
dim(duplicates_in_mapping_3)

View(duplicates_in_mapping_3)

dim(data_GSE77938)
dim(data_GSE1516311)


# 
head(data_GSE77938_filtered2)

head(data_GSE1516311)

data_all <- merge(data_GSE1516311, data_GSE77938_filtered2, by.x = "Gene", by.y = "external_gene_name", all.y = TRUE)
dim(data_GSE151631)
dim(data_GSE77938_filtered2)
dim(data_all)

head(data_all)
rownames(data_all) <- data_all$ensembl_gene_id
data_all$Gene <- NULL
data_all$ensembl_gene_id <- NULL
#sum(is.na(data_all$ensembl_gene_id))

combined_metadata2 <- rbind(metadata, metadata_rep, metadata_GSE151631)

combined_metadata3 <- combined_metadata2[match(colnames(data_all), rownames(combined_metadata2)), ]


data_all[is.na(data_all)] <- 0

# (recommended for edgeR) ensure integer matrix
data_all <- as.matrix(data_all)


dge <- DGEList(counts = data_all)
dge <- calcNormFactors(dge)  # Normalize

dim(combined_metadata)
dim(dge)
combined_metadata <- combined_metadata3
# Use voom to estimate the mean-variance relationship and log2 transformation
v <- voom(dge, design = model.matrix(~ Condition + Batch, data = combined_metadata))

# Apply removeBatchEffect to remove batch effects
v_no_batch <- removeBatchEffect(v$E, batch = combined_metadata$Batch, design = model.matrix(~ Condition, data = combined_metadata))

# Differential expression analysis using limma
fit <- lmFit(v_no_batch, design = model.matrix(~ Condition, data = combined_metadata))
fit <- eBayes(fit)

colnames(fit$coefficients)
# Get the results
results <- topTable(fit, coef = "ConditionKC", number = Inf)

# View results
View(results)

results_df <- as.data.frame(results)
results_df$EnsemblGeneID <- rownames(results_df)
ensembl_ids <- results_df$EnsemblGeneID

dim(results_df)
data_with_hugo <- merge(results_df, gene_mapping, by.x = "EnsemblGeneID", by.y = "ensembl_gene_id", all.x = TRUE)
dim(data_with_hugo)

mmp_all <- gene_mapping[grepl("^MMP\\d{1,2}([A-OR-Z]?\\d?)?$", gene_mapping$external_gene_name, ignore.case = TRUE), ]$external_gene_name
mmp_all

mmp_all_all <- subset( data_with_hugo,  external_gene_name %in% mmp_all)
#^MMP\\d{1,2}$", data_with_hugo$external_gene_name, ignore.case = TRUE), ]

View(mmp_all_all)

all_filtered <- mmp_all_all %>%
  arrange('log2FoldChange')
# View the filtered rows
View(all_filtered)

df <- all_filtered

# Set gene names as rownames
rownames(df) <- df$external_gene_name

# Create the logFC matrix for pheatmap
logFC_matrix <- as.matrix(df[ , "logFC", drop = FALSE])
#logFC_matrix <- as.matrix(df["log2FoldChange"])
#colnames(logFC_matrix) <- "log2FoldChange"

# Create a matrix of significance labels (asterisks)
get_sig_label <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
sig_labels <- matrix(get_sig_label(df$adj.P.Val), ncol = 1)
rownames(sig_labels) <- rownames(df)
colnames(sig_labels) <- "logFC"

library(pheatmap)
# Range of your matrix
min_val <- min(logFC_matrix, na.rm = TRUE)
max_val <- max(logFC_matrix, na.rm = TRUE)

# How many colors total / per side
n_total <- 200
n_neg   <- 100
n_pos   <- 100

# Colors: blue→white for negatives, white→red for positives
cols_neg <- colorRampPalette(c("dodgerblue4", "white"))(n_neg)
cols_pos <- colorRampPalette(c("white", "red"))(n_pos)
my_colors <- c(cols_neg, cols_pos)

# Breaks: min→0 (n_neg+1 steps), 0→max (n_pos+1 steps)
br_neg <- seq(min_val, 0, length.out = n_neg + 1)
br_pos <- seq(0,      max_val, length.out = n_pos + 1)
my_breaks <- c(br_neg, br_pos[-1])  # drop duplicate 0

#ord <- order(as.numeric(sub("MMP", "", rownames(logFC_matrix))))
#logFC_matrix_ordered <- logFC_matrix[ord, , drop = FALSE]

MMP_dict <- list(
  Collagenases = c("MMP1", "MMP8", "MMP13"),
  Gelatinases  = c("MMP2", "MMP9"),
  Stromelysins = c("MMP3", "MMP10", "MMP11"),
  Matrilysins  = c("MMP7", "MMP26"),
  MT_MMPs      = c("MMP14", "MMP15", "MMP16", "MMP24", "MMP17", "MMP25"),
  Elastase     = c("MMP12"),
  Other        = c("MMP19", "MMP20", "MMP21", "MMP23", "MMP23A", "MMP23B", "MMP27", "MMP28")
)

# Flatten to one ordered vector
MMP_order <- unlist(MMP_dict, use.names = FALSE)
keep <- MMP_order[MMP_order %in% rownames(logFC_matrix)]

# reorder logFC_matrix rows
logFC_matrix <- logFC_matrix[keep, , drop = FALSE]


pvec <- df[rownames(logFC_matrix), "adj.P.Val"]  # df_first has unique GeneSymbols

# 2) build the label matrix with matching dimnames
get_sig_label <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*",  "")))
}

sig_labels <- matrix(get_sig_label(pvec),
                     nrow = nrow(logFC_matrix), ncol = 1,
                     dimnames = list(rownames(logFC_matrix), "logFC"))

# optional: replace NAs with blank
sig_labels[is.na(sig_labels)] <- ""

pheatmap(
  logFC_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  breaks = my_breaks,
  display_numbers = sig_labels,   # your "*" matrix
  number_color = "black",         # asterisks in black
  fontsize_number = 12,
  cellwidth = 40,
  cellheight = 14,
  border_color = "black",         # cell borders in black
  main = "Bulk RNA-seq Log2 Fold Change\n(MMP/TIMP Genes)",
  legend = TRUE
)

####################################### Boxplots ###########################
expr_mat <- v_no_batch
combined_metadata$SampleID <- rownames(combined_metadata)
expr_mat <- expr_mat[, combined_metadata$SampleID]
genes_to_plot <- c("ENSG00000102996")  # example gene IDs
genes_to_plot <- intersect(genes_to_plot, rownames(expr_mat))  # keep valid ones


# ----------- 5. Melt expression matrix to long format ------------


meta_sorted<- combined_metadata[match(colnames(expr_mat), combined_metadata$SampleID), ]
colnames(expr_mat) <- meta_sorted$Condition


expr_mat <- as.data.frame(expr_mat)
expr_mat$gene_ensemlid <- rownames(expr_mat)

expr_mat_long <- as.data.frame(expr_mat) %>%
  pivot_longer(cols = starts_with("KC") | starts_with("KR"),
               names_to = "Condition",
               values_to = "Expression") 


expr_mat_long_gene <- subset(expr_mat_long, gene_ensemlid %in% genes_to_plot)
expr_mat_long_gene$Condition <- factor(expr_mat_long_gene$Condition, levels = c("KC", "KR"))


# ----------- 7. Plot boxplots ------------
expr_mat_long_gene$Condition <- factor(expr_mat_long_gene$Condition,
                                       levels = c("KR", "KC"))


ggplot(expr_mat_long_gene, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2.5, alpha = 1) +
  stat_compare_means(
    comparisons = list(c( "KR", "KC")),
    method = "wilcox.test", label = "p.signif"
  ) +
  scale_fill_manual(values = c("#337FC2", "#C03A30")) +
  labs(title = "Expression of Selected Genes", y = "Normalized Expression", x = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"))



