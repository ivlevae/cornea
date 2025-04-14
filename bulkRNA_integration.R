file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_replication_gene_counts.txt.gz"
count_data_GSE77938_rep <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)


file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_discovery_gene_counts.txt.gz"
count_data_GSE77938 <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)


file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE151631_Raw.counts.csv.gz"
count_data_GSE151631 <- read.csv(file_name, header = TRUE)



############################## Combine count_data_GSE77938 study ############################

combined_count_matrix_2_3 <- cbind(count_data_GSE77938, count_data_GSE77938_rep)

combined_count_matrix <- combined_count_matrix_2_3
View(combined_count_matrix)


metadata2 <- data.frame(
  row.names = colnames(count_data),  
  Condition = c(rep("KR", 8), rep("KC", 8)) 
)
metadata2$Condition <- as.factor(metadata2$Condition)
metadata2$Condition <- relevel(metadata2$Condition, ref = "KR")
metadata2$Batch <- 'Batch_2'


metadata_rep <- data.frame(
  row.names = colnames(count_data_rep),  
  Condition = c( rep("KC", 17), rep("KR", 17)) 
)

metadata_rep$Condition <- as.factor(metadata_rep$Condition)
metadata_rep$Condition <- relevel(metadata_rep$Condition, ref = "KR")
metadata3 <- metadata_rep
metadata3$Batch = 'Batch_1'

# Combine metadata into one
combined_metadata <- rbind(metadata2, metadata3)

View(combined_metadata)

# Ensure metadata is in the correct format
combined_metadata$Condition <- factor(combined_metadata$Condition)
combined_metadata$Batch <- factor(combined_metadata$Batch)
combined_metadata$Condition <- relevel(combined_metadata$Condition, ref = "KR")

dge <- DGEList(counts = combined_count_matrix)
dge <- calcNormFactors(dge)  # Normalize

dim(combined_count_matrix)
dim(dge)
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


# Connect to Ensembl biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  # Attributes we want to retrieve
  filters = "ensembl_gene_id",  # Filter based on Ensembl gene IDs
  values = ensembl_ids,  # Ensembl gene IDs to look up
  mart = ensembl
)

data_with_hugo <- merge(results_df, gene_mapping, by.x = "EnsemblGeneID", by.y = "ensembl_gene_id", all.x = TRUE)
#View(data_with_hugo)

mmp_all <- data_with_hugo[grepl("^MMP\\d{1,2}$", data_with_hugo$external_gene_name, ignore.case = TRUE), ]

mmp_filtered <- mmp_all %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(logFC)
# View the filtered rows
View(mmp_filtered)


write.csv(mmp_filtered, 'mmp_filtered_GSE77938_whole_study.csv')




#########################################################

count_data11_ensemlid1


count_data_GSE151631 <- cbind(RowName = rownames(count_data11_ensemlid1), count_data11_ensemlid1)
View(count_data_GSE151631)
count_data_GSE77938 <- cbind(RowName = rownames(combined_count_matrix_2_3), combined_count_matrix_2_3)
View(count_data_GSE77938)
# Perform inner join by RowName
combined_count_matrix <- merge(count_data_GSE151631, count_data_GSE77938, by = "RowName", all = FALSE)

# Set row names back to the merged data frame (using the 'RowName' column)
rownames(combined_count_matrix) <- combined_count_matrix$RowName
combined_count_matrix <- combined_count_matrix[, -1]  



# Create the metadata (this should include information about the condition and batch)
metadata1 <-  metadata_GSE151631
metadata1$Batch <- 'Batch_3'



# Combine metadata into one
combined_metadata <- rbind(metadata1, combined_metadata)

View(combined_metadata)

combined_metadata <- combined_metadata %>%
  mutate(Condition = recode(Condition, "KR" = "control", "KC" = "case"))
# Ensure metadata is in the correct format
combined_metadata$Condition <- factor(combined_metadata$Condition)
combined_metadata$Batch <- factor(combined_metadata$Batch)
combined_metadata$Condition <- relevel(combined_metadata$Condition, ref = "control")
# Create a DESeqDataSet or DGEList if necessary (optional step based on your workflow)
# For example, using DESeq2:
# dds <- DESeqDataSetFromMatrix(countData = combined_count_matrix, colData = combined_metadata, design = ~ Condition + Batch)

# Alternatively, if using edgeR:
# dge <- DGEList(counts = combined_count_matrix)
# dge$samples$group <- combined_metadata$Condition
dim(combined_count_matrix)
dim(combined_metadata)
# Normalize the count data using `voom` (limma approach)
library(edgeR) 
library(limma)
dge <- DGEList(counts = combined_count_matrix)
dge <- calcNormFactors(dge)  # Normalize

dim(combined_metadata)
dim(dge)
# Use voom to estimate the mean-variance relationship and log2 transformation
v <- voom(dge, design = model.matrix(~ Condition + Batch, data = combined_metadata))

# Apply removeBatchEffect to remove batch effects
v_no_batch <- removeBatchEffect(v$E, batch = combined_metadata$Batch, design = model.matrix(~ Condition, data = combined_metadata))

# Differential expression analysis using limma
fit <- lmFit(v_no_batch, design = model.matrix(~ Condition, data = combined_metadata))
fit <- eBayes(fit)

colnames(fit$coefficients)
# Get the results
results <- topTable(fit, coef = "Conditioncase", number = Inf)

# View results
View(results)

results_df <- as.data.frame(results)
results_df$EnsemblGeneID <- rownames(results_df)
ensembl_ids <- results_df$EnsemblGeneID


# Connect to Ensembl biomart

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  # Attributes we want to retrieve
  filters = "ensembl_gene_id",  # Filter based on Ensembl gene IDs
  values = ensembl_ids,  # Ensembl gene IDs to look up
  mart = ensembl
)

data_with_hugo <- merge(results_df, gene_mapping, by.x = "EnsemblGeneID", by.y = "ensembl_gene_id", all.x = TRUE)
View(data_with_hugo)

mmp_all <- data_with_hugo[grepl("^MMP\\d{1,2}$", data_with_hugo$external_gene_name, ignore.case = TRUE), ]

mmp_filtered <- mmp_all %>%
  filter(P.Value < 0.05) %>%
  arrange(logFC)
# View the filtered rows
View(mmp_filtered)

head(data_with_hugo)
mmp_filtered$external_gene_name <- factor(mmp_filtered$external_gene_name,
                                          levels = mmp_filtered$external_gene_name[order(mmp_filtered$logFC)])

# Plot
ggplot(mmp_filtered, aes(x = external_gene_name, y = logFC)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(x = "Gene", y = "logFC", title = "MMP Genes Sorted by logFC")

View(results)
View(data_with_hugo)

abels <- mmp_filtered$external_gene_name
rownames(data_with_hugo) <- data_with_hugo$external_gene_name


# Step 3: Plot with EnhancedVolcano
EnhancedVolcano(data_with_hugo,
                lab = as.character(gene_labels),                               # all gene names for labels
                selectLab = as.character(mmp_filtered$external_gene_name),     # show labels only for MMP genes
                x = 'logFC',
                y = 'P.Value',
                col = c("grey30", "forestgreen", "royalblue", "darkorange"),
              #  colCustom = gene_colors,
              pCutoff = 0.05,
              FCcutoff = 0.5,         # show everything
               max.overlaps = Inf,
              drawConnectors = T,
              colAlpha = 0.7,
                title = "Volcano Plot: Highlighting MMP Genes") +
  theme_classic() +
  xlim(-7, 4) +
  ylim(0, 18)  

View(data_with_hugo)

# pCutoff = 1e-05,
# pCutoffCol = y,
# FCcutoff = 1,


View(data_with_hugo)
expr_mat <- v_no_batch
combined_metadata$SampleID <- rownames(combined_metadata)
expr_mat <- expr_mat[, combined_metadata$SampleID]

# ----------- 4. Select genes to plot ------------
genes_to_plot <- c("ENSG00000102996")  # example gene IDs
genes_to_plot <- intersect(genes_to_plot, rownames(expr_mat))  # keep valid ones
genes_to_plot
# ----------- 5. Melt expression matrix to long format ------------





meta111 <- combined_metadata[match(colnames(expr_mat), combined_metadata$SampleID), ]
colnames(expr_mat) <- meta111$Condition

View(expr_mat1)
expr_mat1 <- as.data.frame(expr_mat)

expr_mat1$gene_ensemlid <- rownames(expr_mat1)

data_long22 <- as.data.frame(expr_mat1) %>%
  pivot_longer(cols = starts_with("case") | starts_with("control"),
               names_to = "Condition",
               values_to = "Expression") #%>%


View(data_long22)

data_long22_gene <- subset(data_long22, gene_ensemlid %in% genes_to_plot)
View(data_long22_gene)
data_long22_gene$Condition <- factor(data_long22_gene$Condition, levels = c("control", "case"))
# ----------- 7. Plot boxplots ------------


ggplot(data_long22_gene, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2.5, alpha = 1) +
  # facet_wrap(~ gene_ensemlid, scales = "free_y") +
  scale_fill_manual(values = c("#337FC2", "#C03A30")) +
  labs(title = "Expression of Selected Genes", y = ylab, x = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

