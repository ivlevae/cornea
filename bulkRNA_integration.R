### GSE151631_meta.csv 
count_data_GSE151631 <- count_data11_ensemlid
count_data_GSE77938  <- count_data 
count_data_GSE77938_rep  <- count_data_rep 
dim(count_data_GSE151631)
dim(count_data_GSE77938)
count_matrix1 <- count_data_GSE151631
count_matrix2 <- count_data_GSE77938
count_matrix3 <- count_data_GSE77938_rep


View(count_data_rep)
# Perform inner join by RowName
combined_count_matrix_2_3 <- cbind(count_matrix2, count_matrix3)

combined_count_matrix <- combined_count_matrix_2_3
View(combined_count_matrix)

metadata2$Batch <- 'Batch_2'
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
  filter(P.Value < 0.05) %>%
  arrange(logFC)
# View the filtered rows
View(mmp_filtered)







#########################################################

count_data11_ensemlid1


count_matrix1 <- cbind(RowName = rownames(count_data11_ensemlid1), count_data11_ensemlid1)
View(count_matrix1)
count_matrix2 <- cbind(RowName = rownames(combined_count_matrix_2_3), combined_count_matrix_2_3)
View(count_matrix2)
# Perform inner join by RowName
combined_count_matrix <- merge(count_matrix1, count_matrix2, by = "RowName", all = FALSE)

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

head(mmp_filtered)
ggplot(mmp_filtered, aes(x = external_gene_name, y = logFC)) + geom_bar(stat = 'identity')
