
library("DESeq2")
library(tidyverse)

###############              GSE77938       17 controls + 17 KC

file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_replication_gene_counts.txt.gz"
count_data_rep <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

metadata_rep <- data.frame(
  row.names = colnames(count_data_rep),  
  Condition = c( rep("Keratoconus", 17), rep("Control", 17)) 
)

metadata_rep$Condition <- as.factor(metadata_rep$Condition)
metadata_rep$Condition <- relevel(metadata_rep$Condition, ref = "Control")


dds <- DESeqDataSetFromMatrix(  countData = count_data_rep, 
                                colData = metadata_rep,      
                                design = ~ Condition)

dds <- DESeq(dds)

results <- results(dds, contrast = c("Condition", "Keratoconus", "Control"))
results_df <- as.data.frame(results)
View(results_df)

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
  filter(pvalue < 0.05) %>%
  arrange(log2FoldChange)
# View the filtered rows
View(mmp_filtered)


####################################   GSE77938       8 controls + 8 KC


file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_discovery_gene_counts.txt.gz"

count_data <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
View(count_data)


metadata2 <- data.frame(
  row.names = colnames(count_data),  
  Condition = c(rep("KR", 8), rep("KC", 8)) 
)
metadata2$Condition <- as.factor(metadata2$Condition)
metadata2$Condition <- relevel(metadata2$Condition, ref = "KR")


dds <- DESeqDataSetFromMatrix(  countData = count_data, 
                                colData = metadata2,      
                                design = ~ Condition)

dds <- DESeq(dds)

results <- results(dds, contrast = c("Condition", "KC", "KR"))
results_df <- as.data.frame(results)
View(results_df)
plotMA(results, ylim = c(-2, 2))

# Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot')




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
mmp_all <- data_with_hugo[grepl("^MMP\\d{1,2}$", data_with_hugo$external_gene_name, ignore.case = TRUE), ]

mmp_filtered <- mmp_all %>%
                filter(pvalue < 0.05) %>%
                arrange(log2FoldChange)
# View the filtered rows
View(mmp_filtered)



##############################################  GSE151631 ##################################

file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE151631_Raw.counts.csv.gz"
count_data1 <- read.csv(file_name, header = TRUE)
View(count_data1)

duplicates <- count_data1 |>
  group_by(Gene) |>
  filter(n() > 1) |>
  ungroup()
View(duplicates)

count_data11 <- count_data1 |>
  group_by(Gene) |>
  summarise(
    across(where(is.numeric), ~ as.integer(floor(mean(.x, na.rm = TRUE)))), # Round down and format as integer
    across(where(~ !is.numeric(.)), ~ first(.x)),                           # Preserve the first value for non-numeric columns
    .groups = "drop"
  )
View(count_data11)
dim(count_data1)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
gene_mapping1 <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  # Attributes to retrieve
  mart = ensembl
)
View(gene_mapping1)

gene_mapping1_fill <- gene_mapping1 %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "", ensembl_gene_id, external_gene_name))

gene_mapping1_clean <- gene_mapping1_fill %>%
  group_by(external_gene_name) %>%
  filter(row_number() == 1) %>%
  ungroup()


View(gene_mapping1_fill)


duplicates_in_mapping <- gene_mapping1_clean[duplicated(gene_mapping1_clean$external_gene_name), ]

# View the duplicated entries (if any)
View(duplicates_in_mapping)


count_data11_ensemlid <- count_data11

count_data11_ensemlid <- merge(count_data11_ensemlid, gene_mapping1_clean, by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)
count_data11_ensemlid1 <- subset(count_data11_ensemlid, !(is.na(ensembl_gene_id) | ensembl_gene_id == ""))

rownames(count_data11_ensemlid1) <- count_data11_ensemlid1$ensembl_gene_id
count_data11_ensemlid1$ensembl_gene_id <- NULL
count_data11_ensemlid1$Gene <- NULL


View(count_data11_ensemlid)
View(count_data11_ensemlid1)
dim(count_data11_ensemlid1)


raw_meta <- read.csv("/home/bnvlab2/Documents/Kate/Alzheimer/GSE151631_meta.csv", sep='')
View(raw_meta)
# raw_meta$Sample_ID <- factor(raw_meta$Sample_ID, levels = colnames(count_data1)[-1])
# 
# # Optionally, reorder the rows of df based on the new ordering of sample_id
# raw_meta <- raw_meta[order(raw_meta$Sample_ID), ]
# View(raw_meta)
# 
# cases <- as.character(subset(raw_meta, Group == 'case')$Sample_ID)
# controls <-  as.character(subset(raw_meta, Group == 'control')$Sample_ID)


count_data11 <- as.data.frame(count_data11)
rownames(count_data11) <- count_data11$Gene
count_data11$Gene <- NULL
View(count_data11)

metadata_GSE151631 <- data.frame(
  row.names = raw_meta$Sample_ID,  
  Condition = raw_meta$Group  
)
metadata_GSE151631$Condition <- as.factor(metadata_GSE151631$Condition)
metadata_GSE151631$Condition <- relevel(metadata_GSE151631$Condition, ref = "control")

View(metadata)
View(count_data112)
count_data112 <- count_data11[, as.character(raw_meta$Sample_ID)]


dds <- DESeqDataSetFromMatrix(
  countData = count_data112,  
  colData = metadata,      
  design = ~ Condition     
)

dds <- DESeq(dds)

results <- results(dds, contrast = c("Condition", "case", "control"))


results_df <- as.data.frame(results)
View(results_df1)
results_df$gene_name <- rownames(results_df)

mmp_all <- results_df[grepl("^MMP\\d{1,2}$", results_df$gene_name, ignore.case = TRUE), ]

mmp_filtered <- mmp_all %>%
  filter(pvalue < 0.05) %>%
  arrange(log2FoldChange)
# View the filtered rows
View(mmp_filtered)



library(EnhancedVolcano)
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot')




###################### TPM ##############################################


# 
file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_discovery_gene_tpm.txt.gz"
data <- fread(file_name)
colnames(data)
mmp_all <- data[grepl("^MMP\\d{1,2}$", data$`Gene symbol`, ignore.case = TRUE), ]$`Gene symbol`
View(mmp_all)

data1 <- data[data$`Gene symbol` %in% mmp_all]
View(data1)
data1 <- data1 %>%
  mutate(GeneSymbols_num = factor(`Gene symbol`,
                                  levels = unique(`Gene symbol`[order(as.numeric(gsub("MMP", "", `Gene symbol`)))])))

View(data1)
colnames(data1)
data1$GeneID <- NULL
data1$`Gene biotype` <- NULL

colnames(data1)[2:17] <- c( rep("KR", 8), rep("KC", 8))
#log_ex <- log2(ex + 1) 


data_long <- data1 %>%
  pivot_longer(cols =  starts_with("KR") | starts_with("KC"),  # select columns that start with 'KC' or 'KR'
               names_to = "Condition",  # New column 'Condition' to store KC or KR
               values_to = "Expression")  # New column 'Expression' for values
mean_diff

data_long$Expression_log <- log2(data_long$Expression + 1)
View(data_long)


# Create a boxplot for each Gene symbol between KC and KR samples
ggplot(data_long, aes(x = Condition, y = Expression_log, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~`Gene symbol`, scales = "free_y") +  # Facet by each Gene symbol
  labs(title = "Expression_log of Genes between KC and KR",
       x = "Condition",
       y = "Expression_log Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



###################### 


mean_expression <- data_long %>%
  group_by(`Gene symbol`, Condition) %>%
  summarise(mean_expression = mean(Expression_log, na.rm = TRUE), .groups = 'drop')

# Step 4: Calculate the difference in means between KR and KC for each gene
mean_diff <- mean_expression %>%
  spread(key = Condition, value = mean_expression) %>%
  mutate(mean_diff = KR - KC) %>%
  arrange(desc(mean_diff))  # Order by the difference in means


order_MMP <- mean_diff$`Gene symbol`  # Store ordered gene symbols in a vector

# Step 6: Reorder `Gene symbol` in the `data_long` dataframe using `order_MMP`
data_long <- data_long %>%
  mutate(`Gene symbol` = factor(`Gene symbol`, levels = order_MMP))  # Reorder based on the extracted order_MMP

# Step 7: Create a boxplot for each Gene symbol, colored by Condition, keeping the custom order
ggplot(data_long, aes(x = Condition, y = Expression_log, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~`Gene symbol`, scales = "free_y") +  # Facet by each Gene symbol
  labs(title = "Expression_log of Genes between KC and KR",
       x = "Condition",
       y = "Expression_log Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# # 
file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE77938_replication_gene_tpm.txt.gz"

data <- fread(file_name)
colnames(data)
mmp_all <- data[grepl("^MMP\\d{1,2}$", data$`Gene symbol`, ignore.case = TRUE), ]$`Gene symbol`
View(mmp_all)

data1 <- data[data$`Gene symbol` %in% mmp_all]
#View(data1)
data1 <- data1 %>%
  mutate(GeneSymbols_num = factor(`Gene symbol`,
                                  levels = unique(`Gene symbol`[order(as.numeric(gsub("MMP", "", `Gene symbol`)))])))


length(colnames(data1))

data1$GeneID <- NULL
data1$`Gene biotype` <- NULL
View(data1)
colnames(data1)
colnames(data1)[2:35] <- c( rep("KC", 17), rep("KR", 17))

data_long <- data1 %>%
  pivot_longer(cols = starts_with("KC") | starts_with("KR"),  # select columns that start with 'KC' or 'KR'
               names_to = "Condition",  # New column 'Condition' to store KC or KR
               values_to = "Expression")  # New column 'Expression' for values
data_long$Expression_log <- log2(data_long$Expression + 1)
# Create a boxplot for each Gene symbol between KC and KR samples
ggplot(data_long, aes(x = Condition, y = Expression_log, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~`Gene symbol`, scales = "free_y") +  # Facet by each Gene symbol
  labs(title = "Expression of Genes between KC and KR",
       x = "Condition",
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





mean_expression <- data_long %>%
  group_by(`Gene symbol`, Condition) %>%
  summarise(mean_expression = mean(Expression_log, na.rm = TRUE), .groups = 'drop')

# Step 4: Calculate the difference in means between KR and KC for each gene
mean_diff <- mean_expression %>%
  spread(key = Condition, value = mean_expression) %>%
  mutate(mean_diff = KR - KC) %>%
  arrange(desc(mean_diff))  # Order by the difference in means


order_MMP <- mean_diff$`Gene symbol`  # Store ordered gene symbols in a vector

# Step 6: Reorder `Gene symbol` in the `data_long` dataframe using `order_MMP`
data_long <- data_long %>%
  mutate(`Gene symbol` = factor(`Gene symbol`, levels = order_MMP))  # Reorder based on the extracted order_MMP

# Step 7: Create a boxplot for each Gene symbol, colored by Condition, keeping the custom order
ggplot(data_long, aes(x = Condition, y = Expression_log, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~`Gene symbol`, scales = "free_y") +  # Facet by each Gene symbol
  labs(title = "Expression_log of Genes between KC and KR",
       x = "Condition",
       y = "Expression_log Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#################### GSE241715   mormilized data only 
# 
# file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE241715_normal_data.csv.gz"
# count_data3 <- read.csv(file_name, header = TRUE)
# count_data31 <- count_data3[, c(1,  27, 28, 31, 32)] #6:16,
# View(count_data31)
# colnames(count_data31)
# mmp_all <- count_data31[grepl("^MMP\\d{1,2}$", count_data31$Symbol, ignore.case = TRUE), ]
# View(mmp_all)



######################################## GSE112155 
# 
# file_name <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE112155_counts_table.txt.gz"
# count_data <- read.table(file_name, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# count_data1 <- count_data[, 1:20]
# 
# metadata_myop <- data.frame(
#   row.names = colnames(count_data1),  
#   Condition = c(  rep("Control", 10), rep("Keratoconus", 10)) 
# )
# 
# metadata_myop$Condition <- as.factor(metadata_myop$Condition)
# metadata_myop$Condition <- relevel(metadata_myop$Condition, ref = "Control")
# 
# 
# dds <- DESeqDataSetFromMatrix(  countData = count_data1, 
#                                 colData = metadata_myop,      
#                                 design = ~ Condition)
# 
# dds <- DESeq(dds)
# 
# results <- results(dds, contrast = c("Condition", "Keratoconus", "Control"))
# results_df <- as.data.frame(results)
# View(results_df)
# 
# results_df$EnsemblGeneID <- rownames(results_df)
# ensembl_ids <- results_df$EnsemblGeneID
# 
# 
# # Connect to Ensembl biomart
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# gene_mapping <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name"),  # Attributes we want to retrieve
#   filters = "ensembl_gene_id",  # Filter based on Ensembl gene IDs
#   values = ensembl_ids,  # Ensembl gene IDs to look up
#   mart = ensembl
# )
# 
# data_with_hugo <- merge(results_df, gene_mapping, by.x = "EnsemblGeneID", by.y = "ensembl_gene_id", all.x = TRUE)
# View(data_with_hugo)
# 
# mmp_all <- data_with_hugo[grepl("^MMP\\d{1,2}$", data_with_hugo$external_gene_name, ignore.case = TRUE), ]
# 
# mmp_filtered <- mmp_all %>%
#  filter(pvalue < 0.05) %>%
#   arrange(log2FoldChange)
# 
# # View the filtered rows
# View(mmp_filtered)

