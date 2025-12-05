library(GEOquery)
library(annotate)
library(org.Hs.eg.db)
library(data.table)
library(reshape2)
library(limma)
library(dplyr)
library(EnhancedVolcano)
library(ggpubr)

## load expressions ##
ex <- read.csv("GSE204791_expression_log2_normalized.csv", row.names = 1)
ex <- as.matrix(ex)                     
# Load annotation file
annotation_data <- fread("/home/bnvlab2/Downloads/GPL21185_noParents.an.txt.gz", fill = TRUE, header = TRUE, skip = 8)
annotation_data1 <- annotation_data[, 1:2]


### extract MMP genes

mmp <- annotation_data[grepl("^MMP\\d{1,2}([A-OR-Z]?\\d?)?$", annotation_data$GeneSymbol, ignore.case = TRUE), ]$GeneSymbol



subset_combined <- c(mmp) #, subset_timp_2)



file_path <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE204791_meta.csv"
meta <- read_csv(file_path, show_col_types = FALSE)



# ---- Function to process a tissue type ----
tissue_name <- 'corneal stroma'
tissue_name <- 'corneal epithelium'

######### to calculate all MMP DEG ##################
ex2 <- merge(results, annotation_data1, by = 'ProbeName', all.x = TRUE)




meta1 <- meta %>% filter(Tissue == tissue_name)
common_samples <- intersect(colnames(ex), meta1$SampleID)
meta1 <- meta1[meta1$SampleID %in% common_samples, ]
ex1 <- ex[, common_samples]

meta1 <- meta1[match(colnames(ex1), meta1$SampleID), ]
colnames(ex1) <- meta1$Condition
colnames(ex1)
# Prepare design and grouping
gs <- factor(meta1$Condition)
gs <- relevel(gs, ref = "control")  
design <- model.matrix(~ gs + 0)
colnames(design) <- levels(gs)


# Fit the linear model
fit <- lmFit(ex1, design)

contrast_name <- paste(levels(gs)[2], levels(gs)[1], sep = "-")
cont.matrix <- makeContrasts(contrasts = contrast_name, levels = design)

# Apply contrast and empirical Bayes moderation
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, adjust = "fdr", number = Inf)
results$ProbeName <- rownames(results)
tT1 <- merge(results, annotation_data1, by = 'ProbeName', all.x = TRUE)
results10 <- subset(tT1, GeneSymbols != '' )
rownames(results10) <- results10$ProbeName 
tT_mmp <- subset(results10, GeneSymbols %in% subset_combined)
            
tT_mmp_epi <- tT_mmp
#tT_mmp_stroma <- tT_mmp

##########   save the results ########################  
  
# write.csv(tT_mmp, "microarrays_stroma_DE.csv", quote = FALSE, row.names = FALSE)
  
#### plot volcano plot ######
filtered_genes <- tT_mmp_filtered$GeneSymbols                
# Volcano Plot
EnhancedVolcano(results10,
                lab = results10$GeneSymbols,
                selectLab = as.character(tT_mmp_filtered$GeneSymbols),
                x = 'logFC',
                y = 'P.Value',
                col = c("grey30",  "forestgreen", "royalblue", "darkorange"),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                #max.overlaps = Inf,
                drawConnectors = TRUE,
                colAlpha = 0.7,
                title = paste("Volcano Plot: Highlighting MMPs in", tissue_name)) +
  theme_classic() +
  xlim(-5, 4) +   # xlim(-10, 8) + ##stroma
  ylim(0, 10)    #ylim(0, 13)  ##stroma
  




######################## boxplot expression #############################


tissue_name <- 'corneal stroma'
tissue_name <- 'corneal epithelium'


expression_data_long <- melt(ex, variable.name = "Var2", value.name = "Expression")
colnames(expression_data_long)[1] <- "ProbeName"


meta1 <- meta %>% filter(Tissue == tissue_name)

# Merge with metadata
expression_data_long <- merge(expression_data_long, meta1, by.x = "Var2", by.y = "SampleID")
expression_data_long1 <- merge(expression_data_long, annotation_data, by = "ProbeName")

expression_data_mmps <- subset(expression_data_long1, GeneSymbols %in% mmp)
######### plot heatmap 
head(tT_mmp)
df <- tT_mmp_epi[c('GeneSymbols', 'logFC', 'adj.P.Val')]
# Set gene names as rownames
df_first <- df[!duplicated(df$GeneSymbols), , drop = FALSE]

rownames(df_first) <- df_first$GeneSymbols

# Create the logFC matrix for pheatmap
logFC_matrix <- as.matrix(df_first[ , "logFC", drop = FALSE])
#logFC_matrix <- as.matrix(df["log2FoldChange"])
#colnames(logFC_matrix) <- "log2FoldChange"
# Flatten to one ordered vector
MMP_order <- unlist(MMP_dict, use.names = FALSE)
keep <- MMP_order[MMP_order %in% rownames(logFC_matrix)]

# reorder logFC_matrix rows
logFC_matrix <- logFC_matrix[keep, , drop = FALSE]

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
min_val <- -9.578429
print(min_val) ### -9.578429
max_val <- max(logFC_matrix, na.rm = TRUE)
print(max_val) ### 1.704696
max_val <- 1.704696
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



pvec <- df_first[rownames(logFC_matrix), "adj.P.Val"]  # df_first has unique GeneSymbols

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

# 3) plot
pheatmap(
  logFC_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  breaks = my_breaks,
  display_numbers = sig_labels,   # now same dim/order as mat
  number_color = "black",
  fontsize_number = 12,
  cellwidth = 40,
  cellheight = 14,
  border_color = "black",
  main = "Bulk RNA-seq Log2 Fold Change\n(MMP/TIMP Genes)",
  legend = TRUE
)

##############################################################

expression_data_mmp15 <- subset(expression_data_long1, GeneSymbols == 'MMP15')

######################## plot boxplot for mmp 15 ##########################


ggplot(expression_data_mmp15, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2.5, alpha = 1)  +
  stat_compare_means(
    comparisons = list(c("KC", "control")),
    method = "wilcox.test", label = "p.signif"
  ) +
  scale_fill_manual(values = c("#337FC2", "#C03A30")) +
  labs(title = "Expression of Selected Genes", y = "Normalized Expression", x = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  coord_cartesian(ylim = c(5, 10))

