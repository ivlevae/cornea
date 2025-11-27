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
regex <- sprintf("^%s\\d{1,2}([A-OR-Z]?\\d?)?$", 'MMP') 
subset_mmp_2 <- grep(regex, annotation_data$GeneSymbols, value = TRUE)
### extract TIMP genes
regex <- sprintf("^%s\\d{1,2}([A-OR-Z]?\\d?)?$", 'TIMP') 
subset_timp_2 <- grep(regex, annotation_data$GeneSymbols, value = TRUE)


subset_combined <- c(subset_mmp_2) #, subset_timp_2)



file_path <- "/home/bnvlab2/Documents/Kate/Cornea_bulk/GSE204791_meta.csv"
meta <- read_csv(file_path, show_col_types = FALSE)



# ---- Function to process a tissue type ----
tissue_name <- 'corneal stroma'
tissue_name <- 'corneal epithelium'

######### to calculate all MMP DEG ##################

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
            

View(tT_mmp)
tT_mmp_filtered <- tT_mmp %>% filter(P.Value < 0.05) %>% arrange(logFC)



tT_mmp_filtered <- tT_mmp %>% filter(P.Value < 0.05 )

##########   save the results ########################  
  
# write.csv(tT_mmp, "microarrays_stroma_DE.csv", quote = FALSE, row.names = FALSE)
deg_table <- read_csv('microarrays_stroma_DE.csv')  
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


expression_data_mmp15 <- subset(expression_data_long1, GeneSymbols == 'MMP15')

######################## plot boxplot for mmp 15 ##########################

p_mmp15 <- deg_table[deg_table$GeneSymbols == "MMP15", "adj.P.Val"][1]

ggplot(expression_data_mmp15, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.16, size = 3.5, alpha = 1) +
  stat_compare_means(
  comparisons = list(c("KC", "control")),
  method = "wilcox.test", label = "p.signif"
   ) +
  annotate("text",
           x = 1.5, y = 9.8,
           label = paste0("p_value_adj = ", signif(p_mmp15, 2))) +
  scale_fill_manual(values = c("#337FC2", "#C03A30")) +
  labs(title = "Expression of Selected Genes", y = "Normalized Expression", x = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  coord_cartesian(ylim = c(5, 10))

