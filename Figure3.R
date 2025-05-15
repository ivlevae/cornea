###########  Endothelium 
FeaturePlot(cornea, 
            features= c('POU6F2',  'COL4A3', 'COL8A1', 'COL4A4' , 'CA12', 'SLC4A11', 'SLC4A4'), 
            raster = T,
            order = T, reduction = 'scVI') &NoAxes() & NoLegend()



###########  Keratocytes
FeaturePlot(cornea, 
            features= c( 'KERA', 'NNMT', 'DCN', 'TIMP2'), raster = T, reduction = 'scVI',
            order = T) &NoAxes() & NoLegend()


###################### Epitheliumn 
FeaturePlot(cornea, 
            features= c(
              'PAX6',    'KRT5', 'KRT24', 'MYH14', 'DSP', 'KRT15' ), 
            raster = F, reduction = 'scVI'#,
           # order = T
            ) &NoAxes() & NoLegend()


################## All together ###################

FeaturePlot(cornea, 
            features= c( 'POU6F2',  'COL4A4', 'COL8A1',  'MMP17',  
                         'DCN', 'COL6A2',  'COL12A1', 'MMP3', ## 'COL6A3', 'COL1A2', 
                         'PAX6',    'KRT3', 'KRT5', #'COL7A1', 
                         'ADAMTS14'), 
            raster = F,
            ncol = 4, reduction = 'scVI', order = T
            ) &NoAxes() & NoLegend()




###################################### Pheatmaps ########################################
desired_order <- c("Corneal Superficial",
                   "Corneal Wing",
                   "Corneal Basal",
                   "TAC",
                   #'Limbal Suprabasal',
                   "Limbal Basal",
                   "Keratocytes",
                   "Myofibroblasts",
                   'Limbal fibroblasts',
                   "Endothelium") 


extract_genes <- function(gene_list, pattern) {
  regex <- sprintf("^%s\\d{1,2}([A-OR-Z]?\\d?)?$", pattern) 
  pattern_genes <- grep(regex, gene_list, value = TRUE)
  
  gene_expression <- FetchData(cornea, vars = c(pattern_genes))
  mean_expression <- colMeans(gene_expression, na.rm = TRUE)
  zero_expr_genes <- names(mean_expression[mean_expression == 0])
  pattern_genes_filtered <- pattern_genes[!pattern_genes %in% zero_expr_genes]
  
  return(pattern_genes_filtered)
}

all_MMP <- extract_genes(rownames(cornea), 'MMP')
all_MMP

all_adam <- extract_genes(rownames(cornea), 'ADAMTS')
all_adam

all_timp <- extract_genes(rownames(cornea), 'TIMP')
all_timp

all_cols <- extract_genes(rownames(cornea), 'COL')
all_cols

all_krt<- extract_genes(rownames(cornea), 'KRT')
all_krt


ecm_genes <- c(all_MMP, all_adam, all_timp)


##############################    Heatmap calculation and plotting #################################
cornea$leiden_annot_V3 <- factor(cornea$leiden_annot_V3, levels = desired_order)
cornea_healthy <- subset(cornea, condition_detailed == 'Healthy')



############ change gene_list for ecm_genes, all_krt, all_cols
gene_list <- all_krt  
DefaultAssay(cornea) <- 'RNA'


avg_table <- AverageExpression(cornea_healthy,  assays = 'RNA', features =  gene_list, group.by = c('leiden_annot_V3'))
avg_table <- as.data.frame(avg_table)
avg_table_transposed <- as.data.frame(t(avg_table))

scaled_gene_expression <- avg_table_transposed %>%
  mutate(across(all_of(gene_list), ~ (.-min(.)) / (max(.) - min(.))))


rownames(scaled_gene_expression) <- desired_order
scaled_gene_expression <- scaled_gene_expression[rev(desired_order), ]


pheatmap_result <- pheatmap(
  scaled_gene_expression,
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = T,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0,
  angle_col = 90
)


################# order for diagonal 

dim(scaled_gene_expression)
length(ordered_cols2)
ordered_ecm <- c(   "MMP24", "ADAMTS16", "ADAMTS13", "ADAMTS19",
                     "MMP17", "ADAMTS6", "TIMP3", "MMP11", "MMP14", "ADAMTS10", "ADAMTS2",
                     "ADAMTS8", "MMP19", "ADAMTS5", "ADAMTS3", "MMP27", "MMP23B", "TIMP4",
                     "MMP16", "ADAMTS15", "ADAMTS7", "MMP2", "TIMP2" , "MMP21", "MMP25",   "ADAMTS1", "ADAMTS12", "ADAMTS4",
                     "ADAMTS9",  "MMP13", "TIMP1", "MMP8", "MMP10", "MMP1", "MMP3", "MMP12",  "MMP9",
                     "MMP7", "MMP20", "MMP15",   "MMP28", "MMP26",  "ADAMTS20",   "ADAMTS18" ,"ADAMTS17", "ADAMTS14", "MMP23A" 
)

pheatmap(
  scaled_gene_expression[, ordered_ecm],
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0,
  angle_col = 90
)


dim(scaled_gene_expression)
length(ordered_cols2)

ordered_cols2 <- c(  "COL8A1", "COL4A4", "COL4A3",   "COL4A6", "COL4A5", "COL27A1",  "COL26A1",
                      "COL5A2", "COL8A2","COL19A1",  "COL20A1",  "COL11A2", "COL2A1",'COL6A5',
                     "COL24A1",  "COL15A1", "COL1A1", "COL6A6", "COL10A1",
                     "COL9A3", "COL5A1", "COL11A1", "COL14A1", "COL16A1", "COL1A2", "COL3A1",
                    "COL25A1",  "COL23A1","COL6A3",    "COL6A2",  "COL13A1",
                     "COL6A1", "COL18A1",  "COL9A1",  "COL4A1", "COL4A2",   "COL9A2", "COL12A1","COL17A1",
                     "COL7A1", "COL22A1",  "COL5A3", "COL21A1",  "COL28A1"
                     )
pheatmap(
  scaled_gene_expression[, ordered_cols2],
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0,
  angle_col = 90
)

dim(scaled_gene_expression)
length(ordered_krts)


ordered_krts <- c( "KRT25",  "KRT79", "KRT1", "KRT36", "KRT222",
                   "KRT26", "KRT28","KRT85", "KRT10",  "KRT81",  "KRT9", "KRT15", 
                   "KRT6C", "KRT32", "KRT75",  "KRT19", "KRT23", "KRT33B", "KRT35", 
                   "KRT6B", "KRT6A", "KRT16", "KRT17", "KRT80", "KRT34",   "KRT31", 
                   "KRT18","KRT5",  "KRT7", "KRT14", "KRT13",  "KRT8", "KRT86", "KRT39",  "KRT37",
                    "KRT71",  "KRT74", "KRT77" ,"KRT76", 
                    "KRT72",     "KRT2",   "KRT40", "KRT33A","KRT20", "KRT12", "KRT27",  "KRT78", "KRT73","KRT3","KRT4",
                     "KRT24",     "KRT84")


pheatmap(
  scaled_gene_expression[, ordered_krts],
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0,
  angle_col = 90
)


