###########  Endothelium 
FeaturePlot(prefinalObj, 
            features= c('POU6F2',  'COL4A3', 'COL8A1', 'COL4A4'), 
            raster = F,
            order = T) &NoAxes() & NoLegend()



###########  Keratocytes
FeaturePlot(prefinalObj, 
            features= c( 'KERA', 'COL6A3', 'COL1A2', 'TIMP2'), raster = F,
            order = T) &NoAxes() & NoLegend()


###################### Epitheliumn 
FeaturePlot(prefinalObj, 
            features= c(
              'ELF3', 'TACSTD2'), 
            raster = F,
            order = T) &NoAxes() & NoLegend()


################## All together ###################

FeaturePlot(prefinalObj, 
            features= c( 'POU6F2',  'COL4A4', 'COL8A1',  'MMP17',
                         'KERA', 'COL6A3', 'COL1A2', 'TIMP2',         
                         'TACSTD2',    'KRT3', 'KRT24', #'COL7A1', 
                         'ADAMTS14'), 
            raster = F,
            ncol = 4,
            order = T) &NoAxes() & NoLegend()




###################################### Pheatmaps ########################################



markers_subset_krt <- c("KRT1", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A", "KRT6B",  "KRT7", "KRT8", "KRT9", "KRT10", "KRT12", "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", "KRT18", "KRT19", "KRT20", "KRT23", "KRT24", "KRT27",  "KRT77", "KRT78", "KRT79", "KRT80")
markers_subset_col <- c("COL1A1", "COL1A2",  "COL3A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "COL5A1", "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6", "COL7A1", "COL8A1", "COL8A2", "COL9A1",  "COL9A3", "COL10A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", "COL17A1", "COL18A1",  "COL21A1", "COL22A1",  "COL27A1", "COL28A1")

mmp_genes <- c( 'MMP1', 'MMP2', 'MMP3', 'MMP7', 'MMP8', 'MMP9', 'MMP10', 
                'MMP11',  'MMP12',  'MMP13', 'MMP14', 'MMP15',  'MMP16', 
                'MMP17',  'MMP19',   'MMP20',  'MMP21',  
                'MMP24', 'MMP25', 'MMP26', 'MMP27', 
                'MMP28' ,
                "TIMP1", "TIMP2", "TIMP3", 
                'ADAMTS1', 'ADAMTS2', 'ADAMTS3', 'ADAMTS4', 'ADAMTS5',  'ADAMTS6', 'ADAMTS7', 'ADAMTS8', 'ADAMTS9', 
                'ADAMTS10', 'ADAMTS12', 'ADAMTS13', 'ADAMTS14','ADAMTS17')

############ change gene_list markers_subset_col 
gene_list <- markers_subset_col
DefaultAssay(prefinalObj) <- 'RNA'
gene_expression <- FetchData(prefinalObj,
                             vars = c(gene_list, "detailed_annot", 'condition_detailed'))


avg_table <- AverageExpression(prefinalObj,  assays = 'RNA', features =  gene_list, group.by = c('detailed_annot'))
avg_table <- as.data.frame(avg_table)



# 
avg_table_transposed <- t(avg_table)

avg_table_transposed <-  as.data.frame(avg_table_transposed)
scaled_gene_expression <- avg_table_transposed %>%
  mutate(across(all_of(gene_list), ~ (.-min(.)) / (max(.) - min(.))))


desired_order <- c("Corneal Superficial",
                   "Corneal Wing",
                   "Corneal Basal",
                   "TAC",
                   'Limbal Suprabasal',
                   "Limbal Basal",
                   "Keratocytes",
                   "Myofibroblasts",
                   'Limbal fibroblasts',
                   "Endothelium") 

rownames(scaled_gene_expression) <- desired_order
scaled_gene_expression <- scaled_gene_expression[match(rev(desired_order), rownames(scaled_gene_expression)), ]


pheatmap_result <- pheatmap(
  scaled_gene_expression,
  #ordered_avg_table_df_selected ,
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0
)


ordered_cols_cols <- c("COL10A1", "COL4A5",
                       "COL6A5",  "COL5A2" , "COL8A2" , "COL27A1", "COL4A6" ,"COL8A1" , "COL4A3" , "COL4A4" ,
                       "COL15A1",
                       "COL6A1" , "COL6A2" , "COL16A1", "COL9A1",  "COL1A1",  "COL3A1",  "COL1A2" , "COL6A6" , "COL5A1" , "COL14A1",
                       "COL9A3" , "COL11A2" ,"COL6A3" ,  "COL18A1", "COL4A1" , "COL4A2","COL5A3" ,"COL12A1", "COL13A1",
                       "COL7A1" , "COL17A1","COL22A1",
                       "COL21A1", "COL28A1" )

pheatmap_result <- pheatmap(
  scaled_gene_expression[, ordered_cols_cols],
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0
)


ordered_cols_krt <-c(  "KRT1" ,"KRT20", "KRT9" , "KRT14",  "KRT15" , "KRT8" ,
                       "KRT12", "KRT5",  "KRT18", "KRT10", "KRT7" , "KRT19", "KRT77","KRT16", "KRT3" , 
                       "KRT17", "KRT6A", 
                       "KRT13", "KRT6B" ,"KRT2"  ,"KRT80", "KRT23", "KRT27", "KRT78", "KRT4" , "KRT24", "KRT79")

pheatmap_result <- pheatmap(
  scaled_gene_expression[, ordered_cols_krt],
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0
)



ordered_mmps <- c("MMP17" , "ADAMTS6" , "ADAMTS13",   "MMP24"  ,    "MMP13" ,  "ADAMTS3" , 
                   "TIMP2" ,  "MMP19"  , 
                   "ADAMTS2",  "MMP21" ,   "MMP16" ,  "ADAMTS5" ,  "TIMP3" ,   "MMP2"   ,  "MMP14"   ,  "ADAMTS8",
                   "ADAMTS10"  ,  "MMP11"  , 
                   "ADAMTS9", "ADAMTS4","ADAMTS1" , "ADAMTS12",
                   "MMP27"  ,
                   
                   "MMP10",    "ADAMTS7" ,   "TIMP1"  , "MMP1"   , "MMP3"    , "MMP12"  , 
                   "MMP25" ,  "MMP9"  ,
                   "MMP28"  ,"MMP26"  ,      "MMP20"  , 
                   "ADAMTS17" ,   "MMP8"   ,"MMP7"  ,   "MMP15" ,"ADAMTS14"   )


pheatmap_result <- pheatmap(
  scaled_gene_expression[, ordered_cols1],
  cellwidth = 12,
  cellheight = 12,
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "none"  ,       # No scaling
  treeheight_col = 0
)


