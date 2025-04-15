
###################### cutting of Conjunctiva cells left from the previous step ############################

cornea150_relabeled_filtered_30_2 <- subset(cornea150_relabeled_filtered_30, idents = setdiff(levels(cornea150_relabeled_filtered_30), c('Conjunctiva')))
cornea150_relabeled_filtered_30_2@meta.data$annot_V1_1 <- droplevels(cornea150_relabeled_filtered_30_2@meta.data$annot_V1_1)
cornea150_relabeled_filtered_30_2 <- subset(cornea150_relabeled_filtered_30_2, idents = unique(cornea150_relabeled_filtered_30_2@meta.data$annot_V1_1))

DimPlot(cornea150_relabeled_filtered_30_2, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


SaveH5Seurat(cornea150_relabeled_filtered_30_2,  'cornea150_relabeled_filtered_30_filtered.h5Seurat', overwrite = TRUE)


DefaultAssay(cornea150_relabeled_filtered_30_2) <- 'integrated'


ProcessInt <- function(data.integrated, npcs, res=1){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = npcs, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:npcs)
  data.integrated <- FindClusters(data.integrated, resolution = res)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:npcs)
}


################ trying 150 and 30 PCs ###################
cornea150_relabeled_filtered_30_2_30 <- ProcessInt(cornea150_relabeled_filtered_30_2, 30, 1)
cornea150_relabeled_filtered_30_2_150 <- ProcessInt(cornea150_relabeled_filtered_30_2, 150, 1)


DimPlot(cornea150_relabeled_filtered_30_2_30, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

DimPlot(cornea150_relabeled_filtered_30_2_150, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


################### 30 PC was chosen for the future analysis #######################################

DefaultAssay(cornea150_relabeled_filtered_30_2_30) <- 'RNA'


p <- DotPlot(cornea150_relabeled_filtered_30_2_30,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )





markers_subset_2 <- c( 
                      # "PAX6", "KRT12", "TACSTD2", ##corneal epithelium
                       "MUC16", "WFDC21P", "KRT24",  "MACC1", "BCAS1", "HOPX", "NECTIN4",  # 'Corneal Superficial'
                       "HES5", "KRT3", "DIO2",  "KC6", "FOXP2", "GALNT18", "CACNA1E", ##Corneal Wing
                       "NKAIN2", "TENM2","LAMA3",  "CDH13", "IVNS1ABP", "MIR205HG", #'Corneal Basal' 
                       "MIR924HG", "DIAPH3", "BRIP1", "CENPP", "ANLN", "ATAD2", "RIMS2", "POLQ",  "TOP2A", "MKI67", 'S100A2',#TAC
                      'CXCL14', 
                      "CSRP2", 
                       "KRT15", "LGR6",  "KRT14",  #"CLDN4",
                      "GJA1", ## limbal basal
                      "PLAUR", 
                      "DCN", "LUM",   "KERA", "ABCA6", "ITGBL1" , "COL6A3",  ##keratocytes
                      'TAGLN', 'RGS5', 'TPM2', 'TPM1', 'ID4', 'SPARCL1', 'NR2F2',
                      
                      'MYL9', 'MYH11', 'ACTA2', 'MYLK' , ### myofibroblasts
                      "FBLN1", 'COL3A1',  "COL1A1", 'COL1A2','SOCS3','SRPX',  'SFRP2',  'ELN','MMP2', ##limbal fibroblasts 
                       "POU6F2", "FAM155A", "CA3" ##corneal endo
                   )             


p <- DotPlot(cornea150_relabeled_filtered_30_2_30_1,
             features = markers_subset_2,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )





cornea150_relabeled_filtered_30_2_30_1@active.ident <- cornea150_relabeled_filtered_30_2_30_1$seurat_clusters
cornea150_relabeled_filtered_30_2_30_1 <-RenameIdents(cornea150_relabeled_filtered_30_2_30_1,    
                                      '0' = 'Keratocytes',
                                      '1' = 'Corneal Wing',
                                      '2' = 'Corneal Basal',
                                      '3' = 'Keratocytes',
                                      '4' = 'Limbal Basal',
                                      '5' = 'Corneal Basal',
                                      '6' = 'Corneal Superficial',
                                      '7' = 'Corneal Wing',
                                      '8' = 'Corneal Wing',
                                      '9' = 'Keratocytes',
                                      '10' = 'Endothelium',
                                      '11' = 'Keratocytes',
                                      '12' = 'Corneal Wing',
                                      '13' = 'Corneal Superficial',
                                      '14' = 'Keratocytes',
                                      '15' = 'Keratocytes',
                                      '16' = 'Corneal Superficial',
                                      '17' = 'Myofibroblasts',
                                      '18' = 'Keratocytes',
                                      '19' = 'TAC',
                                      '20' = 'Corneal Superficial',
                                      '21' = 'Limbal Basal',
                                      '22' = 'Limbal Basal',
                                      '23' = 'Limbal Basal',
                                      '24' = 'Limbal Suprabasal',
                                      '25' = 'Limbal fibroblasts'
                                      )

cornea150_relabeled_filtered_30_2_30_1$annot_iter4 <- cornea150_relabeled_filtered_30_2_30_1@active.ident
SaveH5Seurat(cornea150_relabeled_filtered_30_2_30_1,  'cornea150_relabeled_filtered_30_2_30_1_iter4.h5Seurat', overwrite = TRUE)





