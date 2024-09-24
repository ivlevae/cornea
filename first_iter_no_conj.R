
cornea150_relabeled <- LoadH5Seurat('cornea150_res1_changed_annot_V2.h5Seurat')
cornea150_relabeled <- cornea150_res1
DefaultAssay(cornea150_relabeled) <- 'integrated'


DimPlot(cornea150_relabeled, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

cornea150_relabeled <- RenameIdents(cornea150_relabeled, 
                                   "Corneal endothelium/Keratocytes" = "Corneal endothelium")

cornea150_relabeled_filtered <- subset(cornea150_relabeled, idents = setdiff(levels(cornea150_relabeled), c('Conjunctiva', 'Melanocytes', 
                                                                                                   'Corneal endothelium/Keratocytes')))
cornea150_relabeled_filtered@meta.data$annot_V1 <- droplevels(cornea150_relabeled_filtered@meta.data$annot_V1)
cornea150_relabeled_filtered <- subset(cornea150_relabeled_filtered, idents = unique(cornea150_relabeled_filtered@meta.data$annot_V1))

DimPlot(cornea150_relabeled_filtered, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 




ProcessInt <- function(data.integrated, npcs, res=1){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = npcs, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:npcs)
  data.integrated <- FindClusters(data.integrated, resolution = res)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:npcs)
}

DefaultAssay(cornea150_relabeled_filtered) <- 'integrated'

cornea150_relabeled_filtered <- ProcessInt(cornea150_relabeled_filtered, 30, 1)
cornea150_relabeled_filtered_30 <- cornea150_relabeled_filtered


SaveH5Seurat(cornea150_relabeled_filtered_30,  'cornea150_relabeled_filtered_30.h5Seurat', overwrite = TRUE)


cornea150_relabeled_filtered_30_res3 <- ProcessInt(cornea150_relabeled_filtered, 30, 3)
cornea150_relabeled_filtered_30_res2 <- ProcessInt(cornea150_relabeled_filtered, 30, 2)
DimPlot(cornea150_relabeled_filtered_30, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
table(cornea150_relabeled_filtered_30_res2$seurat_clusters)
############# 30 PC

DefaultAssay(cornea150_relabeled_filtered_30) <- 'RNA'

source("subset_cell_markers.R")
library(ggplot2)

p <- DotPlot(cornea150_relabeled_filtered_30,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )




res2 <- AddScores(cornea150_relabeled_filtered_30, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

cornea150_relabeled_filtered_30@active.ident <- cornea150_relabeled_filtered_30$seurat_clusters
cornea150_relabeled_filtered_30 <-RenameIdents(cornea150_relabeled_filtered_30,    "0" = "Keratocytes",
                                               "1" = "Corneal wing epithelium",
                                               "2" = "Epi/kera",
                                               "3" = "Activated keratocytes",
                                               "4" = "Corneal basal epithelium",
                                               "5" = "Corneal basal epithelium",
                                               "6" = "Keratocytes",
                                               "7" = "Corneal superficial epithelium",
                                               "8" = "Limbal fibroblasts",
                                               "9" = "Corneal wing epithelium",
                                               "10" = "Corneal endothelium",
                                               "11" = "Corneal wing epithelium",
                                               "12" = "Corneal wing epithelium",
                                               "13" = "Keratocytes",
                                               "14" = "Corneal superficial epithelium",
                                               "15" = "Corneal wing epithelium/conjunctiva",
                                               "16" = "Myofibroblasts",
                                               "17" = "TAC",
                                               "18" = "Keratocytes",
                                               "19" = "Conj/Corneal basal epithelium",
                                               "20" = "Corneal wing epithelium",
                                               "21" = "Limbal stem cells",
                                               "22" = "Corneal superficial epithelium",
                                               "23" = "Limbal fibroblasts",
                                               "24" = "Corneal wing epithelium"
)
                                               
                                                                                 
cornea150_relabeled_filtered_30$annot <- cornea150_relabeled_filtered_30@active.ident
DimPlot(cornea150_relabeled_filtered_30, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
table(cornea150_relabeled_filtered_30$annot)

DimPlot(cornea150_relabeled_filtered_30, reduction = "umap", group.by = 'annot',
        split.by = "annot", raster = TRUE, ncol=6, pt.size = 4)  + NoLegend()

kera_epi <- subset(cornea150_relabeled_filtered_30, annot == 'Epi/kera')
table(kera_epi$annot_V1)


subset_3 <- subset(cornea150_relabeled_filtered_30, annot == 'Conj/Corneal basal epithelium')
table(subset_3$annot_V1)

subset_19 <- subset(cornea150_relabeled_filtered_30, annot == 'Corneal wing epithelium/conjunctiva')
table(subset_19$annot_V1)
#######

cornea150_relabeled_filtered_150 <- ProcessInt(cornea150_relabeled_filtered, 150, 1)

DimPlot(cornea150_relabeled_filtered_30, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

DimPlot(cornea150_relabeled_filtered_150, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


###################    19 - 'Corneal wing epithelium/conjunctiva'
DefaultAssay(subset_19) <- 'integrated'

subset_19_proccesed <- ProcessInt(subset_19, 150, 3)

dim(subset_19_proccesed)
DimPlot(subset_19_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

DefaultAssay(subset_19_proccesed) <- 'RNA'


p <- DotPlot(subset_19_proccesed,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(subset_19_proccesed, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

subset_19_proccesed@active.ident <- subset_19_proccesed$seurat_clusters

subset_19_proccesed <-RenameIdents(subset_19_proccesed,  "0" = "Corneal superficial epithelium",
                                   "1" = "Corneal wing epithelium",
                                   "2" = "Corneal wing epithelium",
                                   "3" = "Corneal superficial epithelium",
                                   "4" = "Corneal wing epithelium",
                                   "5" = "Corneal wing epithelium",
                                   "6" = "Corneal basal epithelium",
                                   "7" = "Corneal wing epithelium",
                                   "8" = "Corneal wing epithelium",
                                   "9" = "Conjunctiva",
                                   "10" = "Corneal wing epithelium",
                                   "11" = "Corneal wing epithelium",
                                   "12" = "Corneal wing epithelium",
                                   "13" = "Corneal wing epithelium",
                                   "14" = "Limbal fibroblasts",
                                   "15" = "Corneal wing epithelium",
                                   "16" = "Corneal wing epithelium"
)


subset_19_proccesed$annot<- subset_19_proccesed@active.ident

table(subset_19_proccesed$annot)

DimPlot(subset_19_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

cornea150_relabeled_filtered_30$annot_V1_1 <- cornea150_relabeled_filtered_30$annot
levels(cornea150_relabeled_filtered_30$annot_V1_1) <- c(levels(cornea150_relabeled_filtered_30$annot_V1_1), "Conjunctiva")



for (cell_type in names(table(subset_19_proccesed$annot))){
  print(cell_type)
  subsetObj <- subset(subset_19_proccesed, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_relabeled_filtered_30@meta.data[subset_cells, "annot_V1_1"] <- cell_type
  cornea150_relabeled_filtered_30@active.ident <- cornea150_relabeled_filtered_30$annot_V1_1
}





###################    subset_3 <- subset(cornea150_relabeled_filtered_30, annot == 'Conj/Corneal basal epithelium')

DefaultAssay(subset_3) <- 'integrated'
dim(subset_3)
subset_3_proccesed <- ProcessInt(subset_3, 150, 2)


table(subset_3_proccesed$seurat_clusters)
DimPlot(subset_3_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


DefaultAssay(subset_3_proccesed) <- 'RNA'


p <- DotPlot(subset_3_proccesed,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(subset_3_proccesed, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

subset_3_proccesed@active.ident <- subset_3_proccesed$seurat_clusters

subset_3_proccesed <-RenameIdents(subset_3_proccesed,    "0" = "Conjunctiva",
                                  "1" = "Corneal basal epithelium",
                                  "2" = "Corneal basal epithelium",
                                  "3" = "Corneal basal epithelium",
                                  "4" = "Conjunctiva",
                                  "5" = "Corneal basal epithelium",
                                  "6" = "Keratocytes",
                                  "7" = "Keratocytes",
                                  "8" = "Conjunctiva",
                                  "9" = "Corneal basal epithelium",
                                  "10" = "Corneal basal epithelium"
)


subset_3_proccesed$annot<- subset_3_proccesed@active.ident

table(subset_3_proccesed$annot)

DimPlot(subset_3_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


#levels(cornea150_relabeled_filtered_30$annot_V1_1) <- c(levels(cornea150_relabeled_filtered_30$annot_V1_1), "Conjunctiva")

#cornea150_relabeled_filtered_30$annot_V1_1 <- cornea150_relabeled_filtered_30$annot


for (cell_type in names(table(subset_3_proccesed$annot))){
  print(cell_type)
  subsetObj <- subset(subset_3_proccesed, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_relabeled_filtered_30@meta.data[subset_cells, "annot_V1_1"] <- cell_type
  cornea150_relabeled_filtered_30@active.ident <- cornea150_relabeled_filtered_30$annot_V1_1
}


DimPlot(cornea150_relabeled_filtered_30, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

table(cornea150_relabeled_filtered_30@active.ident)


