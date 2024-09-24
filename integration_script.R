### SeuObjProcessedList was a dictionary 
### Active assay RNA 
library(dplyr)
integration_list <- unname(unlist(SeuObjProcessedList))
print(SeuObjProcessedList)
print(length(SeuObjProcessedList))  
print(names(SeuObjProcessedList))

gc()
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}

features <- SelectIntegrationFeatures(object.list = integration_list, nfeatures = 2000)
gc()
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.anchors
data.anchorsrpca <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features, reduction = 'rpca')
data.anchorsrpca
gc()
cornea <- IntegrateData(anchorset = data.anchors)
gc()
cornea <- ProcessInt(cornea)

ProcessInt <- function(data.integrated){
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}


DefaultAssay(cornea30) <- 'integrated'
DefaultAssay(cornea150) <- 'integrated'

options(future.globals.maxSize = 100000 * 1024^2)
cornea150 <- ProcessInt(cornea)

markerscornea <- FindAllMarkers(cornea30, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
markerscornea30_top25 <- markerscornea %>%
  group_by(cluster) %>%
  slice_max(n=25, order_by = avg_log2FC)

gc()
markerscornea150 <- FindAllMarkers(cornea150, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
DimPlot(cornea150, label = T, repel = T, label.box = T, raster = T)

table(cornea150$seurat_clusters)
markerscornea150_top25 <- markerscornea150 %>%
  group_by(cluster) %>%
  slice_max(n=25, order_by = avg_log2FC)

View(markerscornea150_top25)

write.csv(markerscornea30_top25, '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/cornea30_markers.csv')
write.csv(markerscornea150, 'cornea150_markers_all.csv')

DefaultAssay(cornea30) <- 'RNA'
DefaultAssay(cornea150) <- 'RNA'

library(SeuratDisk)
SaveH5Seurat(cornea150, 'cornea150.h5Seurat', overwrite = TRUE)
SaveH5Seurat(cornea, 'cornea30.h5Seurat', overwrite = TRUE)
saveRDS(cornea150, 'cornea150.rds')
saveRDS(cornea30, 'cornea30.rds')


ProcessInt <- function(data.integrated){
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 2)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}

cornea150_res2 <- ProcessInt(cornea30)
DimPlot(cornea150_res2, label = T, repel = T, label.box = T, raster = T)
