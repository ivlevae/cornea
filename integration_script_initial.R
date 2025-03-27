## SeuObjProcessedList - list of Seurat objects
### Active assay RNA 
integration_list <- unname(unlist(SeuObjProcessedList))

gc()
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T) #, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score")
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
features <- SelectIntegrationFeatures(object.list = integration_list, nfeatures = 2000)
gc()
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.anchors
data.anchorsrpca <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features, reduction = 'cca')
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
DefaultAssay(cornea) <- 'integrated'
options(future.globals.maxSize = 100000 * 1024^2)
cornea150 <- ProcessInt(cornea)
markerscornea <- FindAllMarkers(cornea, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
gc()
markerscornea150 <- FindAllMarkers(cornea150, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
DimPlot(cornea150, label = T, repel = T, label.box = T, raster = T)
test <- markerscornea150 %>%
  group_by(cluster) %>%
  slice_max(n=25, order_by = avg_log2FC)

write.csv(test, '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/cornea150_markers.csv')
DefaultAssay(cornea) <- 'RNA'
DefaultAssay(cornea150) <- 'RNA'

library(SeuratDisk)
SaveH5Seurat(cornea150, '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/cornea150.h5Seurat', overwrite = TRUE)
SaveH5Seurat(cornea, '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/cornea30.h5Seurat', overwrite = TRUE)
saveRDS(cornea150, '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/cornea150.rds')
saveRDS(cornea, '/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset/cornea30.rds')