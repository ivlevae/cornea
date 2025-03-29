markers_subset <- list(
  'Corneal Wing' = c("HES5", "KRT3", "DIO2", "CSRP2", "GALNT18", "CACNA1E"),
  'Corneal Basal' = c("NKAIN2", "TENM2", "LAMA3", "CDH13", "IVNS1ABP", "MOXD1", "MIR205HG"),
  'Corneal Superficial' = c("MUC16", "WFDC21P", "KRT24", "MACC1", "BCAS1", "HOPX", "NECTIN4"),
  'Conjunctiva epithelium' = c("KRT13",  "KRT19", "S100A8", "S100A9"), ##conjunctiva epithelium
  'Cornea epithelium' = c("PAX6", "KRT12", "TACSTD2"), ##Cornea epithelium
  'Keratocytes' = c( "DCN", "COL6A3", "CEMIP", "LSAMP", "ABCA6", "PTPRG", "MME",  "LUM", "KERA") ##keratocytes
)

source("subset_cell_markers.R")

cornea150_res1_copy <- LoadH5Seurat('/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/cornea150_res1.h5Seurat')
DefaultAssay(cornea150_res1_copy)

table(cornea150_res1$seurat_clusters)
table(cornea150_res1$annot_V1)
DefaultAssay(cornea150_res1)
cornea150 <-RenameIdents(cornea150, '0' = 'Corneal endothelium',
                         '1' = 'Keratocytes',
                         '2' = 'Keratocytes + epi',
                         '3' = 'Limbal fibroblasts',
                         '4' = 'Corneal superficial epithelium',
                         '5' = 'Corneal wing epithelium',
                         '6' = 'Corneal wing epithelium',
                         '7' = 'Corneal basal epithelium',
                         '8' = 'Corneal basal epithelium',
                         '9' = 'Corneal wing epithelium',
                         '10' = 'Activated keratocytes',
                         '11' = 'Conjunctiva',
                         '12' = 'Corneal wing epithelium',
                         '13' = 'Keratocytes +',
                         '14' = 'Corneal wing epithelium',
                         '15' = 'Corneal wing epithelium',
                         '16' = 'Myofibroblasts + endo',
                         '17' = 'Keratocytes', 
                         '18' = 'TAC',
                         '19' = 'Melanocytes')


###################    '2' = 'Keratocytes + epi',

subset_1 <- subset(cornea150_res1, seurat_clusters == '2')
dim(subset_1)
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
subset_2_proccesed <- ProcessInt(subset_1)
DimPlot(subset_2_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DefaultAssay(subset_2_proccesed) <- 'RNA'


p <- DotPlot(subset_2_proccesed,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(subset_2_proccesed, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

subset_2_proccesed@active.ident <- subset_2_proccesed$seurat_clusters
subset_2_proccesed <-RenameIdents(subset_2_proccesed,  "0" = "Corneal basal epithelium",
"1" = "Keratocytes",
"2" = "Corneal basal epithelium",
"3" = "Corneal basal epithelium",
"4" = "Keratocytes",
"5" = "Corneal basal epithelium",
"6" = "Corneal basal epithelium",
"7" = "Kera/Epi",
"8" = "Kera/Epi",
"9" = "Corneal basal epithelium",
"10" = "TAC",
"11" = "Corneal basal epithelium",
"12" = "Corneal basal epithelium",
"13" = "Corneal basal epithelium",
"14" = "Corneal basal epithelium",
"15" = "Corneal basal epithelium",
"16" = "Corneal basal epithelium",
"17" = "Limbal fibroblasts",
"18" = "Conjunctiva",
"19" = "Keratocytes",
"20" = "Keratocytes"
)
subset_2_proccesed$annot<- subset_2_proccesed@active.ident


DimPlot(subset_2_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

for (cell_type in names(table(subset_2_proccesed$annot))){
  print(cell_type)
  subsetObj <- subset(subset_2_proccesed, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_res1@meta.data[subset_cells, "annot_V1"] <- cell_type
  cornea150_res1@active.ident <- cornea150_res1$annot_V1
}

########## '13' keratocytes +
subset_13 <- subset(cornea150_res1, seurat_clusters == 13)
dim(subset_13)
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
subset13_proccesed <- ProcessInt(subset_13)
DimPlot(subset13_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DefaultAssay(subset13_proccesed) <- 'RNA'


p <- DotPlot(subset13_proccesed,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )





res2 <- AddScores(subset13_proccesed, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

subset13_proccesed@active.ident <- subset13_proccesed$seurat_clusters
subset13_proccesed <-RenameIdents(subset13_proccesed,  "0" = "Keratocytes",
                                  "1" = "Keratocytes",
                                  "2" = "Keratocytes",
                                  "3" = "Keratocytes",
                                  "4" = "Kera/Epi",
                                  "5" = "Corneal basal epithelium",
                                  "6" = "Corneal basal epithelium",
                                  "7" = "Conjunctiva",
                                  "8" = "Limbal fibroblasts",
                                  "9" = "Kera/Epi",
                                  "10" = "Corneal basal epithelium",
                                  "11" = "Keratocytes",
                                  "12" = "Limbal fibroblasts",
                                  "13" = "Keratocytes",
                                  "14" = "Corneal basal epithelium",
                                  "15" = "Keratocytes",
                                  "16" = "Corneal basal epithelium",
                                  "17" = "Keratocytes",
                                  "18" = "Keratocytes"

)
subset13_proccesed$annot<- subset13_proccesed@active.ident
DimPlot(subset13_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
table(subset13_proccesed$annot)

levels(cornea150_res1$annot_V1) <- c(levels(cornea150_res1$annot_V1), "Kera/Epi")
for (cell_type in names(table(subset13_proccesed$annot))){
  print(cell_type)
  subsetObj <- subset(subset13_proccesed, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_res1@meta.data[subset_cells, "annot_V1"] <- cell_type
  cornea150_res1@active.ident <- cornea150_res1$annot_V1
}

###    16 Myofibroblasts + endo

subset_16 <- subset(cornea150_res1, seurat_clusters == 16)
dim(subset_16)

subset16_proccesed <- ProcessInt(subset_16)
DimPlot(subset16_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DefaultAssay(subset16_proccesed) <- 'RNA'


p <- DotPlot(subset16_proccesed,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

res2 <- AddScores(subset16_proccesed, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

subset16_proccesed@active.ident <- subset16_proccesed$seurat_clusters
subset16_proccesed <-RenameIdents(subset16_proccesed,  "0" = "Myofibroblasts",
                                  "1" = "Limbal fibroblasts",
                                  "2" = "Myofibroblasts",
                                  "3" = "Kera/Epi",
                                  "4" = "Corneal endothelium",
                                  "5" = "Corneal endothelium",
                                  "6" = "Keratocytes",
                                  "7" = "Myofibroblasts",
                                  "8" = "Limbal fibroblasts",
                                  "9" = "Corneal basal epithelium",
                                  "10" = "Keratocytes")

subset16_proccesed$annot<- subset16_proccesed@active.ident
DimPlot(subset16_proccesed, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
table(subset16_proccesed$seurat_clusters)

for (cell_type in names(table(subset16_proccesed$annot))){
  print(cell_type)
  subsetObj <- subset(subset16_proccesed, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_res1@meta.data[subset_cells, "annot_V1"] <- cell_type
  cornea150_res1@active.ident <- cornea150_res1$annot_V1
}



###################    '0' - endothelium + smth 
DefaultAssay(cornea150_res1) <- 'integrated'
subset_0 <- subset(cornea150_res1, seurat_clusters == '0')
dim(subset_0)
ProcessInt_150 <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 3)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}


subset_0_proccesed_150_res3 <- ProcessInt_150(subset_0)
DimPlot(subset_0_proccesed_159, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DimPlot(subset_0_proccesed_150_res3, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DefaultAssay(subset_0_proccesed_150_res3) <- 'RNA'


p <- DotPlot(subset_0_proccesed_150_res3,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(subset_0_proccesed_150_res3, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

subset_0_proccesed_150_res3@active.ident <- subset_0_proccesed_150_res3$seurat_clusters

subset_0_proccesed_150_res3 <-RenameIdents(subset_0_proccesed_150_res3,  "0" = "Keratocytes",
                                  "1" = "Kera/Epi",
                                  "2" = "Corneal endothelium",
                                  "3" = "Keratocytes/endothelium",
                                  "4" = "Keratocytes",
                                  "5" = "Keratocytes",
                                  "6" = "Keratocytes",
                                  "7" = "Keratocytes",
                                  "8" = "Keratocytes/endothelium",
                                  "9" = "Keratocytes",
                                  "10" = "Keratocytes",
                                  "11" = "Keratocytes",
                                  "12" = "Limbal fibroblasts",
                                  "13" = "Limbal fibroblasts",
                                  "14" = "Keratocytes/endothelium",
                                  "15" = "Limbal fibroblasts",
                                  "16" = "Corneal endothelium",
                                  "17" = "Corneal endothelium",
                                  "18" = "Kera/Epi",
                                  "19" = "Kera/Epi"
)
                                  
                                  
                                  
subset_0_proccesed_150_res3$annot<- subset_0_proccesed_150_res3@active.ident

table(subset_0_proccesed_150_res3$annot)

DimPlot(subset_0_proccesed_150_res3, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
levels(cornea150_res1$annot_V1) <- c(levels(cornea150_res1$annot_V1), "Keratocytes/endothelium")

for (cell_type in names(table(subset_0_proccesed_150_res3$annot))){
  print(cell_type)
  subsetObj <- subset(subset_0_proccesed_150_res3, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_res1@meta.data[subset_cells, "annot_V1"] <- cell_type
  cornea150_res1@active.ident <- cornea150_res1$annot_V1
}





######################## kera/endo

DefaultAssay(cornea150_res1) <- 'integrated'
subset_kera_endo <- subset(cornea150_res1, annot_V1 == 'Keratocytes/endothelium')
dim(subset_kera_endo)
ProcessInt_150 <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}


subset_kera_endo_res3 <- ProcessInt_150(subset_kera_endo)
DimPlot(subset_kera_endo_res3, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

DefaultAssay(subset_kera_endo_res3) <- 'integrated'


p <- DotPlot(subset_kera_endo_res3,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(subset_kera_endo_res3, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


table(subset_kera_endo_res3$seurat_clusters)
subset_kera_endo_res3@active.ident <- subset_kera_endo_res3$seurat_clusters


subset_kera_endo_res3 <-RenameIdents(subset_kera_endo_res3,  
                                           "0" = "Keratocytes",
                                           "1" = "Corneal endothelium",
                                           "2" = "Corneal endothelium/Keratocytes",
                                           "3" = "Corneal endothelium/Keratocytes"
)



subset_kera_endo_res3$annot<- subset_kera_endo_res3@active.ident

table(subset_kera_endo_res3$annot)

DimPlot(subset_kera_endo_res3, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


levels(cornea150_res1$annot_V1) <- c(levels(cornea150_res1$annot_V1), "Corneal endothelium/Keratocytes")

for (cell_type in names(table(subset_kera_endo_res3$annot))){
  print(cell_type)
  subsetObj <- subset(subset_kera_endo_res3, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_res1@meta.data[subset_cells, "annot_V1"] <- cell_type
  cornea150_res1@active.ident <- cornea150_res1$annot_V1
}


########## check the whole UMAP

table(cornea150_res1$annot_V1)
dim(cornea150_res1)

DimPlot(cornea150_res1, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DimPlot(cornea150_res1, reduction = "umap", group.by = 'annot_V1',
        split.by = "annot_V1", raster = TRUE, ncol=6, pt.size = 4)  + NoLegend()

DefaultAssay(cornea150_res1) <- 'RNA'

p <- DotPlot(cornea150_res1,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(cornea150_res1, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


table(cornea150_res1$annot_V1)

DimPlot(cornea150_res1, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DimPlot(cornea150_res1_copy, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

DimPlot(cornea150_res1, reduction = "umap", group.by = 'annot_V1',
        split.by = "annot_V1", raster = TRUE, ncol=6, pt.size = 4)  + NoLegend()
DimPlot(cornea150_res1_copy, reduction = "umap", group.by = 'annot_V1',
        split.by = "annot_V1", raster = TRUE, ncol=6, pt.size = 4)  + NoLegend()


SaveH5Seurat(cornea150_res1, 'cornea150_res1_changed_annot.h5Seurat', overwrite = TRUE)

SaveH5Seurat(subset_2_proccesed, 'subset_2_proccesed.h5Seurat', overwrite = TRUE)
SaveH5Seurat(subset13_proccesed, 'subset13_proccesed.h5Seurat', overwrite = TRUE)
SaveH5Seurat(subset16_proccesed, 'subset16_proccesed.h5Seurat', overwrite = TRUE)

DefaultAssay(cornea150_res1) <- 'RNA'
p <- DotPlot(cornea150_res1,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )

res2 <- AddScores(cornea150_res1, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


######################## kera/endo

DefaultAssay(cornea150_res1) <- 'integrated'
subset_kera_epi <- subset(cornea150_res1, annot_V1 == 'Kera/Epi')
dim(subset_kera_epi)
ProcessInt_150 <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 3)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}


subset_kera_epi_res3<- ProcessInt_150(subset_kera_epi)
DimPlot(subset_kera_epi_res3, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 

DefaultAssay(subset_kera_epi_res3) <- 'RNA'


p <- DotPlot(subset_kera_epi_res3,
             features = markers_subset,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


res2 <- AddScores(subset_kera_epi_res3, all_cell_markers_all_subset)

p <- DotPlot(res2$data,
             features = res2$features,
             assay = NULL, cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )


table(subset_kera_epi_res3$seurat_clusters)
subset_kera_epi_res3@active.ident <- subset_kera_epi_res3$seurat_clusters


subset_kera_epi_res3 <-RenameIdents(subset_kera_epi_res3,  
                                    "0" = "Corneal wing epithelium",
                                    "1" = "Corneal endothelium",
                                    "2" = "Corneal wing epithelium",
                                    "3" = "Corneal basal epithelium",
                                    "4" = "Limbal fibroblasts",
                                    "5" = "Keratocytes",
                                    "6" = "Corneal superficial epithelium",
                                    "7" = "Corneal basal epithelium",
                                    "8" = "Myofibroblasts",
                                    "9" = "Conjunctiva",
                                    "10" = "Keratocytes"
)



subset_kera_epi_res3$annot<- subset_kera_epi_res3@active.ident

table(subset_kera_epi_res3$annot)

DimPlot(subset_kera_epi_res3, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


levels(cornea150_res1$annot_V1) <- c(levels(cornea150_res1$annot_V1), "Corneal endothelium/Keratocytes")

for (cell_type in names(table(subset_kera_epi_res3$annot))){
  print(cell_type)
  subsetObj <- subset(subset_kera_epi_res3, annot %in% c(cell_type))
  print(dim(subsetObj))
  subset_cells <- Cells(subsetObj)
  cornea150_res1@meta.data[subset_cells, "annot_V1"] <- cell_type
  cornea150_res1@active.ident <- cornea150_res1$annot_V1
}





########## check the whole UMAP

table(cornea150_res1$annot_V1)
dim(cornea150_res1)

DimPlot(cornea150_res1, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 
DimPlot(cornea150_res1, reduction = "umap", group.by = 'annot_V1',
        split.by = "annot_V1", raster = TRUE, ncol=6, pt.size = 4)  + NoLegend()


SaveH5Seurat(cornea150_res1, 'cornea150_res1_changed_annot_V2.h5Seurat', overwrite = TRUE)
