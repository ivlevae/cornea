### SeuObjProcessedList was a dictionary 
### Active assay: RNA in each Seurat object

library(dplyr)

# Convert the (named) list of processed Seurat objects into an unnamed list
integration_list <- unname(unlist(SeuObjProcessedList))

# Sanity checks on input list
print(SeuObjProcessedList)
print(length(SeuObjProcessedList))  
print(names(SeuObjProcessedList))

# Garbage collection to free memory
gc()

# ---------------------------
# Function 1: Integration post-processing (30 PCs)
# Scales data, regresses out confounders, runs PCA, neighbors, clusters, UMAP
# ---------------------------
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(
    data.integrated,
    verbose = TRUE,
    vars.to.regress = c("percent.mt", "percent.rb", "S.Score", "G2M.Score")
  )
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = TRUE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
  return(data.integrated)
}

# Select integration features from the list of Seurat objects
features <- SelectIntegrationFeatures(
  object.list = integration_list,
  nfeatures = 2000
)

gc()

# Find integration anchors across all objects
data.anchors <- FindIntegrationAnchors(
  object.list = integration_list,
  anchor.features = features
)
data.anchors

gc()

# Integrate data into a single Seurat object
cornea <- IntegrateData(anchorset = data.anchors)
gc()

# Run the 30-PC processing pipeline on the integrated object
cornea <- ProcessInt(cornea)

# ---------------------------
# Function 2: Re-define ProcessInt for 150 PCs (no scaling/regression here)
# Note: this overrides the previous definition of ProcessInt.
# ---------------------------
ProcessInt <- function(data.integrated){
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = TRUE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
  return(data.integrated)
}

# Set default assay for downstream analyses
DefaultAssay(cornea30)  <- "integrated"
DefaultAssay(cornea150) <- "integrated"

# Increase allowed future size for parallelization (if using future-based methods)
options(future.globals.maxSize = 100000 * 1024^2)

# Apply 150-PC pipeline to the integrated object
# (result stored in cornea150; assumes `cornea` already exists from above)
cornea150 <- ProcessInt(cornea)

# ---------------------------
# Marker discovery for 30-PC solution
# ---------------------------
markerscornea <- FindAllMarkers(
  cornea30,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markerscornea30_top25 <- markerscornea %>%
  group_by(cluster) %>%
  slice_max(n = 25, order_by = avg_log2FC)

gc()

# ---------------------------
# Marker discovery for 150-PC solution
# ---------------------------
markerscornea150 <- FindAllMarkers(
  cornea150,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# UMAP plot of 150-PC clustering
DimPlot(
  cornea150,
  label = TRUE,
  repel = TRUE,
  label.box = TRUE,
  raster = TRUE
)

# Quick overview of cluster sizes
table(cornea150$seurat_clusters)

markerscornea150_top25 <- markerscornea150 %>%
  group_by(cluster) %>%
  slice_max(n = 25, order_by = avg_log2FC)

View(markerscornea150_top25)

# Save marker tables
write.csv(
  markerscornea30_top25,
  "/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/cornea30_markers.csv"
)
write.csv(
  markerscornea150,
  "cornea150_markers_all.csv"
)

# Switch default assay back to RNA for downstream expression-based work
DefaultAssay(cornea30)  <- "RNA"
DefaultAssay(cornea150) <- "RNA"

# ---------------------------
# Save Seurat objects (HDF5 and RDS formats)
# ---------------------------
library(SeuratDisk)
SaveH5Seurat(cornea150, "cornea150.h5Seurat", overwrite = TRUE)
SaveH5Seurat(cornea,    "cornea30.h5Seurat", overwrite = TRUE)

saveRDS(cornea150, "cornea150.rds")
saveRDS(cornea30,  "cornea30.rds")

# ---------------------------
# Function 3: Re-define ProcessInt again for resolution = 2 (still 150 PCs)
# This overrides the previous ProcessInt definition.
# ---------------------------
ProcessInt <- function(data.integrated){
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = TRUE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 2)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
  return(data.integrated)
}

# Re-run processing with higher clustering resolution (2.0)
# Note: here you pass cornea30 into a function tuned to 150 PCs
cornea150_res2 <- ProcessInt(cornea30)

# Plot UMAP of the higher-resolution clustering
DimPlot(
  cornea150_res2,
  label = TRUE,
  repel = TRUE,
  label.box = TRUE,
  raster = TRUE
)
