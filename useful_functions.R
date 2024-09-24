library(Seurat)


AddScores <- function(dataset, cell_dict){
  # Function that calculates Module Scores and add it to the Seurat object
  # input:
  # dataset - input Seurat Object
  # cell_dict - dictionary of cell biomarkers 
  #                           like  list(cell1 = c(gene1, gene2), cell2 = c(gene3, gene4))
  # output - list with calculated scores for different cell types
  
  dataset_score <- dataset
  for (cell_marker in names(cell_dict)) {
    print(cell_marker)
    safe_name <- gsub(" ", "_", cell_marker)
    dataset_score <- AddModuleScore(  object = dataset_score,
                                      features = list(cell_dict[[cell_marker]]),
                                      ctrl = 5,
                                      name = safe_name)
    
  }  

  n_cols <- length(names(cell_dict))
  n_all_cols <- length(names(dataset_score@meta.data))
  n_start <- n_all_cols-n_cols+1
  cols_names <-colnames(dataset_score@meta.data[n_start:n_all_cols])
  return_list <- list("data" = dataset_score, "features" = cols_names)
  return(return_list) 
}  


##### How to use #####
# res <- AddScores(cornea30_noconj, subset_cell_markers) 
# p <- DotPlot(res$data, 
#              features = res$features,
#              assay = NULL, cols = c("lightgrey", "blue")) 
# p + theme(axis.text.x = element_text(angle = 90) ) 




ProcessInt <- function(data.integrated, nPCs, res = 1) {
  if (nPCs <= 30) {
    data.integrated <- ScaleData(data.integrated, vars.to.regress = c('percent.mt', 'percent.rb', 'S.Score', 'G2M.Score'), verbose = TRUE)
    data.integrated <- RunPCA(data.integrated, npcs = nPCs, verbose = TRUE)
  } else {
    data.integrated <- RunPCA(data.integrated, npcs = nPCs, verbose = TRUE, vars.to.regress = c('percent.mt', 'percent.rb', 'S.Score', 'G2M.Score'))
  }
  
  data.integrated <- FindNeighbors(data.integrated, dims = 1:nPCs)
  data.integrated <- FindClusters(data.integrated, resolution = res)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:nPCs)
  
  return(data.integrated)
}

##### How to use
# SeuObjIntegrated <- ProcessInt(SeuObj, 30, 1) 





  