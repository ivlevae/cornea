
library(Seurat)
library(SeuratDisk)

#SaveH5Seurat(cornea, paste0(main_path, 'cornea30.h5Seurat'), overwrite = TRUE)

#cornea30 <- LoadH5Seurat('/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/cornea30.h5Seurat')
#DimPlot(cornea30, label = T, repel = T, label.box = T, raster = T)

#table(cornea30$seurat_clusters)
#DefaultAssay(cornea30) <- 'RNA'

markers_subset <- c( 'ZNF90', 'MTRNR2L10',##activated fibroblastic stroma
                                        "MIR924HG", "DIAPH3", "BRIP1", "CENPP", #TAC
                                        'ACTA2', ##myofibroblasts
                                        "FBLN1",  "LEPR",  ##limbal fibroblasts 
                                        "POU6F2", "FAM155A", "CA3", ##corneal endo
                                        
                                        "PAX6", "KRT12", "TACSTD2", ##corneal epithelium
                                        "MUC16", "WFDC21P", "KRT24",  # 'Corneal Superficial'
                                        "HES5", "KRT3", "DIO2",  "KC6", "FOXP2", ##Corneal Wing
                                       "NKAIN2", "TENM2","LAMA3",  #'Corneal Basal' 
                                        "DCN", "LUM",   "KERA", "ABCA6", "ITGBL1" , "COL6A3", ##keratocytes
                                      
                                        'HOMER3', 'CPVL', 'TP63',##stem
                                        'AQP5', 'KRT13', 'IGFBP3', "KRT5", "KLF5", "AQP3") ## conjunctiva

new_cell_stem <- c( 
  'HOMER3', 'CPVL', 'TP63',
  'LY6D', 'KRT6A', 'FABP5','CSRP2','CXCL14',
  'S100A2','MT1X','NCOA7','SLC6A6',
  'SAA1','DCN','KERA','ANGPTL7','PTGDS','LUM','ITGBL1') 



# DimPlot(cornea30, reduction = "umap", raster = TRUE,  repel = T)
# DimPlot(cornea30, reduction = "umap", group.by = 'seurat_clusters', 
#         split.by = "seurat_clusters", raster = TRUE, ncol=6, pt.size = 4)  + NoLegend()
# table(cornea30$seurat_clusters)
# DimPlot(cornea30, label = T, repel = T, label.box = T, raster = T) + NoLegend()
# 
# limbal_progenitor_biomarkers <- c("p63", "ABCG2", "KRT14", "KRT15", "Notch Signaling", "LGR5", "CD34", "E-cadherin", "Sox2")
# DefaultAssay(cornea30) <- 'integrated'
# p <- DotPlot(cornea30,
#              features = c("GPHA2", "S100A2", "KRT15", "NCOA7", "SLC6A6"),
#              assay = NULL, cols = c("lightgrey", "blue"))
# p + theme(axis.text.x = element_text(angle = 90) )
# 
# p <- DotPlot(cornea30,
#              features = markers_subset,
#              assay = NULL, cols = c("lightgrey", "blue"))
# p + theme(axis.text.x = element_text(angle = 90) )
# 
# DefaultAssay(cornea30)

all_cell_markers_all_subset <- list(
  'Activated Keratocytes' = c('ZNF90', 'MTRNR2L10', 'HIST3H2A',  'MTRNR2L7',
                              'MTRNR2L6', 'PCOLCE', 'TIMP2', 'AL136454.1', 'HPR'), # activated fibroblastic stroma
  'TAC' = c("MIR924HG", "DIAPH3", "BRIP1", "CENPP", "ANLN", "ATAD2", "RIMS2", "POLQ", "MOXD1",  "TOP2A", "MKI67"),
  'Myofibroblasts' = c('ACTA2'),
  'Limbal Fibroblasts' = c("FBLN1", "MGP", "LEPR", "IGFBP5", "C1S", "CFD"),
  'Corneal Endothelium' = c("POU6F2", "NRXN3", "COL8A1", "SLC4A4", "COL4A3", "CA12", 
                            "CA3", "COL4A4"),
  'Keratocytes' = c( "DCN", "COL6A3", "CEMIP", "LSAMP", "ABCA6", "PTPRG", "MME",  "LUM", "KERA"),
  "Corneal epithelium" = c('CPAMD8', 'PCSK5', 'MYH14', 'ARHGEF4', "PAX6", "KRT12", "TACSTD2"), ##corneal epithelium
  
  'Corneal Wing' = c("HES5", "KRT3", "DIO2", "CSRP2", "GALNT18", "CACNA1E"),
  'Corneal Basal' = c("NKAIN2", "TENM2", "LAMA3", "CDH13", "IVNS1ABP", "MOXD1", "MIR205HG"),
  'Corneal Superficial' = c("MUC16", "WFDC21P", "KRT24", "MACC1", "BCAS1", "HOPX", "NECTIN4"),
  'Conjunctiva epithelium' = c("KRT13",  "KRT19", "S100A8", "S100A9", "KRT5", "KLF5", "AQP5", "AQP3"), ##conjunctiva epithelium

  'Limbal Basal' = c("KRT15", "LGR6", "PDGFC", "LAMA3"),
  'Limbal Wing' = c("KRT15", "KRT12",  "KC6", "FOXP2", "KRT5"),
  'limbal epithelium stem cells niche' = c("CAV1", "HOMER3", "CPVL" ),
  'Limbal_Neural_Crest_Progenitors' = c("CXCL14", "COL17A1",  "S100A2"),
  'Limbal_Progenitor_cells' = c("GPHA2", "S100A2", "KRT15", "NCOA7", "SLC6A6"),
  'Corneal_Stromal_Stem_Cells' = c("SAA1", "DCN", "KERA", "ANGPTL7", "PTGDS", "LUM", "ITGBL1"),
  'Stem Cells' = c("KRT15", "ITGB1",  "ITGB4", 'TP63' , 'CXCL14'),
  'Corneal_Stem_Stroma' = c('ABCG2', 'PAX6', 'NT5E', 'THY1', 'ENG'),
  'Corneal_Stem_Stroma_2' = c(
    "CD34",
    "THY1",     # for CD90
    "ENG",      # for CD105
    "SSEA4",    # for SSEA-4
    "KERA",     # for Keratocan
    "ALDH1A1",  # for ALDH
    "NES"       # for Nestin
  ))




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

# res2 <- AddScores(cornea30, all_cell_markers_all_subset)
# View(res2$data)
# p <- DotPlot(res2$data,
#              features = res2$features,
#              assay = NULL, cols = c("lightgrey", "blue"))
# p + theme(axis.text.x = element_text(angle = 90) )


cells_of_interest <- c(
  "Corneal limbal stem cells",
  "Corneal transit amplifying cells",
  "Corneal superficial epithelium",
  "Corneal basal epithelium",
  "Corneal wing epithelium",
  "Corneal regenerating epithelium",
  "Keratocytes",
  "Myofibroblasts",
  "Limbal fibroblasts",
  "Corneal endothelium"
)

# table(cornea30$seurat_clusters)
# cornea30@active.ident <- cornea30$seurat_clusters
# cornea30 <-RenameIdents(cornea30, '0' = 'Keratocytes',
#                                  '1' = 'Corneal wing epithelium',
#                                  '2' = 'Activated keratocytes',
#                                  '3' = 'Corneal basal epithelium',
#                                  '4' = 'Corneal wing epithelium',
#                                  '5' = 'Corneal superficial epithelium',
#                                  '6' = 'Conjunctiva + corneal epi',
#                                  '7' = 'Corneal wing epithelium',
#                                  '8' = 'Corneal basal epithelium',
#                                  '9' = 'Corneal wing epithelium',
#                                  '10' = 'Limbal fibroblasts',
#                                  '11' = 'Corneal endothelium',
#                                  '12' = 'Keratocytes',
#                                  '13' = 'Corneal wing epithelium',
#                                  '14' = 'Keratocytes +',
#                                  '15' = 'Myofibroblasts',
#                                  '16' = 'Keratocytes',
#                                  '17' = 'Corneal superficial epithelium',
#                                  '18' = 'TAC',
#                                  '19' = 'Limbal epithelium stem cells',
#                                  '20' = 'Keratocytes',
#                                  '21' = 'Corneal basal epithelium',
#                                  '22' = 'Conjunctiva',
#                                  '23' = 'Corneal wing epithelium',
#                                  '24' = 'Limbal fibroblasts'
#                                  
# )
# 
# cornea30$annot_V1 <- cornea30@active.ident
# DimPlot(cornea30, reduction = "umap", raster = TRUE, label = T, label.box = T,
#         repel = T) + NoLegend() 


ProcessInt <- function(data.integrated, npcs, res=1){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = npcs, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:npcs)
  data.integrated <- FindClusters(data.integrated, resolution = res)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:npcs)
}
