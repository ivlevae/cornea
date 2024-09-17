library(Seurat)
library(SeuratDisk)
library(dplyr)
cornea150 <- LoadH5Seurat("/home/bnvlab2/Documents/Kate/cornea150.h5Seurat")

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


cornea150_annot_V3 <-RenameIdents(cornea150, '0'	= 'Corneal wing epithelium',
                                  '2'	= 'Keratocytes'	,
                                  '4'	= 'Corneal basal epithelium'	,
                                  '12'	= 'Keratocytes',
                                  '18' =	'Mitochondrial enriched cells',
                                  '19'=	'Corneal regenerating epithelium',
                                  '20' =	'Keratocytes',
                                  '21'	=	'Blood vessel endothelium',
                                  '23'	= 'Myofibroblasts',
                                  '25'	= 'Corneal transit amplifying cells',
                                  '26'	 = 'Myeloid',
                                  '29' =	'Lymphatic vessel endothelium'	,
                                  '30' = 'Corneal basal epithelium'	,
                                  '27' =	'Lymphocytes',
                                  '1'	= 'Conjunctival basal epithelium',
                                  '3' =	'Corneal basal epithelium',
                                  '5'	= 'Keratocytes' ,
                                  '6'	= 'Corneal superficial epithelium'	,
                                  '7'	= 'Corneal wing epithelium',
                                  '8'	 = 'Corneal basal epithelium',
                                  '9'	= 'Conjunctival non-basal epithelium/wing',
                                  '10' =	'Keratocytes',
                                  '11'	= 'Conjunctival basal epithelium' ,
                                  '13'	= 'Limbal fibroblasts'	,
                                  '14'	= 'Corneal wing epithelium',
                                  '15'	= 'Keratocytes'	, #####was excluded in further analysis 5673 
                                  '16' = 'Corneal limbal stem cells',
                                  '17'	= 'Corneal endothelium',
                                  '22' = 'Conjunctiva melanocytes'	,
                                  '24' =	'Conjunctival basal epithelium'	,
                                  '28'	= 'Non-myelinating Schwann cells'
)

new_cell_markers_best_list_subset <- c( 'ZNF90', 'MTRNR2L10',##activated fibroblastic stroma
                                        "MIR924HG", "DIAPH3", "BRIP1", "CENPP", #TAC
                                        'ACTA2', ##myofibroblasts
                                        "FBLN1",  "LEPR",  ##limbal fibroblasts 
                                        "POU6F2", "FAM155A", "CA3", ##corneal endo
                                        "PAX6", "KRT12", "TACSTD2", ##corneal epithelium
                                        "DCN", "LUM",   "KERA", "ABCA6", "ITGBL1" , "COL6A3", ##keratocytes
                                        "COL6A1", "COL12A1",
                                        'HOMER3', 'CPVL', 'TP63',##stem
                                        'AQP5', 'KRT13', 'IGFBP3') 
DefaultAssay(cornea150_annot_V3) <- 'RNA'
p <- DotPlot(cornea150_annot_V3, 
             features = new_cell_markers_best_list_subset,
             assay = NULL, cols = c("lightgrey", "blue")) 
p + theme(axis.text.x = element_text(angle = 90) ) 
table(cornea150_annot_V3$annot_V3)

DimPlot(cornea150_annot_V3, reduction = "umap", raster = TRUE,
        label = T, label.box = T,     repel = T) + NoLegend()

cornea150_annot_V3$annot_V3 <- cornea150_annot_V3@active.ident
table(cornea150$annot_V3)

metadata <- cornea150_annot_V3@meta.data
cluster_batch_summary <- metadata %>%
  group_by(annot_V3, batch) %>%
  summarize(cell_count = n())
batch_totals <- metadata %>%
  group_by(batch) %>%
  summarize(total_cells = n())

batch_totals_ncells <- metadata %>%
  filter(annot_V3 %in% cells_of_interest) %>%
  group_by(batch) %>%
  summarize(total_cells = n())

View(summary_with_ncells_sorted)
summary_with_ratios <- cluster_batch_summary %>%
  left_join(batch_totals, by = "batch") %>%
  mutate(ratio = cell_count / total_cells*100)
summary_with_ratios_sorted <-
  summary_with_ratios %>% arrange(batch, desc(ratio))  %>%
  select(annot_V3, batch, ratio, cell_count)
summary_with_ncells_sorted <-
  summary_with_ratios %>% arrange(batch, desc(cell_count))  %>%
  select(annot_V3, batch, cell_count)

View(summary_with_ncells_sorted)

summary_with_ratios_sorted_selected <- subset(summary_with_ratios_sorted, annot_V3 %in% cells_of_interest)
# Create the bar plot
ggplot(summary_with_ratios_sorted_selected, aes(x = batch, y = ratio, fill = annot_V3)) +
  geom_bar(stat = "identity", color = "black") +  # Create bars with black borders
  theme_minimal() +  # Use a minimal theme for clarity
  labs(x = "Category", y = "Ratio", fill = "Annot V3") +  # Label the axes and legend
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    legend.position = "right"  # Position legend to the right
  )


cornea150_annot_V3_selected <- subset(cornea150_annot_V3,  annot_V3 %in% cells_of_interest)
dim(cornea150_annot_V3)
dim(cornea150_annot_V3_selected)


ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 3000)
  Seurat <- ScaleData(Seurat, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
  Seurat <- RunPCA(Seurat, npcs = 100)
  Seurat <- FindNeighbors(Seurat, dims = 1:100)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:100)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

SeuObjProcessedList <- list()
SeuObjList <- list()

table(cornea150_annot_V3_selected$batch)

DefaultAssay(cornea150_annot_V3_selected) <- 'RNA'

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

for (batch_id in unique(cornea150_annot_V3_selected$batch)) {
  print(batch_id)
  SeuObjBatch <- subset(cornea150_annot_V3_selected, batch == batch_id)
  SeuObjList[[batch_id]] <- SeuObjBatch
  ## "Batch5", "Batch6", "Batch7", "Batch9","Batch10", "Batch11" are missing in final version 
  if (!batch_id %in% c("Batch3", "Batch5", "Batch6", "Batch7", "Batch9","Batch10", "Batch11")){
    Obj_proccesed <- ProcessSeu(SeuObjBatch)
  }
  else{
    SeuObjBatch <- CellCycleScoring(SeuObjBatch, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    Obj_proccesed <- ProcessSeu(SeuObjBatch)
  }
  SeuObjProcessedList[[batch_id]] <- Obj_proccesed
}



