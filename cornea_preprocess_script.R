library(Seurat)
library(DoubletFinder)
library(SeuratDisk)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters)
  nExp_poi <- round(0.05*length(colnames(tmp)))  ## Assuming 5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp)
}

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

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


###### Check that sample_id is changed!!!!!! ########
sample_id <- 'GSM5651520'
data_dir <- paste0('/home/bnvlab2/Documents/Kate/Cornea/', sample_id, '/')
list.files(data_dir)
cornea <- Read10X(data.dir = data_dir)


cornea_obj = CreateSeuratObject(counts = cornea)

cornea_obj[["percent.rb"]] <- PercentageFeatureSet(cornea_obj, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA') 
cornea_obj[["percent.mt"]] <- PercentageFeatureSet(cornea_obj, pattern = "^MT-") 
cornea_obj <- CellCycleScoring(cornea_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 12) 
VlnPlot(cornea_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
cornea_obj <- subset(cornea_obj, subset =  nCount_RNA < 100000)

cornea_obj <- subset(cornea_obj, subset = nCount_RNA > 300 & nCount_RNA < 27000 & nFeature_RNA > 400 & nFeature_RNA < 4000
                     & percent.mt < 30 & percent.rb < 40)




cornea_obj <- ProcessSeu(cornea_obj)
cornea_obj <- RDoublet(cornea_obj)
dim(cornea_obj)

DimPlot(cornea_obj)

cornea_obj <- subset(cornea_obj, cells = colnames(cornea_obj )[which(cornea_obj [[]][12] == 'Singlet')])
cornea_obj <- subset(cornea_obj , cells = colnames(cornea_obj )[which(cornea_obj [[]][13] == 'Singlet')])

cornea_obj <- ProcessSeu(cornea_obj)

DimPlot(cornea_obj)
cornea_obj$GSM <- sample_id 
dim(cornea_obj)

SaveH5Seurat(cornea_obj, paste0('/home/bnvlab2/Documents/Kate/Cornea/Output/GSE147979', sample_id, '.h5Seurat'), overwrite = TRUE)
saveRDS(cornea_obj, paste0('/home/bnvlab2/Documents/Kate/Cornea/Output/GSE147979',  sample_id, '.rds'))







# library(Seurat)
# library(dplyr)
# library(Matrix)
# library(data.table)
# 
# # "GSE155683_counts_adult_cornea.txt.gz"
# 
# file_path <- paste0("/home/bnvlab2/Documents/Kate/Cornea/","GSM4451263_CountMatrix_Cornea76.csv.gz")
# 
# # my_df <- read(gzfile(file_path))
# ###### Check that sample_id is changed!!!!!! ########
# sample_id <- 'GSM4451263'
# 
# 
# my_df1 <- gzfile(file_path)
# #file_path_csv <- paste0("/home/bnvlab2/Documents/Kate/Cornea/","GSM5360358_CountMatrix_Cornea73.csv.gz")
# # write.csv(gzfile(file_path), file = file_path_csv, row.names=T)
# count_matrix2 <- fread(file_path)  
# count_matrix2 [1:10, 1:5]
# rownames(count_matrix2) <- count_matrix2$V1
# names(count_matrix2)[1]
# count_matrix2$V1<- NULL
# rownames(count_matrix2)
# colnames(count_matrix2)[1:5]
# mydata1 <- CreateSeuratObject(counts = count_matrix2)





# GSM_ids <- c("GSM5962410", "GSM5962413", "GSM5962416", "GSM5962423", "GSM5962430", "GSM5962435")
# sample_names <- c("Hu216Cornea", "Hu216ODCornea", "Hu218Cornea", "Hu220Co", "Hu235Cornea", "Pt14Cornea")
# 
# # Create DataFrame
# df <- data.frame(GSMID = GSM_ids, SampleName = sample_names)
# subset(df, sample_names == 'Hu216Cornea')$GSMID 
# 
# 
# levels(cornea_obj$orig.ident)
# length(levels(cornea_obj$orig.ident))
# 
# for (patient_id in levels(cornea_obj$orig.ident)){
#   print(patient_id )
#   GSM_id <- subset(df, SampleName == patient_id)$GSMID 
#   cat('GSM number is:', GSM_id, '\n')
# }
# 
# 
# head(cornea_obj_subset)
# cornea_obj_subset <- subset(cornea_obj, subset = orig.ident == "Hu216Cornea")
# dim(cornea_obj)
# 
# for (patient_id in levels(cornea_obj$orig.ident)){
#   
#   print(patient_id)
#   cornea_obj_subset <- subset(cornea_obj, subset = orig.ident == patient_id)
#   cornea_obj_subset  <- ProcessSeu(cornea_obj_subset )
#   cornea_obj_subset  <- RDoublet(cornea_obj_subset )
#   dim(cornea_obj)
#   
#   DimPlot(cornea_obj_subset )
#   #FeaturePlot(cornea_obj, features = c('CD74','NEFL','EMCN','APOE'))
#   cornea_obj_subset  <- subset(cornea_obj_subset , cells = colnames(cornea_obj_subset )[which(cornea_obj_subset  [[]][12] == 'Singlet')])
#   cornea_obj_subset  <- subset(cornea_obj_subset  , cells = colnames(cornea_obj_subset )[which(cornea_obj_subset [[]][13] == 'Singlet')])
#   
#   cornea_obj_subset <- ProcessSeu(cornea_obj_subset )
#   
#   DimPlot(cornea_obj_subset )
#   GSM_id <- subset(df, SampleName == patient_id)$GSMID 
#   cornea_obj_subset$GSE <- 'GSE199013' 
#   cornea_obj_subset$GSM <- GSM_id
#   dim(cornea_obj_subset)
#   
#   cat('GSM number is:', GSM_id, '\n')
#   
#   SaveH5Seurat(cornea_obj_subset, paste0('/home/bnvlab2/Documents/Kate/Cornea/Output/GSE199013_', GSM_id, '.h5Seurat'), overwrite = TRUE)
#   saveRDS(cornea_obj_subset, paste0('/home/bnvlab2/Documents/Kate/Cornea/Output/GSE199013_',  GSM_id, '.rds'))
# }        