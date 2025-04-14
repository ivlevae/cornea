library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratDisk)

table(prefinalObj_filtered_healthy$batch)

DefaultAssay(prefinalObj_filtered_healthy) <- 'RNA'

prefinalObj_sex_male_healthy <- subset(prefinalObj_filtered_healthy, (sex == c( 'M')))
prefinalObj_sex_female_healthy  <- subset(prefinalObj_filtered_healthy, (sex == c( 'F')))


table(cornea150_relabeled_filtered_30_2_30_1$detailed_annot)

cornea_subset_male_2<- subset(prefinalObj_sex_male_healthy, 
                idents = intersect(levels(prefinalObj_sex_male_healthy),
                                   c('Keratocytes', 'Myofibroblasts', 'Limbal Basal')))
                                 #c('Keratocytes', 'Myofibroblasts', 'Corneal Basal', 'Corneal Wing')))
cornea_subset_male_2$detailed_annot <- as.factor(cornea_subset_male_2$detailed_annot)
cornea_subset_male_2@meta.data$detailed_annot <- droplevels(cornea_subset_male_2@meta.data$detailed_annot)
cornea_subset_male_2 <- subset(cornea_subset_male_2, idents = unique(cornea_subset_male_2@meta.data$detailed_annot))
table(cornea_subset_male_2$detailed_annot)


cornea_subset_female_2<- subset(prefinalObj_sex_female_healthy, 
                              idents = intersect(levels(prefinalObj_sex_female_healthy),
                                                 c('Keratocytes', 'Myofibroblasts', 'Limbal Basal')))
                                                #c('Keratocytes', 'Myofibroblasts', 'Corneal Basal', 'Corneal Wing')))
cornea_subset_female_2$detailed_annot <- as.factor(cornea_subset_female_2$detailed_annot)
cornea_subset_female_2@meta.data$detailed_annot <- droplevels(cornea_subset_female_2@meta.data$detailed_annot)
cornea_subset_female_2 <- subset(cornea_subset_female_2, idents = unique(cornea_subset_female_2@meta.data$detailed_annot))

table(cornea_subset_female_2$detailed_annot)



###################################################################################

cellchat_male_2 <- createCellChat(object = cornea_subset_male_2, group.by = "detailed_annot")

CellChatDB <- CellChatDB.human
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- CellChatDB
cellchat_male_2@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_male_2 <- subsetData(cellchat_male_2) # This step is necessary even if using the whole database
cellchat_male_2 <- identifyOverExpressedGenes(cellchat_male_2)
cellchat_male_2 <- identifyOverExpressedInteractions(cellchat_male_2)
# project gene expression data onto PPI (Optional: when running it, USER should set 
# `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_male_2 <- smoothData(cellchat_male_2, adj = PPI.human)
cellchat_male_2 <- computeCommunProb(cellchat_male_2, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)


cellchat_male_2 <- filterCommunication(cellchat_male_2, min.cells = 10)
cellchat_male_2 <- computeCommunProbPathway(cellchat_male_2)
cellchat_male_2 <- aggregateNet(cellchat_male_2)
cellchat_male_2 <- netAnalysis_computeCentrality(cellchat_male_2, slot.name = 'netP')

saveRDS(cellchat_male_2, paste0('/home/bnvlab2/Documents/Kate/Cornea/cellchat_male_subset2', ', .rds'))

cellchat_male_2 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_male_subset2, .rds')
cellchat_male_1 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_male_subset1, .rds')
##############################################################################

cellchat_female_2 <- createCellChat(object = cornea_subset_female_2, group.by = "detailed_annot")

CellChatDB <- CellChatDB.human
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- CellChatDB
cellchat_female_2@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_female_2 <- subsetData(cellchat_female_2) # This step is necessary even if using the whole database
cellchat_female_2 <- identifyOverExpressedGenes(cellchat_female_2)
cellchat_female_2 <- identifyOverExpressedInteractions(cellchat_female_2)
# project gene expression data onto PPI (Optional: when running it, USER should set 
# `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_female_2 <- smoothData(cellchat_female_2, adj = PPI.human)
cellchat_female_2 <- computeCommunProb(cellchat_female_2, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)


cellchat_female_2 <- filterCommunication(cellchat_female_2, min.cells = 10)
cellchat_female_2 <- computeCommunProbPathway(cellchat_female_2)
cellchat_female_2 <- aggregateNet(cellchat_female_2)
cellchat_female_2 <- netAnalysis_computeCentrality(cellchat_female_2, slot.name = 'netP')

saveRDS(cellchat_female_2, paste0('/home/bnvlab2/Documents/Kate/Cornea/cellchat_female_subset2', ', .rds'))
cellchat_female_1 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_female_subset1, .rds')
cellchat_female_2 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_female_subset2, .rds')




##############################################################################
prefinalObj_filtered_kera <- subset(prefinalObj_filtered, condition_detailed == 'Keratoconus')
dim(prefinalObj_filtered_kera)
table(prefinalObj_filtered_kera$detailed_annot)
prefinalObj_filtered_kera_2<- subset(prefinalObj_filtered_kera, 
                                        idents = intersect(levels(prefinalObj_filtered_kera),
                                                           ### no limbal Suprabasal and limbal fibroblasts
                                                           c('Keratocytes', 'Myofibroblasts', 'Limbal Basal')))
                                                           #c('Keratocytes', 'Myofibroblasts', 'Corneal Basal', 'Corneal Wing')))
prefinalObj_filtered_kera_2$detailed_annot <- as.factor(prefinalObj_filtered_kera_2$detailed_annot)
prefinalObj_filtered_kera_2@meta.data$detailed_annot <- droplevels(prefinalObj_filtered_kera_2@meta.data$detailed_annot)
prefinalObj_filtered_kera_2 <- subset(prefinalObj_filtered_kera_2, idents = unique(prefinalObj_filtered_kera_2@meta.data$detailed_annot))

table(prefinalObj_filtered_kera_2$detailed_annot)


cellchat_kera_2 <- createCellChat(object = prefinalObj_filtered_kera_2, group.by = "detailed_annot")

CellChatDB <- CellChatDB.human
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- CellChatDB
cellchat_kera_2@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_kera_2 <- subsetData(cellchat_kera_2) # This step is necessary even if using the whole database
cellchat_kera_2 <- identifyOverExpressedGenes(cellchat_kera_2)
cellchat_kera_2 <- identifyOverExpressedInteractions(cellchat_kera_2)
# project gene expression data onto PPI (Optional: when running it, USER should set 
# `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_kera_2 <- smoothData(cellchat_kera_2, adj = PPI.human)
cellchat_kera_2 <- computeCommunProb(cellchat_kera_2, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)


cellchat_kera_2 <- filterCommunication(cellchat_kera_2, min.cells = 10)
cellchat_kera_2 <- computeCommunProbPathway(cellchat_kera_2)
cellchat_kera_2 <- aggregateNet(cellchat_kera_2)
cellchat_kera_2 <- netAnalysis_computeCentrality(cellchat_kera_2, slot.name = 'netP')

saveRDS(cellchat_kera_2, paste0('/home/bnvlab2/Documents/Kate/Cornea/cellchat_kera_subset2', ', .rds'))

cellchat_kera_2 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_kera_subset2, .rds')
cellchat_kera_1 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_kera_subset1, .rds')
##############################################################################

prefinalObj_filtered_healthy_2<- subset(prefinalObj_filtered_healthy, 
                                idents = intersect(levels(prefinalObj_filtered_healthy),
                                                   c('Keratocytes', 'Myofibroblasts', 'Limbal Basal')))
                                                  # c('Keratocytes', 'Myofibroblasts', 'Corneal Basal', 'Corneal Wing')))
prefinalObj_filtered_healthy_2$detailed_annot <- as.factor(prefinalObj_filtered_healthy_2$detailed_annot)
prefinalObj_filtered_healthy_2@meta.data$detailed_annot <- droplevels(prefinalObj_filtered_healthy_2@meta.data$detailed_annot)
prefinalObj_filtered_healthy_2 <- subset(prefinalObj_filtered_healthy_2, idents = unique(prefinalObj_filtered_healthy_2@meta.data$detailed_annot))


cellchat_healthy_2 <- createCellChat(object = prefinalObj_filtered_healthy_2, group.by = "detailed_annot")

CellChatDB <- CellChatDB.human
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- CellChatDB
cellchat_healthy_2@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_healthy_2 <- subsetData(cellchat_healthy_2) # This step is necessary even if using the whole database
cellchat_healthy_2 <- identifyOverExpressedGenes(cellchat_healthy_2)
cellchat_healthy_2 <- identifyOverExpressedInteractions(cellchat_healthy_2)
# project gene expression data onto PPI (Optional: when running it, USER should set 
# `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_healthy_2 <- smoothData(cellchat_healthy_2, adj = PPI.human)
cellchat_healthy_2 <- computeCommunProb(cellchat_healthy_2, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)


cellchat_healthy_2 <- filterCommunication(cellchat_healthy_2, min.cells = 10)
cellchat_healthy_2 <- computeCommunProbPathway(cellchat_healthy_2)
cellchat_healthy_2 <- aggregateNet(cellchat_healthy_2)
cellchat_healthy_2 <- netAnalysis_computeCentrality(cellchat_healthy_2, slot.name = 'netP')

saveRDS(cellchat_healthy_2, paste0('/home/bnvlab2/Documents/Kate/Cornea/cellchat_healthy_subset2', ', .rds'))

cellchat_healthy_2 <- readRDS('/home/bnvlab2/Documents/Kate/Cornea/cellchat_healthy_subset2, .rds')
###############################################################################################


object.list <- list(male = cellchat_male_2, female = cellchat_female_2)
MERGED_cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
gg1 <- rankNet(MERGED_cellchat, mode = "comparison", comparison = c(1:2), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat, mode = "comparison", comparison = c(1:2),stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Keratocytes") 
#xlims = c(-0.05, 0.05),ylims = c(-0.04, 0.04)) #, signaling.exclude = "MIF")
gg1

gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Limbal Basal")# ,
#xlims = c(-0.025, 0.05),ylims = c(-0.02, 0.04))
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Myofibroblasts")
#xlims = c(-0.005, 0.05),ylims = c(-0.005, 0.04))
gg1



###############################################################################################


object.list_kera_m1 <- list(male_kera = cellchat_kera_1, male_healthy = cellchat_male_1)
MERGED_cellchat_kera_m1 <- mergeCellChat(object.list_kera_m1, add.names = names(object.list_kera_m1), cell.prefix = T)
gg1 <- rankNet(MERGED_cellchat_kera_m, mode = "comparison", comparison = c(1:2), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat_kera_m, mode = "comparison", comparison = c(1:2),stacked = F, do.stat = TRUE)
gg1 + gg2


gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m1, idents.use = "Keratocytes")
gg1

gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Keratocytes",
xlims = c(-0.05, 0.05),ylims = c(-0.04, 0.04)) #, signaling.exclude = "MIF")
gg1

gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Limbal Basal")# ,
#xlims = c(-0.025, 0.05),ylims = c(-0.02, 0.04))
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Myofibroblasts")
#xlims = c(-0.005, 0.05),ylims = c(-0.005, 0.04))
gg1

unique_f_2 <- (cellchat_kera_2@netP$pathways[!cellchat_kera_2@netP$pathways %in% cellchat_male_2@netP$pathways])
unique_m_2 <- (cellchat_male_2@netP$pathways[!cellchat_male_2@netP$pathways %in% cellchat_kera_2@netP$pathways])
common_2 <- (cellchat_kera_2@netP$pathways[cellchat_kera_2@netP$pathways %in% cellchat_male_2@netP$pathways])
length(unique_f_2)
length(unique_m_2)
length(common_2)

unique_f_1 <- (cellchat_kera_1@netP$pathways[!cellchat_kera_1@netP$pathways %in% cellchat_male_1@netP$pathways])
unique_m_1 <- (cellchat_male_1@netP$pathways[!cellchat_male_1@netP$pathways %in% cellchat_kera_1@netP$pathways])
common_1 <- (cellchat_kera_1@netP$pathways[cellchat_kera_1@netP$pathways %in% cellchat_male_1@netP$pathways])
length(unique_f_2)
length(unique_m_2)
length(common_1)


length(unique_m_2[unique_m_2 %in% unique_m_1])
length(unique_m_2[!unique_m_2 %in% unique_m_1])
length(unique_m_1[!unique_m_1 %in% unique_m_2])

length(cellchat_female_2@netP$pathways)
length(cellchat_male_2@netP$pathways)

length(unique_f_1[unique_f_2 %in% unique_f_1])
length(unique_f_1[!unique_f_1 %in% unique_f_2])
length(unique_f_2[!unique_f_2 %in% unique_f_1])





unique_f_2 <- (cellchat_female_2@netP$pathways[!cellchat_female_2@netP$pathways %in% cellchat_male_2@netP$pathways])
unique_m_2 <- (cellchat_male_2@netP$pathways[!cellchat_male_2@netP$pathways %in% cellchat_female_2@netP$pathways])
common_2 <- (cellchat_female_2@netP$pathways[cellchat_female_2@netP$pathways %in% cellchat_male_2@netP$pathways])

unique_f_1 <- (cellchat_female_1@netP$pathways[!cellchat_female_1@netP$pathways %in% cellchat_male_1@netP$pathways])
unique_m_1 <- (cellchat_male_1@netP$pathways[!cellchat_male_1@netP$pathways %in% cellchat_female_1@netP$pathways])
common_1 <- (cellchat_female_1@netP$pathways[cellchat_female_1@netP$pathways %in% cellchat_male_1@netP$pathways])

unique_m_2[unique_m_2 %in% unique_m_1]
length(unique_m_2[!unique_m_2 %in% unique_m_1])
length(unique_m_1[!unique_m_1 %in% unique_m_2])

length(cellchat_female_2@netP$pathways)
length(cellchat_male_2@netP$pathways)


object.list <- list(male = cellchat_male_2, female = cellchat_female_2)
MERGED_cellchat_1 <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
gg1 <- rankNet(MERGED_cellchat_1, mode = "comparison", comparison = c(1:2), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat_1, mode = "comparison", comparison = c(1:2),stacked = F, do.stat = TRUE)

gg1$data <- gg1$data %>%
  group_by(name) %>%
  mutate(contribution_pct = (contribution / sum(contribution)) * 100)
View(gg1$data)
filtered_data <- gg1$data %>%
  filter(group == "male" & (contribution_pct < 40 | contribution_pct > 60)  & pvalues < 0.05)
dim(filtered_data)

filtered_data1 <- gg1$data %>%
  filter(group == "male"  & pvalues < 0.05)
dim(filtered_data1)

gene_2 <- filtered_data$name
gene_21 <- filtered_data1$name

gene_1 <- filtered_data$name
gene_11 <- filtered_data1$name

length(gene_11[!gene_11 %in% gene_21])
length(gene_21[!gene_21 %in% gene_11])
length(gene_21[gene_21 %in% gene_11])

netAnalcellchat_male_2netAnalysis_signalingRole_heatmap(cellchat_female, pattern = "incoming", 
                                  signaling = pathways_of_interest_2,  
                                  width = 6, height = 7, font.size = 10)

netAnalysis_signalingRole_heatmap(cellchat1, pattern = "outgoing", 
                                  signaling = pathways_of_interest_2,  
                                  width = 6, height = 7, font.size = 10)



cellchat_female@netP$pathways[!cellchat_female@netP$pathways %in% cellchat_male@netP$pathways]
cellchat_male@netP$pathways[!cellchat_male@netP$pathways %in% cellchat_female@netP$pathways]
cellchat_male@netP$pathways[cellchat_male@netP$pathways %in% cellchat_female@netP$pathways]
View(cellchat_male@netP$pathways)

table(cellchat_kera@meta$ident)
table(cellchat_healthy@meta$ident)

object.list1 <- list(kera = cellchat_kera, healthy = cellchat_healthy)
MERGED_cellchat1 <- mergeCellChat(object.list1, add.names = names(object.list1), cell.prefix = T)
gg1 <- rankNet(MERGED_cellchat1, mode = "comparison", comparison = c(1:2), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat1, mode = "comparison", comparison = c(1:2),stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Keratocytes",
                                            xlims = c(-0.05, 0.05), ylims = c(-0.04, 0.04)) #, signaling.exclude = "MIF")
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Corneal Wing")#,
                                            xlims = c(-0.05, 0.05),ylims = c(-0.04, 0.04))
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Corneal Basal")# ,
                                           # xlims = c(-0.025, 0.05),ylims = c(-0.02, 0.04))
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat_kera_m, idents.use = "Myofibroblasts")
                                            #xlims = c(-0.005, 0.05),ylims = c(-0.005, 0.04))
gg1






object.list_kera_m <- list(male_kera = cellchat_kera_2, male_healthy = cellchat_male_2)
MERGED_cellchat_kera_m <- mergeCellChat(object.list_kera_m, add.names = names(object.list_kera_m), cell.prefix = T)
table(cellchat_male@meta$ident)
table(cellchat_kera@meta$ident)



gg1 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:2))
gg2 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:2), measure = "weight")
gg1 + gg2
gg1 <- rankNet(MERGED_cellchat_kera_m, mode = "comparison", comparison = c(1:2), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat_kera_m, mode = "comparison", comparison = c(1:2),stacked = F, do.stat = TRUE)
gg1 + gg2
View(gg1)
gg1$data$name
gg <- ggplot(gg1$data, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity", position ="fill")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg
gg1$data <- gg1$data %>%
  group_by(name) %>%
  mutate(contribution_pct = (contribution / sum(contribution)) * 100)
View(gg1$data)
filtered_data <- gg1$data %>%
  filter(group == "male_kera" & pvalues < 0.05)
dim(filtered_data)
pathways_of_interest_2 <- filtered_data$name
pathways_of_interest_2
filtered_data <- gg1$data %>%
  filter(group == "male_kera" & (contribution_pct < 40 | contribution_pct > 60)  & pvalues < 0.05)
dim(filtered_data)
pathways_of_interest_2 <- filtered_data$name
pathways_of_interest_2
# View the filtered data
View(filtered_data)


list1 <- c("PARs", "IL6", "CD46", "TENASCIN", "VTN", "MHC-I", "CD6", "ICAM", "CALCR", "UNC5", "OSM", "CCK", 
           "NRXN", "RANKL", "FLRT", "ACTIVIN", "CXCL", "ANGPT", "MSTN", "LT", "SIRP", "TWEAK", "GDF", "HSPG", 
           "KLK", "VCAM", "LIGHT", "CD137", "THY1", "Desmosterol", "CD23", "NGF", "BTLA", "VEGI", "CD80", "EPO", 
           "CD86", "AGRN", "GDNF", "ApoA", "TAC", "IL2", "BAFF", "NPR2", "GP1BA", "NEGR", "COMPLEMENT", "CSF3", 
           "FLT3", "CCL", "IL1", "CNTN", "PTN", "IGF", "MPZ", "LIFR", "CSF", "BMP", "EPHB", "NOTCH", "WNT", "RECN", 
           "SEMA4", "DHT", "FGF", "EPHA", "PTPR", "DHEA", "APP", "TGFb", "PDGF", "EGF", "NT", "CADM", "HGF", 
           "SEMA3", "NRG", "CDH1", "LAMININ", "NCAM", "CDH5", "PCDH", "IL16", "SLURP", "SEMA5", "SLITRK")

list2 <- c("TENASCIN", "PARs", "MHC-I", "VTN", "IL6", "CD46", "IL10", "CCL", "ANGPT", "ICAM", "RANKL", "UNC5", 
           "CCK", "ACTIVIN", "CXCL", "ESAM", "OSM", "SIRP", "FLRT", "CD6", "CALCR", "MSTN", "TWEAK", "LT", 
           "Glutamate", "NRXN", "FASLG", "Desmosterol", "GDF", "LIGHT", "KLK", "NMU", "L1CAM", "HSPG", "CD23", 
           "NGF", "GDNF", "BTLA", "VEGI", "THY1", "CD80", "EPO", "CD137", "CD86", "CSF3", "GH", "BAFF", "ApoA", 
           "IL2", "AGRN", "NPR2", "GP1BA", "TAC", "CEACAM", "FLT3", "COMPLEMENT", "IL1", "MPZ", "CNTN", "LIFR", 
           "IGF", "EPHB", "BMP", "NOTCH", "RELN", "PLAU", "FGF", "CSF", "DHT", "WNT", "TGFb", "PDGF", "APP", 
           "EPHA", "DHEA", "ANGPTL", "NT", "MIF", "CDH1", "EGF", "LAMININ", "PTPRM", "PCDH", "IL16", "SLURP", 
           "SEMA5", "SLITRK", "ADGRA")

list2[list2 %in% list1]
# Find genes in List 1 that are not in List 2
unique_to_list1 <- setdiff(list1, list2)

# Find genes in List 2 that are not in List 1
unique_to_list2 <- setdiff(list2, list1)

# Output the results
length(unique_to_list1)  # Number of genes in List 1 but not in List 2
length(unique_to_list2)  # Number of genes in List 2 but not in List 1

??netAnalysis_signalingChanges_scatter
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Keratocytes", xlims = c(-0.025, 0.015), , ylims = c(-0.025, 0.005)) 
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Corneal Wing", xlims = c(-0.025, 0.015), , ylims = c(-0.025, 0.005))
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Corneal Basal", xlims = c(-0.025, 0.015), , ylims = c(-0.025, 0.005))
gg1
gg1 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat, idents.use = "Myofibroblasts", xlims = c(-0.05, 0.05), , ylims = c(-0.05, 0.025))
gg1

pathways_of_interest_2 <- filtered_data$name
pathways_of_interest_2
netAnalysis_signalingRole_heatmap(cellchat_female, pattern = "incoming", 
                                  signaling = pathways_of_interest_2, 
                                  width = 6, height = 11, font.size = 8)

netAnalysis_signalingRole_heatmap(cellchat_female, pattern = "outgoing", 
                                  signaling = pathways_of_interest_2,  
                                  width = 6, height = 11, font.size = 8)

netAnalysis_signalingRole_heatmap(cellchat_male, pattern = "incoming", 
                                  signaling = pathways_of_interest_2,  
                                  width = 6, height = 11, font.size = 8)

netAnalysis_signalingRole_heatmap(cellchat_male, pattern = "outgoing", 
                                  signaling = pathways_of_interest_2,  
                                  width = 6, height = 11, font.size = 8)



table(cellchat_male@meta$ident)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))


gg11 <- netAnalysis_signalingChanges_scatter(MERGED_cellchat1, idents.use = "Keratocytes") #, signaling.exclude = "MIF")
gg11



gg1 <- rankNet(MERGED_cellchat, mode = "comparison", signaling = 'TGFb', comparison = c(1:2), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat, mode = "comparison",  signaling = 'TGFb',comparison = c(1:2),stacked = F, do.stat = TRUE)
gg1 + gg2
netVisual_bubble(MERGED_cellchat,  sources.use = 1, targets.use = c(1,2),  comparison = c(1:2), angle.x = 45)
netVisual_bubble(MERGED_cellchat,  sources.use = 2, targets.use = c(1,2),  comparison = c(1:2), angle.x = 45)

all_pathways <- vector("list", length(object.list))
all_pathways
# Loop through each group to collect all unique pathways
for (i in 1:length(object.list)) {
  all_pathways[[i]] <- object.list[[i]]@netP$pathways
}
# Reduce the list of pathways into a single vector containing the union of all pathways
pathway.union_all <- Reduce(union, all_pathways)
pathway.union_all
'EGFR' %in% pathway.union_all
pathway.union <- Reduce(union, selected_pathways)
















table(cornea150_relabeled_filtered_30_2_30_1@active.ident)

length(cellchat1@netP$pathways)
View(cellchat1@netP$pathways)

pathways_of_interest <- c('COLLAGEN', 'LAMININ', 'FN1','VTN','TENASCIN','HSPG','PLAU')
plot_list <- list()
for (pathway in pathways_of_interest) {
  pathways.show <- c(pathway) 
  vertex.receiver <- seq(1, 4) # Ensure this is correct for your specific data
  
  # Generate and store the network plot
  net_plot <- netVisual_aggregate(cellchat1, signaling = pathways.show, vertex.receiver = vertex.receiver)
  #plot_list[[pathway]] <- net_plot
  
  # Generate and store the circular layout plot
  par(mfrow = c(1, 1)) # Ensure this is correct for your plotting environment
  circle_plot <- netVisual_aggregate(cellchat1, signaling = pathways.show, layout = "circle")
  plot_list[[paste0(pathway, "_circle")]] <- circle_plot
}
plot_list[[1]] 
plot_list[[2]] 
plot_list[[3]] 
plot_list[[4]] 
plot_list[[5]] 
plot_list[[6]] 
plot_list[[7]] 


cord_list <- list()
for (pathway in pathways_of_interest) {
  pathways.show <- c(pathway) 
  strwidth <- function(str) {0.3}
  par(mfrow=c(1,1))
  ht <- netVisual_aggregate(cellchat1, signaling = pathways.show, layout = "chord", pt.title = 3) #,
  #title = names(object.list)[i] ,  pt.title = 10)
  cord_list[[pathway]] <- ht
}
cord_list[[1]]


View(cellchat1)
pathways.show <- c("PLAU") 
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot for NRXN pathway
par(mfrow=c(1,1))
netVisual_aggregate(cellchat1, signaling = pathways.show, layout = "circle")


par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat1, signaling = pathways.show, layout = "chord")

length(cellchat@meta$cell_names) 
cellchat@meta$cell_names <- as.factor(cellchat@meta$cluster_annot)





save_recordedplot_pdf <- function(plot, filename, width = 10, height = 6) {
  # Use grDevices::pdf() to avoid namespace issues
  grDevices::pdf(file = filename, width = width, height = height)
  replayPlot(plot)  # Replay the recorded plot
  dev.off()         # Close the graphics device
}

# Ensure no other devices are open
while (dev.cur() > 1) dev.off()

# Iterate over plot_list and save each plot as a PDF with wider dimensions
for (name in names(cord_list)) {
  save_recordedplot_pdf(cord_list[[name]], paste0(name, "_chord.pdf"), width = 8, height = 6)  # Adjust width and height as needed
}


pathways_of_interest_2 <- list(
  "FGF", "EGF", "TGFb", "NOTCH", ##Differentiation
"PLAU",  # Example MMPs

  "COLLAGEN", "LAMININ", "FN1", "TENASCIN", "DESMOSOME", ## ECM
  "VEGF", "CDH" # Other relevant factors
)


netAnalysis_signalingRole_heatmap(cellchat1, pattern = "incoming", 
                                        signaling = pathways_of_interest_2,  
                                  width = 6, height = 7, font.size = 10)

netAnalysis_signalingRole_heatmap(cellchat1, pattern = "outgoing", 
                                  signaling = pathways_of_interest_2,  
                                  width = 6, height = 7, font.size = 10)




groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = c('COLLAGEN'), color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = c('LAMININ'), color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = c('VTN'), color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = c('TENASCIN'), color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = c('HSPG'), color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = c('FN1'), color.heatmap = "Reds")
pathways_of_interest
heatmap_list <- list()
for (pathway in pathways_of_interest) {
  pathways.show <- c(pathway) 
  par(mfrow=c(1,1))
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  #title = names(object.list)[i] ,  pt.title = 10)
  #heatmap_list[[pathway]] <- ht
}

heatmap_list[[0]]





cornea150_no_conj@active.ident  <- cornea150_no_conj$cell_types_main

cornea150_no_4types <- subset(cornea150_no_conj, cell_types_main != 'non-Myelinating Shwann cells')
cornea150_no_4types <- subset(cornea150_no_4types, cell_types_main != 'melanocytes')
cornea150_no_4types <- subset(cornea150_no_4types, cell_types_main != 'vessel')
cornea150_no_conj <- cornea150_no_4types
# cornea150_no_4types <- subset(cornea150_no_conj, !cell_types_main %in% 
#                               c('non-Myelinating Shwann cells', 'melanocytes', 'vessel'))

cornea150_no_conj@active.ident <- droplevels(cornea150_no_conj@active.ident, exclude = setdiff(levels(cornea150_no_conj@active.ident), unique(cornea150_no_conj@active.ident)))
cornea150_no_conj$cell_types_main <- droplevels(cornea150_no_conj$cell_types_main, exclude = setdiff(levels(cornea150_no_conj$cell_types_main), unique(cornea150_no_conj$cell_types_main)))

table(cornea150_no_conj$cell_types_main)

table(cornea150_no_conj$condition_detailed)

condition_sub <- 'healthy'
cornea150_no_conj_sub <- subset(cornea150_no_conj, condition_detailed == condition_sub)

if (identical(cornea30_noconj_V2_res1_annot_act_endo, cornea30_noconj_V2_res1_annot_act_endo_2)) {
  print("The data frames are identical.")
} else {
  print("The data frames are not identical.")
}


cellchat <- createCellChat(object = cornea150_no_conj_sub, group.by = "cell_types_main")

CellChatDB <- CellChatDB.human
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use




# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set 
# `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- smoothData(cellchat, adj = PPI.human)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')

saveRDS(cellchat, paste0('/home/bnvlab2/Documents/Kate/cellchat_no4types', condition_sub ,', .rds'))

#droplevels(subset_df)

########################## Visualization ##########################
table(cornea150_no_conj$condition_detailed)

condition_sub <- 'healthy'

cellchat <- readRDS(paste0('/home/bnvlab2/Documents/Kate/cellchat_no4types', condition_sub ,', .rds'))

list_pathways_healthy  <- cellchat@netP$pathways   
length(list_pathways_healthy)
list_pathways_kera <- cellchat_kera@netP$pathways
length(list_pathways_kera)

length(setdiff(list_pathways_kera, list_pathways_healthy))
length(setdiff(list_pathways_healthy,list_pathways_kera))
length(intersect(list_pathways_healthy,list_pathways_kera))

'EIF2' %in% intersect(list_pathways_healthy,list_pathways_kera)


pathways.show <- c("COLLAGEN") 
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat_kera, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot for NRXN pathway
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_kera, signaling = pathways.show, layout = "circle")

pathways.show <- c("COLLAGEN") 
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot for NRXN pathway
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle") #, title.space = 1)

?netVisual_aggregate



table(cellchat@meta$cell_types_main)

cellchat_kera@meta$condition
object.list <- list(healthy = cellchat, keratoconus = cellchat_kera,
                    dysplasia = cellchat_dysp, cataract = cellchat_kera)

MERGED_cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
gg1 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:4))
gg2 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:4), measure = "weight")
gg1 + gg2
gg1 <- rankNet(MERGED_cellchat, mode = "comparison", comparison = c(1:4), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat, mode = "comparison", comparison = c(1:4),stacked = F, do.stat = TRUE)
gg1 + gg2
gg1 <- rankNet(MERGED_cellchat, mode = "comparison", signaling = 'COLLAGEN', comparison = c(1:4), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat, mode = "comparison",  signaling = 'COLLAGEN',comparison = c(1:4),stacked = F, do.stat = TRUE)
gg1 + gg2



# Calculate average expression for each gene by cell type
mmp_genes <- c(
  'COL21A1',
  'COL6A2',
  'COL5A2',
  'COL1A2',
  'COL12A1',
  'COL6A3',
  'COL24A1',
  'COL6A1',
  'COL8A2',
  'COL23A1'
  
  #'TGFB1', 'TGFBR1', 'TGFBR2', 'TGFB3'
  # "MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11",
  # "MMP12", "MMP13", "MMP14", "MMP15", "MMP16", "MMP17", "MMP19",
  # "MMP20", "MMP21", "MMP23", "MMP24", "MMP25", "MMP26", "MMP27",
  # "MMP28", "MMP29", "MMP30", "MMP31", "MMP32"
)
cornea150_no_conj@meta.data$condition

DotPlot(cornea150_no_conj, assay = "RNA", features = mmp_genes,  
        split.by = 'condition',cols = c("lightblue", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))




library(dplyr)
library(tidyr)
library(ggplot2)
cornea150_no_conj
pbmc_small <- ScaleData(cornea150_no_conj, rownames(cornea150_no_conj))
?DotPlot
pbmc_small

DotPlot(cornea150_no_conj, assay = "RNA", features = mmp_genes, 
        cols = c("lightblue","red"),  split.by = 'condition') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) 

library(Seurat)
library(ggplot2)

# Generate the dot plot with custom color scale
dot_plot <- DotPlot(cornea150_no_conj, assay = "RNA", features = mmp_genes, 
                    cols = c("lightblue", "red"),  # Custom color scale
                    split.by = 'condition') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Rotate x-axis labels
    legend.position = "right"  # Position legend
  ) +
  scale_color_manual(values = c("lightblue", "red"),  # Custom color scale
                     name = "Expression",  # Legend title
                     breaks = c("lightblue", "red"),
                     labels = c("Low Expression", "High Expression")) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Customize legend

# Display the plot
print(dot_plot)

DotPlot(cornea150_no_conj, assay = "RNA", features = mmp_genes, 
        split.by = 'condition') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Rotate x-axis labels
    legend.position = "right"  # Position legend
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  
?DotPlot

DotPlot(cornea150_no_conj, assay = "RNA", features = mmp_genes, 
        split.by = 'condition',cols = c("lightblue", "red")) +
 # Use discrete color scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))





seurat <- cornea150_no_conj
#Idents(object = seurat) <- new.idents #add cell type as idents
seurat$celltype_new<- paste(Idents(seurat),seurat$condition, sep = "_") #add sample subgruops to the cell type idents
new.idents <- seurat@meta.data$celltype_new
Idents(object = seurat) <- new.idents # now just add it to the seurat object as new idents

DotPlot(seurat, dot.scale = 8, features = mmp_genes) + RotatedAxis() + theme(axis.text = element_text(size = 14))  +
scale_colour_gradient(low =c("lightblue"), high =c("red"))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12)  # Remove color legend
  ) 
#+
  #+guides(color = guide_colorbar(title = 'Average Expression'))
# 
# DotPlot(cornea150_no_conj, dot.scale = 8, features = mmp_genes) + 
#   RotatedAxis() + theme(axis.text = element_text(size = 14))  +
#   scale_colour_gradient(low =c("lightblue"), high =c("red"))+
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, size = 12)  # Remove color legend
#   )
