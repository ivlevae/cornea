################################### bimarkers plot##########################################


desired_order <- c("Corneal Superficial",
                   "Corneal Wing",
                   "Corneal Basal",
                   "TAC",
                   'Limbal Suprabasal', 
                   "Limbal Basal",
                   "Keratocytes",
                   "Myofibroblasts",
                   'Limbal fibroblasts',
                   "Endothelium")

table(prefinalObj$detailed_annot)

prefinalObj$detailed_annot <- factor(prefinalObj$detailed_annot, levels = desired_order)
prefinalObj@active.ident <- prefinalObj$detailed_annot
table(prefinalObj$detailed_annot)
rownames(prefinalObj)

markers_subset_2 <- c( 
  # "PAX6", "KRT12", "TACSTD2", ##corneal epithelium
  "MUC16", "KRT24",  "MACC1",  "HOPX", # "NECTIN4", #"WFDC21P", "BCAS1", # 'Corneal Superficial'
  "HES5", "KRT3", "DIO2",  "KC6", #"FOXP2", "GALNT18", "CACNA1E", ##Corneal Wing
  "NKAIN2", "TENM2","LAMA3",  "CDH13", #"IVNS1ABP", "MIR205HG", #'Corneal Basal' 
  "ATAD2",   "TOP2A", "DIAPH3",  "ANLN",#,  'S100A2', #"ATAD2",   "TOP2A", "MKI67", "RIMS2", "POLQ", "MIR924HG","BRIP1",   #TAC
  'CXCL14',  "CSRP2",   "KRT15",  "KRT14", #suprabasal #"CLDN4","LGR6", 
  "GJA1", ## limbal basal
  "PLAUR", 
  "DCN", "LUM",   "KERA", "ABCA6", #, "ITGBL1" , "COL6A3",  ##keratocytes
  'RGS5',  'ID4', 'NR2F2',  'MYL9', # 'TPM2', 'TPM1','MYH11', 'ACTA2', 'MYLK' ,'SPARCL1', 'TAGLN', ### myofibroblasts
  "FBLN1", 'COL3A1',  "COL1A1", #'SOCS3','SRPX', # 'COL1A2', 'SFRP2',  'ELN','MMP2', ##limbal fibroblasts 
  "POU6F2",   "SLC4A4", "COL4A3"#, "NRXN3", "COL8A1", "CA3","CA12", "FAM155A", "COL4A4"  ##corneal endo
)             


DefaultAssay(prefinalObj) <- 'RNA'

p <- DotPlot(prefinalObj,
             features = rev(markers_subset_2) , #markers_subset_2,
             assay = NULL, 
             cols = c("lightgrey", "blue"))
p + theme(axis.text.x = element_text(angle = 90) )