################################### UMAP Figure B  #######################################################
prefinalObj@active.ident <- prefinalObj$detailed_annot

DimPlot(prefinalObj, reduction = "umap", raster = TRUE, label = T, label.box = T,
        repel = T) + NoLegend() 


################################### biomarkers plot Figure C  ##########################################


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


################################### UMAP by condition Figure C  #######################################################


DimPlot(prefinalObj, reduction = "umap", raster = FALSE, split.by = 'condition_detailed') + NoLegend()


################################### metadata line plot of cell percentage by age   Figure D  #######################################################


metadata <- prefinalObj@meta.data


# Count the number of each cell type per age group
age_group_counts <- metadata %>%
  group_by(age_num, detailed_annot) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate the total number of cells per age group
age_group_totals <- metadata %>%
  group_by(age_num) %>%
  summarise(total_cells = n())

# Calculate the percentage for each cell type in each age group
cell_type_percentages <- age_group_counts %>%
  left_join(age_group_totals, by = "age_num") %>%
  mutate(percentage = (count / total_cells) * 100)

cell_type_percentages_1 <- cell_type_percentages %>%
  group_by(detailed_annot) %>%
  summarise(Av_perc = mean(percentage)) %>%
  arrange(Av_perc)

order_for_plot <- cell_type_percentages_1$detailed_annot

cell_type_percentages_2 <- cell_type_percentages %>%
  # Group by age number and detailed annotation
  group_by(age_num) %>%
  arrange(age_num, factor(detailed_annot, levels = order_for_plot)) %>%
  # Summarize the cumulative sum of percentages
  mutate(Cum_perc = cumsum(percentage)) %>%
  # Reorder detailed_annot as specified in order_for_plot
  ungroup()



cell_type_percentages_2$age_character <- as.character(cell_type_percentages_2$age_num)




ggplot(cell_type_percentages_2, aes(x = age_character, y = Cum_perc, fill = detailed_annot, group = detailed_annot)) +
  geom_area(alpha = 0.6, position = "identity") +  # Fill the area under the lines with colors
  geom_line(color = "black", size = 0.8) +  # Draw the lines as black
  labs(
    x = "Age, years", 
    y = "Cell Percentages per Age, (%)"
  ) +
  scale_x_discrete(breaks = unique(cell_type_percentages_2$age_character)) +  # Add ticks for each age_num
  theme_minimal(base_size = 15) +  # Set base text size for better readability
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis labels
    axis.line = element_line(color = "black"),  # Add axis lines
    legend.title = element_blank()  # Remove legend title
  )



################################### nCells by sex and age  Figure E  #######################################################

table(prefinalObj$sex) 

prefinalObj_sex <- subset(prefinalObj, sex %in% c('F', 'M'))

meta <- prefinalObj_sex@meta.data


meta$sex <- as.factor(meta$sex)

ggplot(meta, aes(x = age_num, fill = sex)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("M" = "blue", "F" = "pink")) +
  labs(title = "Age Density Plot by Sex", x = "Age", y = "Density") +
  theme_minimal()

# to calculate  n of batches
# count_data <- meta %>%
#   group_by(age_num, sex) %>%
#   summarise(count = n_distinct(batch), .groups = "drop")


count_data <- meta %>%
  group_by(age_num, sex) %>%
  summarise(count = n(), .groups = "drop")



ggplot(count_data, aes(x = age_num, y = count, color = sex)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = 0, ymax = count, fill = sex), alpha = 0.3) + # Fill under the line
  labs(title = "Number of Samples by Age and Sex", x = "Age", y = "Number of Samples") +
  scale_x_continuous(
    breaks = seq(min(count_data$age_num), max(count_data$age_num), by = 4), # Start from min value with interval of 6
    expand = c(0, 0) # Remove extra space around the axis
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),   
    axis.ticks.x = element_line(linewidth = 0.5), 
    axis.ticks.length = unit(5, "pt") 
  )