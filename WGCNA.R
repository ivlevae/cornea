library(WGCNA)
library(hdWGCNA)
DefaultAssay(prefinalObj) <- 'RNA'


################# the scrip was run for all the cells altogether and for keratocytes only ###################
prefinalObj_kera <- subset(prefinalObj, detailed_annot == 'Keratocytes')
prefinalObj_kera$detailed_annot <- droplevels(as.factor(prefinalObj_kera$detailed_annot))

table(prefinalObj_kera$detailed_annot)

seurat_obj <- SetupForWGCNA(
  prefinalObj_kera,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "test" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("condition_detailed"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'condition_detailed' # set the Idents of the metacell seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj@assays[["RNA"]]@var.features <- seurat_obj@assays[["integrated"]]@var.features
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='condition_detailed')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)


p1 <- DimPlotMetacells(seurat_obj, group.by='condition_detailed') + umap_theme() + ggtitle("condition_detailed")
p1

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Keratoconus", # the name of the group of interest in the group.by column
  group.by='condition_detailed', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'Keratoconus2' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='HUT hdWGCNA Dendrogram')

TOM <- GetTOM(seurat_obj)

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="condition_detailed"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'condition_detailed', group_name = 'Keratoconus'
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
library(patchwork) 
# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$condition_detailed), ' '))[,1]

ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster', axis.label.size=4,
  grid.label.size=4
)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'condition_detailed')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

ModuleNetworkPlot(
  seurat_obj,
  outdir = '/home/bnvlab2/Downloads/HUT_Networks'
)

g1<- HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'yellow', return_graph = TRUE
)

g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=10,
  edge_prop = 0.75,
  mods = c('brown','blue','turquoise', 'yellow') #,'red','green','black')
)


saveRDS(seurat_obj, 'testrun_keratocytes.rds')

#seurat_obj1 <- readRDS('testrun_allcells.rds')
#setwd('/home/bnvlab2/Documents/Kate/Cornea/Cells_Subset2/WGCNA')
