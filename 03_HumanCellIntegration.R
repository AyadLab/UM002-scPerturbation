# devtools::install_github(repo = "satijalab/seurat", ref = "develop")

# Need to account for duplicate cell names...?

library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)

################################################################################

# Load 10X data & Create Seurat Objects...

RdI_DMSO.data <- Read10X(
  data.dir = "CellRangerOuts/RdI_DMSO_GRCh38/filtered_feature_bc_matrix/")
RdII_DMSO.data <- Read10X(
  data.dir = "CellRangerOuts/RdII_DMSO_GRCh38/filtered_feature_bc_matrix/")

RdI_NM002.data <- Read10X(
  data.dir = "CellRangerOuts/RdI_NM002_GRCh38/filtered_feature_bc_matrix/")
RdII_NM002.data <- Read10X(
  data.dir = "CellRangerOuts/RdII_UM002_GRCh38/filtered_feature_bc_matrix/")

################################################################################

# DMSO
######

# Load in human cell barcode list from barnyard analysis
DMSO_human_barcodes <- as.character(read.csv(file = "02_output/DMSO_hu_barcodes.csv")$transfer)
head(DMSO_human_barcodes)

DMSO_human_barcodes <- paste0(DMSO_human_barcodes, "-1")

head(colnames(RdI_DMSO.data))

RdI_DMSO <- CreateSeuratObject(
  counts = RdI_DMSO.data,
  min.cells = 3,
  min.features = 200,
  project = "RdI_DMSO"
)
RdII_DMSO <- CreateSeuratObject(
  counts = RdII_DMSO.data,
  min.cells = 3,
  min.features = 200,
  project = "RdII_DMSO"
)

RdI_DMSO$Arm <- "DMSO"
RdII_DMSO$Arm <- "DMSO"

RdI_DMSO <- RenameCells(
  object = RdI_DMSO,
  new.names = paste0("DMSO_I_", colnames(x = RdI_DMSO))
)
RdII_DMSO <- RenameCells(
  object = RdII_DMSO,
  new.names = paste0("DMSO_II_", colnames(x = RdII_DMSO))
)

head(colnames(RdI_DMSO))
head(DMSO_human_barcodes)
RdI_DMSO <- subset(RdI_DMSO, cells = DMSO_human_barcodes)
RdII_DMSO <- subset(RdII_DMSO, cells = DMSO_human_barcodes)

# Check total # of cells matches number of barcodes...

length(DMSO_human_barcodes)
length(colnames(RdI_DMSO)) + length(colnames(RdII_DMSO))
# Cell Filtering QC

RdI_DMSO$percent.mt <- PercentageFeatureSet(RdI_DMSO, pattern = "^MT-")
RdII_DMSO$percent.mt <- PercentageFeatureSet(RdII_DMSO, pattern = "^MT-")

pdf(
  file = "03_output/RdI_DMSOhu_prefilterQC.pdf"
)
FeatureScatter(object = RdI_DMSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdI_DMSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf(
  file = "03_output/RdII_DMSOhu_prefilterQC.pdf"
)
FeatureScatter(object = RdII_DMSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdII_DMSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

################################################################################
################################################################################
################################################################################

# NM002
#######

# Load in human cell barcode list from barnyard analysis
NM002_human_barcodes <- as.character(read.csv(file = "02_output/NM002_hu_barcodes.csv")$transfer)
head(NM002_human_barcodes)

NM002_human_barcodes <- paste0(NM002_human_barcodes, "-1")

head(colnames(RdI_NM002.data))

RdI_NM002 <- CreateSeuratObject(
  counts = RdI_NM002.data,
  min.cells = 3,
  min.features = 200,
  project = "RdI_NM002"
)
RdII_NM002 <- CreateSeuratObject(
  counts = RdII_NM002.data,
  min.cells = 3,
  min.features = 200,
  project = "RdII_NM002"
)

RdI_NM002$Arm <- "NM002"
RdII_NM002$Arm <- "NM002"

RdI_NM002 <- RenameCells(
  object = RdI_NM002,
  new.names = paste0("NM002_I_", colnames(x = RdI_NM002))
)
RdII_NM002 <- RenameCells(
  object = RdII_NM002,
  new.names = paste0("NM002_II_", colnames(x = RdII_NM002))
)

head(colnames(RdI_NM002))
head(NM002_human_barcodes)
RdI_NM002 <- subset(RdI_NM002, cells = NM002_human_barcodes)
RdII_NM002 <- subset(RdII_NM002, cells = NM002_human_barcodes)

# Check total # of cells matches number of barcodes...

length(NM002_human_barcodes)
length(colnames(RdI_NM002)) + length(colnames(RdII_NM002))
# Cell Filtering QC

RdI_NM002$percent.mt <- PercentageFeatureSet(RdI_NM002, pattern = "^MT-")
RdII_NM002$percent.mt <- PercentageFeatureSet(RdII_NM002, pattern = "^MT-")

pdf(
  file = "03_output/RdI_NM002hu_prefilterQC.pdf"
)
FeatureScatter(object = RdI_NM002, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdI_NM002, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf(
  file = "03_output/RdII_NM002hu_prefilterQC.pdf"
)
FeatureScatter(object = RdII_NM002, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdII_NM002, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
################################################################################
################################################################################
################################################################################

##                         INTEGRATION

################################################################################
################################################################################
################################################################################

CellTagGBM.list <- c(RdI_DMSO, RdII_DMSO, RdI_NM002, RdII_NM002)
rm(RdI_DMSO, RdII_DMSO, RdI_NM002, RdII_NM002, RdI_DMSO.data, RdII_DMSO.data, RdI_NM002.data, RdII_NM002.data)
# And uniformly filter all datasets...

for (i in 1:length(CellTagGBM.list)) {
  CellTagGBM.list[[i]] <- subset(x = CellTagGBM.list[[i]],
    subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
}


# Run scTransform on all objects...

for (i in 1:length(CellTagGBM.list)) {
  CellTagGBM.list[[i]] <- SCTransform(CellTagGBM.list[[i]], verbose = TRUE, return.only.var.genes	= FALSE)
}



# Next, select features for downstream integration, and run PrepSCTIntegration,
# which ensures that all necessary Pearson residuals have been calculated.

obj.features <- SelectIntegrationFeatures(
  object.list = CellTagGBM.list,
  nfeatures = 3000
)

CellTagGBM.list <- PrepSCTIntegration(
  object.list = CellTagGBM.list,
  anchor.features = obj.features,
  verbose = TRUE)

# Next, identify anchors and integrate the datasets. Commands are identical to
# the standard workflow, but make sure to set normalization.method = "SCT":

obj.anchors <- FindIntegrationAnchors(
  object.list = CellTagGBM.list,
  normalization.method = "SCT",
  anchor.features = obj.features,
  verbose = TRUE
  )

to_integrate <- Reduce(intersect, lapply(obj.anchors@object.list, rownames))

print("Preview of genes to integrate... ")
head(to_integrate)

rm(CellTagGBM.list)

obj <- IntegrateData(
  anchorset = obj.anchors,
  features.to.integrate = to_integrate,
  normalization.method = "SCT",
  verbose = TRUE
)

rm(obj.anchors)
rm(obj.features)
rm(to_integrate)

saveRDS(obj, file = "03_output/CellTagHuInt.RDS")

# Now proceed with downstream analysis (i.e. visualization, clustering) on the
# integrated dataset. Commands are identical to the standard workflow, but do
# not run the ScaleData function after integration. This should have helped to
# remove any variation caused by technical batch effect.

################################################################################
################################################################################
################################################################################

obj <- RunPCA(obj, verbose = FALSE)

pdf(
    file = "03_output/CellTagHuInt_ElbowPlot.pdf"
)
ElbowPlot(obj)
dev.off()

obj <- FindNeighbors(obj, dims = 1:15)
# obj <- FindClusters(obj, dims = 1:15, )
obj <- FindClusters(obj, dims = 1:15, resolution = seq(from = 0.9, to = 0.05, by = -0.05), algorithm = 2)

obj <- RunUMAP(
  obj,
  dims = 1:15,
  min.dist = 0.5,
  n.neighbors = 30,
  n.epochs = 1000,
  # umap.method = "uwot",
  metric = "correlation",
  assay = "integrated"
)

pdf(
  file = "03_output/CellTagHuIntUMAP.pdf"
)
DimPlot(obj, group.by = c("Arm"))#, combine = FALSE)
DimPlot(obj, group.by = "ident", label = TRUE)#, combine = FALSE)
dev.off()

# # Identify celltypes - due to homology, some mouse cells are going to align and
# # mess with our data, let's remove them. Our tumor cells should cluster seperately
# # from any mouse cells, and this can be biologically checked.
#
# DefaultAssay(obj) <- "RNA"
#
# pdf(
#   file = "03_output/CellTypeMarkers.pdf"
# )
# FeaturePlot(obj, features = c("EGFR"))
# FeaturePlot(obj, features = c("CD74"))
# FeaturePlot(obj, features = c("CLDN5"))
# FeaturePlot(obj, features = c("PTPRC"))
# FeaturePlot(obj, features = c("PDGFRB"))
# FeaturePlot(obj, features = c("AQP4"))
# FeaturePlot(obj, features = c("MBP"))
# FeaturePlot(obj, features = c("HBB"))
# FeaturePlot(obj, features = c("PDGFRA"))
# FeaturePlot(obj, features = c("CDK4"))
# dev.off()
#
# ################################################################################
# ################################################################################
# ################################################################################
#
# # Split into wha we believe are neoplastic and contaminating mouse cells..
# # Then run CytoTRACE and see what it loosk like...
#
# Idents(obj) <- obj$seurat_clusters
# obj <- RenameIdents(obj,
#    "0" = "Mouse",
#    "1" = "Neoplastic",
#    "2" = "Neoplastic",
#    "3" = "Neoplastic",
#    "4" = "Neoplastic",
#    "5" = "Neoplastic",
#    "6" = "Neoplastic",
#    "7" = "Neoplastic",
#    "8" = "Neoplastic",
#    "9" = "Mouse",
#    "10" = "Mouse",
#    "11" = "Neoplastic",
#    "12" = "Mouse",
#    "13" = "Mouse",
#    "14" = "Mouse",
#    "15" = "Neoplastic"
#  )
# obj$CellType <- Idents(obj)
#
# pdf(file = "03_output/CellTagSpeciesSplitUMAP.pdf")
# DimPlot(obj, group.by = "CellType")
# dev.off()
#
# ## Run CytoTRACE
#
# DefaultAssay(obj) <- "RNA"
# counts_matrix <- as.data.frame(GetAssayData(obj, slot = "counts"))
#
# # Run CytoTRACE
#
# library(CytoTRACE)
# results <- CytoTRACE(counts_matrix, ncores = 24)
# saveRDS(results, file = "03_output/CellTag_CytoTRACE_Results.RDS")
#
# z <- as.character(obj$CellType)
# names(z) <- names(obj$CellType)
#
# # Plot will be saved below
# plotCytoTRACE(results, phenotype = z)
#
#
#
# ################################################################################
#
# # Now subset the neoplastic human cells - and run differential expression
# # between each arm as a whole?
#
# Idents(obj) <-obj$CellType
# obj <- subset(obj, idents = "Neoplastic")
#
# ################################################################################
#
# # CellCycleScoring, Rescale & Regression, re-UMAP, etc.
#
# ################################################################################
#
# # Cell Cycle Scoring
#
# DefaultAssay(obj) <- "RNA"
# obj <- CellCycleScoring(obj,
#   s.features = cc.genes$s.genes,
#   g2m.features = cc.genes$g2m.genes,
#   set.ident = TRUE)
#
# obj$CellCycle <- Idents(obj)
#
# # Rerun dimensionality reduction on human cells
#
# DefaultAssay(obj) <- "SCT"
# obj <- FindVariableFeatures(obj)
# obj <- RunPCA(obj, verbose = FALSE)
#
# pdf(file = "03_output/CellTagTumorElbowPlot_noreg.pdf")
# ElbowPlot(obj)
# dev.off()
#
# obj <- FindNeighbors(obj, dims = 1:15)
# obj <- FindClusters(obj, dims = 1:15)
# obj <- RunUMAP(
#   obj,
#   dims = 1:15,
#   min.dist = 0.5,
#   n_neighbors = 30,
#   umap.method = "uwot",
#   metric = "correlation",
#   assay = "SCT"
# )
# pdf(file = "03_output/CellTagTumorUMAP_NoReg.pdf")
# DimPlot(obj, group.by = "CellCycle")
# dev.off()
#
# # Regress out cell cycle...
#
# DefaultAssay(obj) <- "integrated"
# obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), assay = "integrated")
#
# # Rerun dimensionality reduction...
#
# obj <- FindVariableFeatures(obj)
# obj <- RunPCA(obj, verbose = FALSE, assay = "integrated")
#
# pdf(file = "03_output/CellTagTumorElbowPlot_reg.pdf")
# ElbowPlot(obj)
# dev.off()
#
# obj <- FindNeighbors(obj, dims = 1:15, assay = "SCT")
# obj <- FindClusters(obj, dims = 1:15)
# obj <- RunUMAP(
#   obj,
#   dims = 1:15,
#   min.dist = 0.5,
#   n_neighbors = 30,
#   umap.method = "uwot",
#   metric = "correlation",
#   assay = "integrated"
# )
#
# pdf(file = "03_output/CellTagTumorUMAP_Reg.pdf")
# DimPlot(obj, group.by = "CellCycle")
# dev.off()
#
# ################################################################################
# ################################################################################
# ################################################################################
#
#
# pdf(file = "03_output/subsettedUMAP.pdf")
# DimPlot(obj, group.by = "CellCycle")
# DimPlot(obj, group.by = "Arm")
# dev.off()
# # library(MAST)
#
# Idents(obj) <- obj$Arm
#
# DefaultAssay(obj) <- "RNA"
#
# ArmMarkers <- FindAllMarkers(
#   object = obj,
#   test = "wilcox",
#   assay = "RNA",
#   features = NULL,
#   logfc.threshold = 0.25,
#   slot = "data")
#
# write.csv(
#   ArmMarkers,
#   file="ArmMarkers.csv",
#   row.names=TRUE,
#   col.names=TRUE
# )
#
# obj <- ScaleData(obj, assay = "RNA", vars.to.regress = c("S.Score", "G2M.Score"))
#
# library(pheatmap)
#
# mat <- as.data.frame(GetAssayData(obj, assay = "RNA", slot = "scale.data"))
# head(rownames(mat))
# mat <- subset(mat, rownames(mat) %in% ArmMarkers$gene)
#
# saveRDS(mat, file = "03_output/mat_for_heatmap.RDS")
#
# pdf(file = "03_output/ArmHeatmapDEx.pdf")
# pheatmap(mat, annotation_col = as.data.frame(obj$Arm), scale = "none", show_colnames = FALSE, cluster_cols = TRUE)
# dev.off()
#
# ################################################################################
# ################################################################################
# ################################################################################
#
# # Plot average differential expression between both arms...
#
# library(dplyr)
# library(ggplot2)
# library(cowplot)
# Idents(obj) <- obj$Arm
# head(Idents(obj))
# theme_set(theme_cowplot())
# # DMSO <- subset(obj, idents = "DMSO")
# # Idents(t.cells) <- "stim"
# avg.obj <- log1p(AverageExpression(obj, verbose = FALSE)$RNA)
# avg.obj$gene <- rownames(avg.obj)
#
# # cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
# # Idents(cd14.mono) <- "stim"
# # avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
# # avg.cd14.mono$gene <- rownames(avg.cd14.mono)
#
# ArmMarkers <- read.csv(file = "03_output/ArmMarkers.csv")
# head(ArmMarkers)
# top <- ArmMarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#
# genes.to.label = top$gene
# # genes.to.label = c("AURKA")
# p1 <- ggplot(avg.obj, aes(NM002, DMSO)) + geom_point() + ggtitle("Total sample average")
# p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
# # p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
# # p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
#
# pdf(file = "03_output/DMSO_vs_NM002_avg_expression_plot.pdf")
# plot_grid(p1)
# dev.off()
#
# ################################################################################
# ################################################################################
# ################################################################################
#
# # Are these genes cell type specific? I.e NPC, OPC, MES, AClike... etc.
#
################################################################################

# Split into Suva Subtypes for differential expression analysis.

################################################################################

################################################################################

# Suva Subtypes
print("Assign Suva Subtypes")
################################################################################

Suva.Sigs <- read.csv(
  file = "03_output/SuvaMetaSignatures.csv",
  skip = 4)

################################################################################
################################################################################
################################################################################

DefaultAssay(obj) <- "SCT"

################################################################################
################################################################################
################################################################################

# MES-like2

Suva.MES2.sig <- list(Suva.Sigs$MES2)

obj <- AddModuleScore(
  object = obj,
  features = Suva.MES2.sig,
  assay = "SCT",
  name = "Suva.MES2.Signature")

################################################################################
################################################################################
################################################################################

# MES-like1

Suva.MES1.sig <- list(Suva.Sigs$MES1)

obj <- AddModuleScore(
  object = obj,
  features = Suva.MES1.sig,
  assay = "SCT",
  name = "Suva.MES1.Signature")

################################################################################
################################################################################
################################################################################

# AC-like

Suva.AC.sig <- list(Suva.Sigs$AC)

obj <- AddModuleScore(
  object = obj,
  features = Suva.AC.sig,
  assay = "SCT",
  name = "Suva.AC.Signature")

################################################################################
################################################################################
################################################################################

# OPC-like

Suva.OPC.sig <- list(Suva.Sigs$OPC)

obj <- AddModuleScore(
  object = obj,
  features = Suva.OPC.sig,
  assay = "SCT",
  name = "Suva.OPC.Signature")

################################################################################
################################################################################
################################################################################

# NPC1-like

Suva.NPC1.sig <- list(Suva.Sigs$NPC1)

obj <- AddModuleScore(
  object = obj,
  features = Suva.NPC1.sig,
  assay = "SCT",
  name = "Suva.NPC1.Signature")

################################################################################
################################################################################
################################################################################

# NPC2-like

Suva.NPC2.sig <- list(Suva.Sigs$NPC2)

obj <- AddModuleScore(
  object = obj,
  features = Suva.NPC2.sig,
  assay = "SCT",
  name = "Suva.NPC2.Signature")

################################################################################
################################################################################
################################################################################

# G1/S

Suva.G1_S.sig <- list(Suva.Sigs$G1_S)

obj <- AddModuleScore(
  object = obj,
  features = Suva.G1_S.sig,
  assay = "SCT",
  name = "Suva.G1_S.Signature")

################################################################################
################################################################################
################################################################################

# G2/M

Suva.G2_M.sig <- list(Suva.Sigs$G2_M)

obj <- AddModuleScore(
  object = obj,
  features = Suva.G2_M.sig,
  assay = "SCT",
  name = "Suva.G2_M.Signature")

################################################################################
################################################################################
################################################################################

################################################################################

# Clustree Figure

library(clustree)

pdf(file = "03_output/CellTagHu_Clustree_Figures.pdf")
DimPlot(obj, group.by = "Arm")
clustree(obj, prefix = "integrated_snn_res.")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.OPC.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.NPC1.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.NPC2.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.MES1.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.MES2.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.AC.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.G2_M.Signature1", node_colour_aggr = "mean")
clustree(obj, prefix = "integrated_snn_res.", node_colour = "Suva.G1_S.Signature1", node_colour_aggr = "mean")
dev.off()

###############################################################################
###############################################################################
###############################################################################

# How do we run differential expression between suva Subtypes

Idents(obj) <- "AMBIGUOUS"

# 3) There will be cells expressing multiple signatures, how do we address this?

################################################################################
################################################################################
################################################################################

# MES1

MES1SigScores <- FetchData(
  obj,
  vars = "Suva.MES1.Signature1"
)

MES1SigQuantiles <- quantile(
  x = MES1SigScores$Suva.MES1.Signature1
)

# head(MES1SigQuantiles)

SuvaMES1Cells <- WhichCells(
  obj,
  expression = Suva.MES1.Signature1 > MES1SigQuantiles[4])

Idents(obj, cells = SuvaMES1Cells) <- "MES1"

################################################################################
################################################################################
################################################################################

# MES2

MES2SigScores <- FetchData(
  obj,
  vars = "Suva.MES2.Signature1"
)

MES2SigQuantiles <- quantile(
  x = MES2SigScores$Suva.MES2.Signature1
)

SuvaMES2Cells <- WhichCells(
  obj,
  expression = Suva.MES2.Signature1 > MES2SigQuantiles[4]
)

Idents(obj, cells = SuvaMES2Cells) <- "MES2"

################################################################################
################################################################################
################################################################################

# AC

ACSigScores <- FetchData(
  obj,
  vars = "Suva.AC.Signature1"
)

ACSigQuantiles <- quantile(
  x = ACSigScores$Suva.AC.Signature1
)

SuvaACCells <- WhichCells(
  obj,
  expression = Suva.AC.Signature1 > ACSigQuantiles[4]
)

Idents(obj, cells = SuvaACCells) <- "AC"

################################################################################
################################################################################
################################################################################

# OPC

OPCSigScores <- FetchData(
  obj,
  vars = "Suva.OPC.Signature1"
)

OPCSigQuantiles <- quantile(
  x = OPCSigScores$Suva.OPC.Signature1
)

SuvaOPCCells <- WhichCells(
  obj,
  expression = Suva.OPC.Signature1 > OPCSigQuantiles[4]
)

Idents(obj, cells = SuvaOPCCells) <- "OPC"

################################################################################
################################################################################
################################################################################

# NPC1

NPC1SigScores <- FetchData(
  obj,
  vars = "Suva.NPC1.Signature1"
)

NPC1SigQuantiles <- quantile(
  x = NPC1SigScores$Suva.NPC1.Signature1
)

SuvaNPC1Cells <- WhichCells(
  obj,
  expression = Suva.NPC1.Signature1 > NPC1SigQuantiles[4]
)

Idents(obj, cells = SuvaNPC1Cells) <- "NPC1"

################################################################################
################################################################################
################################################################################

# NPC2

NPC2SigScores <- FetchData(
  obj,
  vars = "Suva.NPC2.Signature1"
)

NPC2SigQuantiles <- quantile(
  x = NPC2SigScores$Suva.NPC2.Signature1
)

SuvaNPC2Cells <- WhichCells(
  obj,
  expression = Suva.NPC2.Signature1 > NPC2SigQuantiles[4]
)

Idents(obj, cells = SuvaNPC2Cells) <- "NPC2"

################################################################################
################################################################################
################################################################################

# G1_S

G1_SSigScores <- FetchData(
  obj,
  vars = "Suva.G1_S.Signature1"
)

G1_SSigQuantiles <- quantile(
  x = G1_SSigScores$Suva.G1_S.Signature1
)

SuvaG1_SCells <- WhichCells(
  obj,
  expression = Suva.G1_S.Signature1 > G1_SSigQuantiles[4]
)

# Idents(obj, cells = SuvaG1_SCells) <- "G1_S"

################################################################################
################################################################################
################################################################################

# G2_M

G2_MSigScores <- FetchData(
  obj,
  vars = "Suva.G2_M.Signature1"
)

G2_MSigQuantiles <- quantile(
  x = G2_MSigScores$Suva.G2_M.Signature1
)

SuvaG2_MCells <- WhichCells(
  obj,
  expression = Suva.G2_M.Signature1 > G2_MSigQuantiles[4]
)

# Idents(obj, cells = SuvaG2_MCells) <- "G2_M"

################################################################################
################################################################################
################################################################################

# And how does this look?
pdf(file = "03_output/IntSubtypeUMAP.pdf")
DimPlot(
  object = obj,
  reduction = "umap",
  pt.size = 1.0
)
dev.off()

################################################################################
################################################################################
################################################################################

obj$Subtype <- Idents(obj)

SuvaModules.list <- c(
  "Suva.MES1.Signature1",
  "Suva.MES2.Signature1",
  "Suva.AC.Signature1",
  "Suva.OPC.Signature1",
  "Suva.NPC1.Signature1",
  "Suva.NPC2.Signature1",
  "Suva.G1_S.Signature1",
  "Suva.G2_M.Signature1"
  )

SuvaModules.list.2 <- c(
  "Suva.MES1.Signature1",
  "Suva.MES2.Signature1",
  "Suva.AC.Signature1",
  "Suva.OPC.Signature1",
  "Suva.NPC1.Signature1",
  "Suva.NPC2.Signature1",
  "Suva.G1_S.Signature1",
  "Suva.G2_M.Signature1"
  )

pdf(file = "03_output/CellTagSubtypeDotPlots.pdf")
DotPlot(
  obj,
  assay = NULL,
  features = SuvaModules.list,
  cols = c("cyan", "maroon"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 9,
  # group.by = "Idents",
  # split.by = NULL,
  scale.by = "radius"
  # scale.min = NA,
  # scale.max = NA
)
#
# obj$ArmSubtype <- paste0(obj$Subtype, "_", obj$Arm)
#
# DotPlot(
#   obj,
#   assay = NULL,
#   features = SuvaModules.list.2,
#   cols = c("cyan", "maroon"),
#   col.min = -2.5,
#   col.max = 2.5,
#   dot.min = 0,
#   dot.scale = 9,
#   group.by = "ArmSubtype",
#   # split.by = NULL,
#   scale.by = "radius"
#   # scale.min = NA,
#   # scale.max = NA
# )
# dev.off()
#
# ################################################################################
#
# #     Investigate subtype proportion shift from DMSO to NM002
#
# ################################################################################
#
# # Get cell counts...
#
# Idents(obj) <- obj$Subtype
# df <- table(Idents(obj), obj$Arm)
# df <- as.data.frame(table(Idents(obj), obj$Arm))
# df <- data.frame(df)
# library(reshape2)
# df <- recast(df, formula = Var1 ~ Var2, mean)
# rownames(df) <- df$Var1
# df$Var1 <- NULL
# DMSOcount <- sum(df$DMSO)
# NM002count <- sum(df$NM002)
#
# # Plottable dataframe... set percentage
#
# df.Arm <- data.frame(table(Idents(obj), obj$Arm))
# df.Arm
# df.Arm$Perc <- "?"
# length(df.Arm$Freq)
# for (i in 1:length(df.Arm$Freq)){
#                                 if (df.Arm$Var2[[i]] == "DMSO") {
#                                                             df.Arm$Perc[[i]] <- round(df.Arm$Freq[[i]] / DMSOcount * 100, digits = 2)
#                                 }
#                                 if (df.Arm$Var2[[i]] == "NM002") {
#                                                           df.Arm$Perc[[i]] <- round(df.Arm$Freq[[i]] / NM002count * 100, digits = 2)
#                                 }
#
# }
# df.Arm
# #
# library(ggplot2)
#
# # Barplot
# bp<- ggplot(df.Arm, aes(x=Var2, y=Perc, fill=Var1))+
# geom_bar(width = 1, stat = "identity")
#
# pdf(file = "03_output/CellTagSubtypeShiftBarplot.pdf")
# bp
# dev.off()
#
# saveRDS(obj, file = "03_output/CellTagInt_Subtype.RDS")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
################################################################################
################################################################################
################################################################################

saveRDS(obj, file = "03_output/CellTagHuInt_UMAP.Rds")

# obj <- readRDS(file = "03_output/obj_UMAP.Rds")

print("Pipeline successfully completed...")

################################################################################
################################################################################
################################################################################
