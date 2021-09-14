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

obj.features <- SelectIntegrationFeatures(
  object.list = CellTagGBM.list,
  nfeatures = 3000
)

CellTagGBM.list <- PrepSCTIntegration(
  object.list = CellTagGBM.list,
  anchor.features = obj.features,
  verbose = TRUE)

obj.anchors <- FindIntegrationAnchors(
  object.list = CellTagGBM.list,
  normalization.method = "SCT",
  anchor.features = obj.features,
  verbose = TRUE
  )
to_integrate <- Reduce(intersect, lapply(obj.anchors@object.list, rownames))
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

################################################################################
################################################################################
################################################################################

saveRDS(obj, file = "03_output/CellTagHuInt_UMAP.Rds")

# obj <- readRDS(file = "03_output/obj_UMAP.Rds")

print("Pipeline successfully completed...")

################################################################################
################################################################################
################################################################################
