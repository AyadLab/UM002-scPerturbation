# devtools::install_github(repo = "satijalab/seurat", ref = "develop")

library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggplotify)
options(future.globals.maxSize = 4000 * 1024^2)

################################################################################

# Load 10X data & Create Seurat Object...

DMSO_RdI.data <- Read10X(data.dir = "CellRangerOuts/RdI_DMSO_barnyard_outs/")
NM002_RdI.data <- Read10X(data.dir = "CellRangerOuts/RdI_NM002_barnyard_outs/")
DMSO_RdII.data <- Read10X(data.dir = "CellRangerOuts/RdII_DMSO_barnyard_outs/filtered_feature_bc_matrix/")
NM002_RdII.data <- Read10X(data.dir = "CellRangerOuts/RdII_UM002_barnyard_outs/filtered_feature_bc_matrix/")

################################################################################

# DMSO
######

RdI_DMSO <- CreateSeuratObject(
  counts = DMSO_RdI.data,
  min.cells = 3,
  min.features = 200,
  project = "DMSO_RdI"
)

RdII_DMSO <- CreateSeuratObject(
  counts = DMSO_RdII.data,
  min.cells = 3,
  min.features = 200,
  project = "DMSO_RdII"
)

RdI_DMSO$Arm <- "DMSO"
RdI_DMSO <- RenameCells(
  object = RdI_DMSO,
  new.names = paste0("DMSO_I_", colnames(x = RdI_DMSO))
)

RdII_DMSO$Arm <- "DMSO"
RdII_DMSO <- RenameCells(
  object = RdII_DMSO,
  new.names = paste0("DMSO_II_", colnames(x = RdII_DMSO))
)

# Cell Filtering QC
# RdI

RdI_DMSO$percent.mt <- PercentageFeatureSet(RdI_DMSO, pattern = "^MT-")
RdI_DMSO[["percent.hg19"]] <- PercentageFeatureSet(RdI_DMSO, pattern = "^hg19-")
RdI_DMSO[["percent.mm10"]] <- PercentageFeatureSet(RdI_DMSO, pattern = "^mm10-")

pdf(
  file = "01_output/RdI_DMSO_prefilterQC.pdf"
)
VlnPlot(RdI_DMSO, features = c("percent.hg19", "percent.mm10"))
FeatureScatter(object = RdI_DMSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdI_DMSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# RdII

RdII_DMSO$percent.mt <- PercentageFeatureSet(RdII_DMSO, pattern = "^MT-")
RdII_DMSO[["percent.hg19"]] <- PercentageFeatureSet(RdII_DMSO, pattern = "^hg19-")
RdII_DMSO[["percent.mm10"]] <- PercentageFeatureSet(RdII_DMSO, pattern = "^mm10-")

pdf(
  file = "01_output/RdII_DMSO_prefilterQC.pdf"
)
VlnPlot(RdII_DMSO, features = c("percent.hg19", "percent.mm10"))
FeatureScatter(object = RdII_DMSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdII_DMSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

################################################################################

# NM002
#######

RdI_NM002 <- CreateSeuratObject(
  counts = NM002_RdI.data,
  min.cells = 3,
  min.features = 200,
  project = "NM002_RdI"
)

RdII_NM002 <- CreateSeuratObject(
  counts = NM002_RdII.data,
  min.cells = 3,
  min.features = 200,
  project = "NM002_RdII"
)

RdI_NM002$Arm <- "NM002"
RdI_NM002 <- RenameCells(
  object = RdI_NM002,
  new.names = paste0("NM002_I_", colnames(x = RdI_NM002))
)

RdII_NM002$Arm <- "NM002"
RdII_NM002 <- RenameCells(
  object = RdII_NM002,
  new.names = paste0("NM002_II_", colnames(x = RdII_NM002))
)

# Cell Filtering QC
# RdI

RdI_NM002$percent.mt <- PercentageFeatureSet(RdI_NM002, pattern = "^MT-")
RdI_NM002[["percent.hg19"]] <- PercentageFeatureSet(RdI_NM002, pattern = "^hg19-")
RdI_NM002[["percent.mm10"]] <- PercentageFeatureSet(RdI_NM002, pattern = "^mm10-")

pdf(
  file = "01_output/RdI_NM002_prefilterQC.pdf"
)
VlnPlot(RdI_NM002, features = c("percent.hg19", "percent.mm10"))
FeatureScatter(object = RdI_NM002, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = RdI_NM002, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# RdII

RdII_NM002$percent.mt <- PercentageFeatureSet(RdII_NM002, pattern = "^MT-")
RdII_NM002[["percent.hg19"]] <- PercentageFeatureSet(RdII_NM002, pattern = "^hg19-")
RdII_NM002[["percent.mm10"]] <- PercentageFeatureSet(RdII_NM002, pattern = "^mm10-")

pdf(
  file = "01_output/RdII_NM002_prefilterQC.pdf"
)
VlnPlot(RdII_NM002, features = c("percent.hg19", "percent.mm10"))
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

# And uniformly filter all datasets...

for (i in 1:length(CellTagGBM.list)) {
  CellTagGBM.list[[i]] <- subset(x = CellTagGBM.list[[i]],
                                 subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
}


# Run scTransform on all objects...

for (i in 1:length(CellTagGBM.list)) {
  CellTagGBM.list[[i]] <- SCTransform(CellTagGBM.list[[i]], verbose = TRUE)
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

obj <- IntegrateData(
  anchorset = obj.anchors,
  #  features.to.integrate = to_integrate,
  normalization.method = "SCT",
  verbose = TRUE
)

rm(obj.anchors)
rm(obj.features)
rm(CellTagGBM.list)
rm(to_integrate)

saveRDS(obj, file = "01_output/CellTagInt.RDS")
obj1 <- readRDS(file = "01_output/CellTagInt.RDS")
dim(obj1)
table(obj1$Arm)
rm(obj1)
# Now proceed with downstream analysis (i.e. visualization, clustering) on the
# integrated dataset. Commands are identical to the standard workflow, but do
# not run the ScaleData function after integration. This should have helped to
# remove any variation caused by technical batch effect.

################################################################################
################################################################################
################################################################################

obj <- RunPCA(obj, verbose = FALSE)

pdf(
  file = "01_output/CellTagInt_ElbowPlot.pdf"
)
ElbowPlot(obj)
dev.off()

obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, dims = 1:15)

obj <- RunUMAP(
  obj,
  dims = 1:15,
  min.dist = 0.5,
  n_neighbors = 30,
  umap.method = "uwot",
  metric = "correlation",
  assay = "SCT"
)

pdf(
  file = "01_output/CellTagIntUMAP.pdf"
)
DimPlot(obj, group.by = c("Arm"), cols = c("azure4", "chartreuse3"))#, combine = FALSE)
DimPlot(obj, group.by = "ident", label = TRUE)#, combine = FALSE)
FeaturePlot(obj, features = "percent.hg19", cols = viridis(100, option = "D"))
FeaturePlot(obj, features = "percent.mm10", cols = viridis(100, option = "D"))

dev.off()

saveRDS(obj, file = "01_output/CellTagInt_Barnyard_UMAP.Rds")

################################################################################

# Generate supplemental figure showing seperation of pdx and mouse cells. 

# source umap

# sourceUMAP <- as.ggplot(DimPlot(obj, group.by = "orig.ident") + ggtitle(label = "Original Identity"))
# armUMAP <- as.ggplot(DimPlot(obj, group.by = c("Arm"), cols = c("azure4", "chartreuse3")) + ggtitle(label = "Treatment"))
saveRDS(obj, file = "01_output/barnyardObj.RDS")

snnUMAP <- as.ggplot(DimPlot(obj, group.by = "ident", label = TRUE, cols = DiscretePalette(20, palette = "alphabet2"), pt.size = 0.01))
hg19UMAP <- as.ggplot(FeaturePlot(obj, features = "percent.hg19", cols = viridis(100, option = "D"), pt.size = 0.01))# + ggtitle(label = "Percent alignment to hg19"))
mm10UMAP <- as.ggplot(FeaturePlot(obj, features = "percent.mm10", cols = viridis(100, option = "D"), pt.size = 0.01))# + ggtitle(label = "Percent alignment to mm10"))

armUMAP <- as.ggplot(DimPlot(obj, group.by = "Arm", label = F, cols = c("Azure4", "Chartreuse"), pt.size = 0.001, shuffle = TRUE))

sample10xUMAP <- as.ggplot(DimPlot(obj, group.by = "orig.ident", label = F, cols = c("yellow", "blue", "red", "violet"), pt.size = 0.01, shuffle = TRUE))
# pl <- as.ggplot(DimPlot(obj, group.by = "orig.ident", label = F, cols = DiscretePalette(2, "alphabet"), pt.size = 0.05))
sample10xUMAP

# Add PieChart of proportions of cells from each treatment arm. 

armdf <- as.data.frame(table(int$Arm))
armdf$perc <- armdf$Freq / sum(armdf$Freq) * 100
armtab <- armdf$perc
names(armtab) <- armdf$Var1
lbls <- lbls <- paste(names(armtab), "\n", armtab, sep = "")

pdf(file = "01_output/arm_proportionPie.pdf")
pie3D(armtab, explode = 0.1, labels = lbls, main="Treatment", start = 2, col = c("lightgrey", "chartreuse"))
dev.off()
# armUMAP
# DimPlo
# 
# DimPlot(obj, cols = DiscretePalette(20, palette = "alphabet2"))
dev.off()
pdf(file = "01_output/mm10_hg19_seperation.pdf", width = 14, height = 14)
plot_grid(plotlist = list(armUMAP,
                          sample10xUMAP,
                          hg19UMAP, 
                          mm10UMAP
                          ), ncol = 2, labels = "AUTO")
dev.off()


################################################################################
################################################################################
################################################################################

print("Pipeline successfully completed...")

################################################################################
################################################################################
################################################################################
