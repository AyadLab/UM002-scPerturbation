
library(Seurat)
library(dplyr)
library(CytoTRACE)
library(viridis)

# Load Seurat object and extract count matrix
obj <- readRDS(file = "03_output/CellTagHuInt_UMAP.Rds")
DefaultAssay(obj) <- "RNA"
counts_matrix <- as.data.frame(GetAssayData(obj, slot = "counts"))

# Save count matrix as mtx and rds files
# saveRDS(counts_matrix, "obj_counts_matrix.RDS")
# write.table(round(counts_matrix, digits = 3), file = "obj_counts.matrix", quote = F, sep = "\t")

# Extract cell-type annotations from Seurat object
obj_cell_type_anno <- as.data.frame(obj$integrated_snn_res.0.35)

# Run CytoTRACE

results <- CytoTRACE(counts_matrix, ncores = 24)
saveRDS(results, file = "04_output/04_CellTagHu_CytoTRACE_Results.RDS")
pheno <- as.character(obj$integrated_snn_res.0.35)
names(pheno) <- names(obj$integrated_snn_res.0.35)
# Get umap embeddings from Seurat...
umapEmb <- as.data.frame(Embeddings(obj, reduction = "umap"))
plotCytoTRACE(results, phenotype = pheno, gene = "TIMP1", emb = umapEmb)

# Integrate CytoTRACE results back into Seurat object.
obj <- AddMetaData(obj, metadata = results$CytoTRACE, col.name = "CytoTRACE")
obj <- AddMetaData(obj, metadata = results$CytoTRACErank, col.name = "CytoTRACErank")

pdf(
  file = "04_CellTagHu_CytoTRACE.pdf"
)
FeaturePlot(obj, features = c("CytoTRACE"), cols = viridis(256, option = "D"), reduction = "umap")
FeaturePlot(obj, features = c("CytoTRACE"), cols = viridis(256, option = "D"), reduction = "umap") + NoLegend()
VlnPlot(obj, features = c("CytoTRACE"), pt.size = 0.001, group.by = "Arm")
dev.off()

saveRDS(obj, file = "04_output/04_CellTagHu_CytoTRACE.RDS")
