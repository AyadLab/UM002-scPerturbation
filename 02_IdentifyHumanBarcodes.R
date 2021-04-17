library(Seurat)

obj <- readRDS(file = "01_output/CellTagInt_Barnyard_UMAP.Rds")

# DimPlot(obj, group.by = "seurat_clusters", label = TRUE)
# DimPlot(obj, cells.highlight = WhichCells(obj, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "14", "15", "16", "18", "19", "20")))

# VlnPlot(obj, features = c("percent.hg19"))
human <- subset(obj, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "14", "15", "16", "18", "19", "20"))

barcodes <- data.frame(barcodes = colnames(human))
for (i in 1:length(barcodes$barcodes)){
  barcodes$arm[i] <- as.character(unlist(strsplit(as.character(barcodes$barcodes[i]), split = "_", fixed = TRUE))[1])
  barcodes$newbarcodes[i] <- as.character(unlist(strsplit(as.character(barcodes$barcodes[i]), split = "_", fixed = TRUE))[2])
  barcodes$transfer[i] <- as.character(unlist(strsplit(as.character(barcodes$barcodes[i]), split = "-", fixed = TRUE))[1])
}

DMSO_barcodes <- subset(barcodes, arm == "DMSO")
NM002_barcodes <- subset(barcodes, arm == "NM002")

DMSO_barcodes2 <- DMSO_barcodes[,c("newbarcodes")]
write.table(DMSO_barcodes2, file = "02_output/DMSO_hu_barcodes.tsv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
NM002_barcodes2 <- NM002_barcodes[,c("newbarcodes")]
write.table(NM002_barcodes2, file = "02_output/NM002_hu_barcodes.tsv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(DMSO_barcodes, file = "02_output/DMSO_hu_barcodes.csv", row.names = FALSE)
write.csv(NM002_barcodes, file = "02_output/NM002_hu_barcodes.csv", row.names = FALSE)

################################################################################
################################################################################
################################################################################

print("Pipeline successfully completed...")

################################################################################
################################################################################
################################################################################
