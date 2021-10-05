
library(Seurat)
library(GSEABase)
library(cowplot)
library(singscore)
library(SummarizedExperiment)
library(dplyr)
library(ggpubr)
library(Rmisc)
library(ggplot2)
library(plotrix)
library(EnhancedVolcano)
library(ggsci)
library(reshape2)
library(clustree)

###########################################################################################

#' Plots a series of barplots and connects them
#' Modified from https://stackoverflow.com/questions/22560850/barplot-with-connected-series

connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {
  b <- barplot(dat, col=color, space = space, ...)

  for (i in seq_len(ncol(dat) - 1)) {
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line

    for (j in seq_len(nrow(dat))) {
      if (j == 1) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(0, dat[j,i], dat[j,i+1], 0),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
      if (j == 2) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
      if (j > 2) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1],
                  colSums(dat[1:(j-1),])[i+1]),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
    }
  }
}

# Add this line

###########################################################################################

int <- readRDS(file = "04_output/04_CellTagHu_CytoTRACE.RDS")

# scale RNA assay

int <- NormalizeData(int, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
DefaultAssay(int) <- "RNA"
int <- CellCycleScoring(int, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
int$CellCycleIdents <- Idents(int)
int <- ScaleData(int, assay = "RNA", features = rownames(int), vars.to.regress = c("G2M.Score", "S.Score", "percent.mt"))
# saveRDS(int, file = "05_output/int.scaled.g2m.s.percent.mt.RDS")
# int <- readRDS(file = "05_output/int.scaled.g2m.s.percent.mt.RDS")
Idents(int) <- int$Arm
NM <- subset(int, idents = "NM002")
DMSO <- subset(int, idents = "DMSO")

#######################################################################################
#######################################################################################
#######################################################################################

# Part I - Seurat module scoring

DefaultAssay(int) <- "RNA"
int <- CellCycleScoring(int, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
int$CellCycleIdents2 <- Idents(int)
CellCycleCounts <- table(int$CellCycleIdents, int$Arm)
DMSOsum <- sum(CellCycleCounts[,"DMSO"])
NMsum <- sum(CellCycleCounts[,"NM002"])
CellCycleCounts2 <- as.data.frame(CellCycleCounts)
colnames(CellCycleCounts2) <- c("Phase", "Treatment", "Freq")

for (i in 1:length(CellCycleCounts2$Freq)){
  if (CellCycleCounts2$Treatment[i] == "DMSO"){
    CellCycleCounts2$Percent[i] <- CellCycleCounts2$Freq[i]/DMSOsum*100
  }
  if (CellCycleCounts2$Treatment[i] == "NM002"){
    CellCycleCounts2$Percent[i] <- CellCycleCounts2$Freq[i]/NMsum*100
  }
}

CellCycleCounts3 <- recast(CellCycleCounts2, formula = Phase ~ Treatment, measure.var = "Percent")
rownames(CellCycleCounts3) <- CellCycleCounts3$Phase
CellCycleCounts3$Phase <- NULL
CellCycleCounts4 <- as.matrix(CellCycleCounts3)

pdf(file = "05_output/UM002_DMSO_CellCycleShift_connectedBarplot.pdf")
connectedBarplot(dat = CellCycleCounts4)
dev.off()

# Make a plot showing the changes in % of assigned cell cycle identities.
CellCycleCounts5 <- as.data.frame(CellCycleCounts4)

CellCycleCounts5$delta <- CellCycleCounts5$NM002 - CellCycleCounts5$DMSO
CellCycleCounts5$Phase <- rownames(CellCycleCounts5)

PhasePercShiftPlot <- ggplot(CellCycleCounts5, aes(x=Phase, y=delta, fill=Phase)) +
  geom_bar(aes(fill = Phase), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  scale_fill_jama(alpha = 0.7) +
  xlab("Assigned Cell Cycle Phase Identity") +
  ylab("Delta % Phase Proportion") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
        axis.title.x = element_text(size = rel(1.6), angle = 0))

pdf(file = "05_output/phasePercentShiftPlot.pdf")
PhasePercShiftPlot
dev.off()
####################################################################################################
pdf(file = "05_output/CellTag_Hu_CellCycleIdents.pdf")
DimPlot(int, group.by = "CellCycleIdents")
dev.off()

# Now, is this significant? growth of G1, loss of G2M and S?

CellCycleScoreDF <- data.frame(G2M.Score = int$G2M.Score, S.Score = int$S.Score, Treatment = int$Arm)
mlt <- melt(CellCycleScoreDF)
colnames(mlt) <- c("Treatment", "Phase", "Score")

CellCyclesum <- summarySE(mlt, measurevar = "Score", groupvars = c("Treatment", "Phase"))
CellCyclesum

pdf(file = "05_output/CellTTagHu_UM002_CellCyclePhaseShiftPairedBarplot.pdf")
# 95% CI instead of SE...
ggplot(CellCyclesum, aes(x=Phase, y=Score, fill=Treatment)) +
  geom_bar(aes(fill = Treatment), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  geom_errorbar(aes(ymin=Score-ci, ymax = Score+ci),
                width = 0.2,
                position = position_dodge(.9)) +
  scale_fill_manual(values = c("gray80", "chartreuse"),
                    name = "Treatment",
                    breaks = c("DMSO", "NM002"),
                    labels = c("DMSO", "UM002")) +
  xlab("Cell Cycle Phase Signature") +
  ylab("Module Score") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
        axis.title.x = element_text(size = rel(1.6), angle = 0))
dev.off()

g2mmlt <- subset(mlt, mlt$Phase == "G2M.Score")
g2mmlt

smlt <- subset(mlt, mlt$Phase == "S.Score")
smlt

g2mres <- t.test(Score ~ Treatment, data = g2mmlt, paired = FALSE)
g2mres

sres <- t.test(Score ~ Treatment, data = smlt, paired = FALSE)
sres

#########################################################################
#########################################################################
#########################################################################


##########################################################################################################



Idents(int) <- int$Arm
bulkArmMarkers <- FindAllMarkers(int, test.use = "MAST", only.pos = TRUE)

##########################################################################################################
##########################################################################################################
##########################################################################################################

# SNN Clustering Analysis

##########################################################################################################
##########################################################################################################
##########################################################################################################

int2 <- int
int2$integrated_snn_res.0.3 <- NULL
int2$integrated_snn_res.0.35 <- NULL
int2$integrated_snn_res.0.4 <- NULL
int2$integrated_snn_res.0.45 <- NULL
int2$integrated_snn_res.0.5 <- NULL
int2$integrated_snn_res.0.55 <- NULL
int2$integrated_snn_res.0.6 <- NULL
int2$integrated_snn_res.0.65 <- NULL
int2$integrated_snn_res.0.7 <- NULL
int2$integrated_snn_res.0.75 <- NULL
int2$integrated_snn_res.0.8 <- NULL
int2$integrated_snn_res.0.85 <- NULL
int2$integrated_snn_res.0.9 <- NULL

snnclustree <- clustree(int2, prefix = "integrated_snn_res.", layout = "tree", node_alpha = 1, node_text_colour = "white", node_text_size = 5, node_label_nudge = -0.7)
clustreeG2M <- clustree(int2, prefix = "integrated_snn_res.", node_colour = "G2M.Score", node_colour_aggr = "mean", layout = "tree", node_alpha = 1, node_text_size = 5, node_text_colour = "white")
clustreeS <- clustree(int2, prefix = "integrated_snn_res.", node_colour = "S.Score", node_colour_aggr = "mean", node_alpha = 1, node_text_size = 5, node_text_colour = "white")

pdf(file = "05_output/clustree.pdf", width = 21)
plot_grid(plotlist = list(snnclustree, clustreeG2M, clustreeS), ncol = 3, labels = "AUTO")
dev.off()

rm(int2)

Idents(int) <- int$integrated_snn_res.0.25
snn.0.25.markers <- FindAllMarkers(int, test.use = "MAST", only.pos = T)
write.csv(snn.0.25.markers, file = "05_output/snn.0.25.markers.csv")
top_snn_markers <- snn.0.25.markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
write.csv(top_snn_markers, file = "05_output/top_snn_markers.csv")

snnMarkersHeatmap <- DoHeatmap(int, features = top_snn_markers$gene, group.by = "integrated_snn_res.0.25") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")

pdf(file = "05_output/SNN.0.25.markers.heatmap.pdf", height = 14, width = 10)
snnMarkersHeatmap
dev.off()

##########################################################################################################

# How do these SNN clusters shift following treatment?

SNNcounts <- table(int$integrated_snn_res.0.25, int$Arm)

DMSOsum <- sum(SNNcounts[,"DMSO"])
NMsum <- sum(SNNcounts[,"NM002"])
SNNcounts <- as.data.frame(SNNcounts)
# SNNcounts
colnames(SNNcounts) <- c("SNN.Cluster", "Treatment", "Freq")

for (i in 1:length(SNNcounts$Freq)){
  if (SNNcounts$Treatment[i] == "DMSO"){
    SNNcounts$Percent[i] <- SNNcounts$Freq[i]/DMSOsum*100
  }
  if (SNNcounts$Treatment[i] == "NM002"){
    SNNcounts$Percent[i] <- SNNcounts$Freq[i]/NMsum*100
  }
}

# SNNcounts

library(reshape2)
SNNcounts <- recast(SNNcounts, formula = SNN.Cluster ~ Treatment, measure.var = "Percent")
rownames(SNNcounts) <- SNNcounts$SNN.Cluster
SNNcounts$SNN.Cluster <- NULL
SNNcounts <- as.matrix(SNNcounts)
# SNNcounts

pdf(file = "05_output/UM002_DMSO_SNNShift_connectedBarplot.pdf")
connectedBarplot(dat = SNNcounts)
dev.off()

pdf(file = "05_output/SNN.0.25_umap.pdf")
DimPlot(int, group.by = "integrated_snn_res.0.25", pt.size = 0.01, shuffle = T, cols = rainbow(8))
dev.off()

####################################################################################################

# Plot the change in percentage of each cluster...

SNNcounts2 <- as.data.frame(SNNcounts)

SNNcounts2$delta <- SNNcounts2$NM002 - SNNcounts2$DMSO
SNNcounts2$Sig <- rownames(SNNcounts2)

statePercShiftPlot <- ggplot(SNNcounts2, aes(x=Sig, y=delta, fill=Sig)) +
  geom_bar(aes(fill = Sig), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  scale_color_jama(alpha = 0.7) +
  xlab("Transcriptional State Signature") +
  ylab("Delta % State Proportion") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
        axis.title.x = element_text(size = rel(1.6), angle = 0))

pdf(file = "05_output/SNN_statePercentShiftPlot.pdf")
statePercShiftPlot
dev.off()

# Now export data to py, visualize effects on cell cycle by SNN cluster (also by transcriptional state, one at a time)

SNN_CellCycleDat <- FetchData(int, vars = c("Arm", "integrated_snn_res.0.25", "G2M.Score", "S.Score"))
for(i in 1:length(rownames(SNN_CellCycleDat))){
  SNN_CellCycleDat$armSNN[i] <- paste0(SNN_CellCycleDat$Arm[i],"_",SNN_CellCycleDat$integrated_snn_res.0.25[i])
  SNN_CellCycleDat$armSNN2[i] <- paste0(SNN_CellCycleDat$Arm[i],"_",SNN_CellCycleDat$integrated_snn_res.0.25[i])
  if (SNN_CellCycleDat$Arm[i] == "DMSO"){
    SNN_CellCycleDat$armSNN2[i] <- "DMSO"
  }
  SNN_CellCycleDat$log.G2M.Score[i] <- log1p(SNN_CellCycleDat$G2M.Score[i])
  SNN_CellCycleDat$log.S.Score[i] <- log1p(SNN_CellCycleDat$S.Score[i])
}


SNN_CellCycleDat <- SNN_CellCycleDat %>% group_by(armSNN)
write.csv(SNN_CellCycleDat, file = "05_output/SNN_CellCycleDat.csv", row.names = F)

SNN_0_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 0)
SNN_1_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 1)
SNN_2_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 2)
SNN_3_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 3)
SNN_4_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 4)
SNN_5_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 5)
SNN_6_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 6)
SNN_7_CellCycleDat <- subset(SNN_CellCycleDat, SNN_CellCycleDat$integrated_snn_res.0.25 == 7)


write.csv(SNN_0_CellCycleDat, file = "05_output/SNN_0_CellCycleDat.csv", row.names = F)
write.csv(SNN_1_CellCycleDat, file = "05_output/SNN_1_CellCycleDat.csv", row.names = F)
write.csv(SNN_2_CellCycleDat, file = "05_output/SNN_2_CellCycleDat.csv", row.names = F)
write.csv(SNN_3_CellCycleDat, file = "05_output/SNN_3_CellCycleDat.csv", row.names = F)
write.csv(SNN_4_CellCycleDat, file = "05_output/SNN_4_CellCycleDat.csv", row.names = F)
write.csv(SNN_5_CellCycleDat, file = "05_output/SNN_5_CellCycleDat.csv", row.names = F)
write.csv(SNN_6_CellCycleDat, file = "05_output/SNN_6_CellCycleDat.csv", row.names = F)
write.csv(SNN_7_CellCycleDat, file = "05_output/SNN_7_CellCycleDat.csv", row.names = F)

##########################################################################################################
##########################################################################################################
##########################################################################################################

# SNN Clustering Analysis

##########################################################################################################
##########################################################################################################
##########################################################################################################

int$armSNN <- paste0(int$Arm, "_", int$integrated_snn_res.0.25)
Idents(int) <- int$armSNN

armMarkers_0 <- FindMarkers(int, ident.1 = "NM002_0", ident.2 = "DMSO_0", test.use = "MAST", logfc.threshold = 0)
armMarkers_1 <- FindMarkers(int, ident.1 = "NM002_1", ident.2 = "DMSO_1", test.use = "MAST", logfc.threshold = 0)
armMarkers_2 <- FindMarkers(int, ident.1 = "NM002_2", ident.2 = "DMSO_2", test.use = "MAST", logfc.threshold = 0)
armMarkers_3 <- FindMarkers(int, ident.1 = "NM002_3", ident.2 = "DMSO_3", test.use = "MAST", logfc.threshold = 0)
armMarkers_4 <- FindMarkers(int, ident.1 = "NM002_4", ident.2 = "DMSO_4", test.use = "MAST", logfc.threshold = 0)
armMarkers_5 <- FindMarkers(int, ident.1 = "NM002_5", ident.2 = "DMSO_5", test.use = "MAST", logfc.threshold = 0)
armMarkers_6 <- FindMarkers(int, ident.1 = "NM002_6", ident.2 = "DMSO_6", test.use = "MAST", logfc.threshold = 0)
armMarkers_7 <- FindMarkers(int, ident.1 = "NM002_7", ident.2 = "DMSO_7", test.use = "MAST", logfc.threshold = 0)

armMarkers_0$cluster <- "0"
armMarkers_1$cluster <- "1"
armMarkers_2$cluster <- "2"
armMarkers_3$cluster <- "3"
armMarkers_4$cluster <- "4"
armMarkers_5$cluster <- "5"
armMarkers_6$cluster <- "6"
armMarkers_7$cluster <- "7"

allArmMarkersByCluster <- rbind(armMarkers_0, armMarkers_1, armMarkers_2, armMarkers_3, armMarkers_4, armMarkers_5, armMarkers_6, armMarkers_7)
write.csv(allArmMarkersByCluster, file = "05_output/allArmMarkersBySNNCluster.csv")

EV_0 <- EnhancedVolcano(armMarkers_0,
                         x = 'avg_log2FC',
                         y = "p_val_adj",
                         lab = rownames(armMarkers_0),
                         title = "UM002 vs DMSO in Cluster 0",
                         titleLabSize = 36,
                         pCutoff = 10e-9,
                         FCcutoff = 0.23,
                         col = c("black", "black", "black", "red3"), colAlpha = 1,
                         cutoffLineType = "twodash",
                         cutoffLineWidth = 0.8,
                         pointSize = 10,
                         legendPosition = 'none',
                         drawConnectors = T,
                         labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                         xlim = c(min(armMarkers_0[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                  max(armMarkers_0[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_0
EV_1 <- EnhancedVolcano(armMarkers_1,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_1),
                        title = "UM002 vs DMSO in Cluster 1",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_1[["avg_log2FC"]], na.rm=TRUE) - 0.5,
                                 max(armMarkers_1[["avg_log2FC"]], na.rm=TRUE) + 0.5)
) + coord_flip()
# EV_1
EV_2 <- EnhancedVolcano(armMarkers_2,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_2),
                        title = "UM002 vs DMSO in Cluster 2",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_2[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                 max(armMarkers_2[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_2
EV_3 <- EnhancedVolcano(armMarkers_3,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_3),
                        title = "UM002 vs DMSO in Cluster 3",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_3[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                 max(armMarkers_3[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_3
EV_4 <- EnhancedVolcano(armMarkers_4,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_4),
                        title = "UM002 vs DMSO in Cluster 4",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_4[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                 max(armMarkers_4[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_4
EV_5 <- EnhancedVolcano(armMarkers_5,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_5),
                        title = "UM002 vs DMSO in Cluster 5",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_5[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                 max(armMarkers_5[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_6
EV_6 <- EnhancedVolcano(armMarkers_6,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_6),
                        title = "UM002 vs DMSO in Cluster 6",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_6[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                 max(armMarkers_6[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_6
EV_7 <- EnhancedVolcano(armMarkers_7,
                        x = 'avg_log2FC',
                        y = "p_val_adj",
                        lab = rownames(armMarkers_7),
                        title = "UM002 vs DMSO in Cluster 7",
                        titleLabSize = 36,
                        pCutoff = 10e-9,
                        FCcutoff = 0.23,
                        col = c("black", "black", "black", "red3"), colAlpha = 1,
                        cutoffLineType = "twodash",
                        cutoffLineWidth = 0.8,
                        pointSize = 10,
                        legendPosition = 'none',
                        drawConnectors = T,
                        labSize = 12, labhjust = 1, boxedLabels = T, caption = NULL, labvjust = 2.5,
                        xlim = c(min(armMarkers_7[["avg_log2FC"]], na.rm=TRUE) - 0.35,
                                 max(armMarkers_7[["avg_log2FC"]], na.rm=TRUE) + 0.35)
) + coord_flip()
# EV_7

pdf(file = "05_output/snn_0.25_volcanoes_1.pdf", width = 35, height = 35*2)
plot_grid(plotlist = list(EV_0, EV_1, EV_2, EV_3),
          ncol = 1, nrow = 4,
          labels = "AUTO")
dev.off()
pdf(file = "05_output/snn_0.25_volcanoes_2.pdf", width = 35, height = 35*2)
plot_grid(plotlist = list(EV_4, EV_5, EV_6, EV_7),
          ncol = 1, nrow = 4,
          labels = "AUTO")
dev.off()

#### Pie charts of cell-cycle phase within SNN clusters.

#Cluster0
Cluster0sub <- subset(int, integrated_snn_res.0.25 == "0")
Cluster0subdf <- data.frame(table(Cluster0sub$CellCycleIdents))
Cluster0subdf$perc <- Cluster0subdf$Freq / sum(Cluster0subdf$Freq) * 100
Cluster0subdf
Cluster0tab <- Cluster0subdf$perc
names(Cluster0tab) <- Cluster0subdf$Var1
lbls <- paste(names(Cluster0tab), "\n", Cluster0tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster0.pdf")
pie3D(Cluster0tab, explode = 0.1, labels = lbls, main="Cluster 0")
dev.off()
rm(Cluster0sub)
#Cluster1
Cluster1sub <- subset(int, integrated_snn_res.0.25 == "1")
Cluster1subdf <- data.frame(table(Cluster1sub$CellCycleIdents))
Cluster1subdf$perc <- Cluster1subdf$Freq / sum(Cluster1subdf$Freq) * 100
Cluster1subdf
Cluster1tab <- Cluster1subdf$perc
names(Cluster1tab) <- Cluster1subdf$Var1
lbls <- paste(names(Cluster1tab), "\n", Cluster1tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster1.pdf")
pie3D(Cluster1tab, explode = 0.1, labels = lbls, main="Cluster 1")
dev.off()
rm(Cluster1sub)
#Cluster2
Cluster2sub <- subset(int, integrated_snn_res.0.25 == "2")
Cluster2subdf <- data.frame(table(Cluster2sub$CellCycleIdents))
Cluster2subdf$perc <- Cluster2subdf$Freq / sum(Cluster2subdf$Freq) * 100
Cluster2subdf
Cluster2tab <- Cluster2subdf$perc
names(Cluster2tab) <- Cluster2subdf$Var1
lbls <- paste(names(Cluster2tab), "\n", Cluster2tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster2.pdf")
pie3D(Cluster2tab, explode = 0.1, labels = lbls, main="Cluster 2")
dev.off()
rm(Cluster2sub)
#Cluster3
Cluster3sub <- subset(int, integrated_snn_res.0.25 == "3")
Cluster3subdf <- data.frame(table(Cluster3sub$CellCycleIdents))
Cluster3subdf$perc <- Cluster3subdf$Freq / sum(Cluster3subdf$Freq) * 100
Cluster3subdf
Cluster3tab <- Cluster3subdf$perc
names(Cluster3tab) <- Cluster3subdf$Var1
lbls <- paste(names(Cluster3tab), "\n", Cluster3tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster3.pdf")
pie3D(Cluster3tab, explode = 0.1, labels = lbls, main="Cluster 3")
dev.off()
rm(Cluster3sub)
#Cluster4
Cluster4sub <- subset(int, integrated_snn_res.0.25 == "4")
Cluster4subdf <- data.frame(table(Cluster4sub$CellCycleIdents))
Cluster4subdf$perc <- Cluster4subdf$Freq / sum(Cluster4subdf$Freq) * 100
Cluster4subdf
Cluster4tab <- Cluster4subdf$perc
names(Cluster4tab) <- Cluster4subdf$Var1
lbls <- paste(names(Cluster4tab), "\n", Cluster4tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster4.pdf")
pie3D(Cluster4tab, explode = 0.05, labels = lbls, main="Cluster 4")
dev.off()
rm(Cluster4sub)
#Cluster5
Cluster5sub <- subset(int, integrated_snn_res.0.25 == "5")
Cluster5subdf <- data.frame(table(Cluster5sub$CellCycleIdents))
Cluster5subdf$perc <- Cluster5subdf$Freq / sum(Cluster5subdf$Freq) * 100
Cluster5subdf
Cluster5tab <- Cluster5subdf$perc
names(Cluster5tab) <- Cluster5subdf$Var1
lbls <- paste(names(Cluster5tab), "\n", Cluster5tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster5.pdf")
pie3D(Cluster5tab, explode = 0.1, labels = lbls, main="Cluster 5")
dev.off()
rm(Cluster5sub)
#Cluster6
Cluster6sub <- subset(int, integrated_snn_res.0.25 == "6")
Cluster6subdf <- data.frame(table(Cluster6sub$CellCycleIdents))
Cluster6subdf$perc <- Cluster6subdf$Freq / sum(Cluster6subdf$Freq) * 100
Cluster6subdf
Cluster6tab <- Cluster6subdf$perc
names(Cluster6tab) <- Cluster6subdf$Var1
lbls <- paste(names(Cluster6tab), "\n", Cluster6tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster6.pdf")
pie3D(Cluster6tab, explode = 0.1, labels = lbls, main="Cluster 6")
dev.off()
rm(Cluster6sub)
#Cluster7
Cluster7sub <- subset(int, integrated_snn_res.0.25 == "7")
Cluster7subdf <- data.frame(table(Cluster7sub$CellCycleIdents))
Cluster7subdf$perc <- Cluster7subdf$Freq / sum(Cluster7subdf$Freq) * 100
Cluster7subdf
Cluster7tab <- Cluster7subdf$perc
names(Cluster7tab) <- Cluster7subdf$Var1
lbls <- paste(names(Cluster7tab), "\n", Cluster7tab, sep = "")
pdf(file = "05_output/stateCellCyclePie_1_Cluster7.pdf")
pie3D(Cluster7tab, explode = 0.1, labels = lbls, main="Cluster 7")
dev.off()
rm(Cluster7sub)

##################################################################################
##################################################################################
##################################################################################
