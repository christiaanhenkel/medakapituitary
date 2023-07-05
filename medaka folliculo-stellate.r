library(RColorBrewer)
library(rcartocolor)
library(Seurat)

sccolp <- colorRampPalette(c(carto_pal(7, "Sunset")))		# gene expression colour scale
sccol <- sccolp(41)
clcol <- c(carto_pal(12, "Vivid")[1:10], brewer.pal(3, "Set1")[1], carto_pal(12, "Pastel")[1:5])

zscale <- function(y) { ifelse(y > 0, y/max(y), ifelse(y < 0, -y/min(y), 0)) }
scplot <- function(x, umap1, umap2, xcol, diverging, pch, cex, add = F) { 
	y <- order(abs(zscale(x)))
	if (diverging == T) {
		ycol = xcol[ceiling(0.00001 + (zscale(x[y]) + 1) * (length(xcol) - 1) * 0.5)]
		} else {
		ycol = xcol[1 + ceiling((length(xcol) - 1) * ((x[y] - min(x)) / max(x)))]
		}
	if (add == T) {
		points(umap1[y], umap2[y], col = ycol, pch = pch, cex = cex)
		} else {
		plot(umap1[y], umap2[y], col = ycol, pch = pch, cex = cex)
		}
	}

# Data from https://identifiers.org/geo:GSE162787 	
fdata <- Read10X(data.dir = "../Female_data/")
fclusters <- read.delim("../Female_metadata_2592x9.tsv", row.names = 1)
fclusters <- fclusters[colnames(fdata),]	# same order
fcells <- CreateSeuratObject(counts = fdata, metadata = fclusters, assay = "female")
fcells@meta.data$sample <- fclusters$sample		# Not managing to import the annotations any other way
#fcells@meta.data$celltype <- fclusters$celltype
fcells <- NormalizeData(fcells, normalization.method = "LogNormalize", scale.factor = 10000)
fcells <- FindVariableFeatures(fcells, selection.method = "vst", nfeatures = 2000)

mdata <- Read10X(data.dir = "../Male_data/")
mclusters <- read.delim("../Male_metadata_3804x9.tsv", row.names = 1)
mclusters <- mclusters[colnames(mdata),]	# same order
mcells <- CreateSeuratObject(counts = mdata, metadata = mclusters, assay = "male")
mcells@meta.data$sample <- mclusters$sample		# Not managing to import the annotations any other way
#mcells@meta.data$celltype <- mclusters$celltype
mcells <- NormalizeData(mcells, normalization.method = "LogNormalize", scale.factor = 10000)
mcells <- FindVariableFeatures(mcells, selection.method = "vst", nfeatures = 2000)

# Integration
anchors <- FindIntegrationAnchors(object.list = c(fcells, mcells), dims = 1:15)
allcells <- IntegrateData(anchorset = anchors, dims = 1:15)
allcells <- ScaleData(allcells)
allcells <- RunPCA(allcells)
allcells <- RunUMAP(allcells, dims = 1:15)
DimPlot(allcells, reduction = "umap")

allcells <- FindNeighbors(allcells, dims = 1:15)
allcells <- FindClusters(allcells, resolution = 0.5)
#allcells@meta.data$seurat_clusters <- as.factor(allexpr$cell_cluster.x)	# replace by previously defined clusters
#allcells@meta.data$integrated_snn_res.0.5 <- as.factor(allexpr$cell_cluster.x)
#allcells@active.ident <- as.factor(allexpr$cell_cluster.x)
#names(allcells@active.ident) <- allexpr$Row.names.x

allsummary <- data.frame(allcells@meta.data, allcells@reductions$umap@cell.embeddings)
fclusters <- cbind(fclusters, paste(rownames(fclusters), rep(c("_1")), sep = ""))
colnames(fclusters)[9] <- c("barcode")
mclusters <- cbind(mclusters, paste(rownames(mclusters), rep(c("_2")), sep = ""))
colnames(mclusters)[9] <- c("barcode")
allsummary <- rbind(merge(allsummary, fclusters, by.x = 0, by.y = "barcode", all = F), merge(allsummary, mclusters, by.x = 0, by.y = "barcode", all = F))
rownames(allsummary) <- allsummary$Row.names
allcells@active.ident <- as.factor(allsummary$cell_cluster)
allcells@meta.data$seurat_clusters <- as.factor(allsummary$cell_cluster)

fexpr <- t(as.matrix(fcells@assays$female@data))
fcounts <- t(as.matrix(fdata))
mexpr <- t(as.matrix(mcells@assays$male@data))
mcounts <- t(as.matrix(mdata))

fexpr <- fexpr[,intersect(colnames(fexpr), colnames(mexpr))]
fcounts <- fcounts[,intersect(colnames(fexpr), colnames(mexpr))]
fexpr <- merge(fclusters, fexpr, by.x = 0, by.y = 0)
fcounts <- merge(fclusters, fcounts, by.x = 0, by.y = 0)
mexpr <- mexpr[,intersect(colnames(fexpr), colnames(mexpr))]	
mcounts <- mcounts[,intersect(colnames(fexpr), colnames(mexpr))]
mexpr <- merge(mclusters, mexpr, by.x = 0, by.y = 0)
mcounts <- merge(mclusters, mcounts, by.x = 0, by.y = 0)

allexpr <- rbind(merge(allsummary, fexpr, by.x = 0, by.y = "barcode", all = F), merge(allsummary, mexpr, by.x = 0, by.y = "barcode", all = F))
rownames(allexpr) <- allexpr$Row.names

# Some UMAP plots
pdf("original_clusters_umap__01072022.pdf", useDingbats = F, width = 10, height = 10)
	allvis <- allexpr[order(allexpr$nFeature_RNA.x),]
	plot(allvis$UMAP_1.x, allvis$UMAP_2.x, col = clcol[allvis$cell_cluster.x], pch = 16, cex = 0.75)
dev.off()

pdf("tshr_ENSORLG00000014222_03042023.pdf", useDingbats = F, width = 10, height = 10)
	scplot(allexpr$ENSORLG00000014222, allexpr$UMAP_1.x, allexpr$UMAP_2.x, sccol, F, 16, 0.75)
dev.off()

pdf("cyp19b_ENSORLG00000005548_03042023.pdf", useDingbats = F, width = 10, height = 10)
	scplot(allexpr$ENSORLG00000005548, allexpr$UMAP_1.x, allexpr$UMAP_2.x, sccol, F, 16, 0.75)
dev.off()


# Marker genes
cl13markers <- FindMarkers(allcells, ident.1 = c(13), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 71
cl14markers <- FindMarkers(allcells, ident.1 = c(14), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 107 
cl15markers <- FindMarkers(allcells, ident.1 = c(15), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 41
cl16markers <- FindMarkers(allcells, ident.1 = c(16), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 101
cl13markers <- cl13markers[cl13markers$p_val_adj<0.05,] # 68 left
cl14markers <- cl14markers[cl14markers$p_val_adj<0.05,] # 107 left
cl15markers <- cl15markers[cl15markers$p_val_adj<0.05,] # 41 left
cl16markers <- cl16markers[cl16markers$p_val_adj<0.05,] # 90 left
fsmarkers1 <- FindMarkers(allcells, ident.1 = c(13, 14, 15, 16), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 104
fsmarkers2 <- FindMarkers(allcells, ident.1 = c(13, 14, 16), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 133
fsmarkers3 <- FindMarkers(allcells, ident.1 = c(13, 14, 16), ident.2 = c(1,2,3,5,6,7,8,9,10,11,12), only.pos = T, min.cells.group = 2, logfc.threshold = 1) # 141

litmarkers <- read.delim("fs_markergenes_literature_april_2023_edit.txt")
annot <- read.delim("../../annotation/Medaka_annotation_check.txt", header = T, sep = " ") # Gene_ID	Gene_name

maxexpr <- max(allexpr[, 28:17581]) # 9.087876
pdf("fs_pituicyte_markes_literature_03042023.pdf", width = 10, height = 20, useDingbats = F)
plot(-1, -1, col = "white", xlim = c(1,31), ylim = c(60, 0))
for (g in 1:55) {
	text(17, g, litmarkers[g,1], cex = 0.5, pos = 4)
	text(20, g, litmarkers[g,2], cex = 0.5, pos = 4)
	text(25, g, litmarkers[g,3], cex = 0.5, pos = 4)
	for (cl in 1:16) {
		expression <- expm1(allexpr[allexpr$cell_cluster.x == cl ,colnames(allexpr) == litmarkers[g,2]])
		expressed <- length(expression[expression > 0]) / length(expression)
		points(cl, g, pch = 16, cex = expressed * 4, col = sccol[ceiling(1 + 40/maxexpr * log1p(mean(expression[expression > 0])))]) 
	}
}

pdf("cluster13_markers_03042023.pdf", width = 10, height = 25, useDingbats = F)
plot(-1, -1, col = "white", xlim = c(1,31), ylim = c(75, 0))
for (g in 1:68) {
	text(17, g, annot$Gene_name[annot$Gene_ID==rownames(cl13markers[g,])], cex = 0.5, pos = 4)
	text(20, g, rownames(cl13markers[g,]), cex = 0.5, pos = 4)
	text(25, g, annot$Description[annot$Gene_ID==rownames(cl13markers[g,])], cex = 0.5, pos = 4)
	for (cl in 1:16) {
		expression <- expm1(allexpr[allexpr$cell_cluster.x == cl ,colnames(allexpr) == rownames(cl13markers[g,])])
		expressed <- length(expression[expression > 0]) / length(expression)
		points(cl, g, pch = 16, cex = expressed * 4, col = sccol[ceiling(1 + 40/maxexpr * log1p(mean(expression[expression > 0])))]) 
	}
}
for (i in 1:16) { text(i, 0, i, cex = 0.5) }
dev.off()

pdf("cluster14_markers_03042023.pdf", width = 10, height = 40, useDingbats = F)
plot(-1, -1, col = "white", xlim = c(1,31), ylim = c(120, 0))
for (g in 1:107) {
	text(17, g, annot$Gene_name[annot$Gene_ID==rownames(cl14markers[g,])], cex = 0.5, pos = 4)
	text(20, g, rownames(cl14markers[g,]), cex = 0.5, pos = 4)
	text(25, g, annot$Description[annot$Gene_ID==rownames(cl14markers[g,])], cex = 0.5, pos = 4)
	for (cl in 1:16) {
		expression <- expm1(allexpr[allexpr$cell_cluster.x == cl ,colnames(allexpr) == rownames(cl14markers[g,])])
		expressed <- length(expression[expression > 0]) / length(expression)
		points(cl, g, pch = 16, cex = expressed * 4, col = sccol[ceiling(1 + 40/maxexpr * log1p(mean(expression[expression > 0])))]) 
	}
}
for (i in 1:16) { text(i, 0, i, cex = 0.5) }
dev.off()

pdf("cluster15_markers_03042023.pdf", width = 10, height = 15, useDingbats = F)
plot(-1, -1, col = "white", xlim = c(1,31), ylim = c(45, 0))
for (g in 1:90) {
	text(17, g, annot$Gene_name[annot$Gene_ID==rownames(cl15markers[g,])], cex = 0.5, pos = 4)
	text(20, g, rownames(cl15markers[g,]), cex = 0.5, pos = 4)
	text(25, g, annot$Description[annot$Gene_ID==rownames(cl15markers[g,])], cex = 0.5, pos = 4)
	for (cl in 1:16) {
		expression <- expm1(allexpr[allexpr$cell_cluster.x == cl ,colnames(allexpr) == rownames(cl15markers[g,])])
		expressed <- length(expression[expression > 0]) / length(expression)
		points(cl, g, pch = 16, cex = expressed * 4, col = sccol[ceiling(1 + 40/maxexpr * log1p(mean(expression[expression > 0])))]) 
	}
}
for (i in 1:16) { text(i, 0, i, cex = 0.5) }
dev.off()

pdf("cluster16_markers_03042023.pdf", width = 10, height = 30, useDingbats = F)
plot(-1, -1, col = "white", xlim = c(1,31), ylim = c(90, 0))
for (g in 1:90) {
	text(17, g, annot$Gene_name[annot$Gene_ID==rownames(cl16markers[g,])], cex = 0.5, pos = 4)
	text(20, g, rownames(cl16markers[g,]), cex = 0.5, pos = 4)
	text(25, g, annot$Description[annot$Gene_ID==rownames(cl16markers[g,])], cex = 0.5, pos = 4)
	for (cl in 1:16) {
		expression <- expm1(allexpr[allexpr$cell_cluster.x == cl ,colnames(allexpr) == rownames(cl16markers[g,])])
		expressed <- length(expression[expression > 0]) / length(expression)
		points(cl, g, pch = 16, cex = expressed * 4, col = sccol[ceiling(1 + 40/maxexpr * log1p(mean(expression[expression > 0])))]) 
	}
}
for (i in 1:16) { text(i, 0, i, cex = 0.5) }
dev.off()
