library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(data.table)
library(SeuratData)
library(SeuratDisk)

# young2 all
redeem = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/Young2.All.Seurat.RDS")


meta = redeem@meta.data
print(redeem)
DimPlot(redeem,group.by = "STD.CellType",reduction = "wnn.umap",label = T,label.size = 5,split.by = "Sample")
DimPlot(redeem,group.by = "STD.CellType",reduction = "umap.rna",label = T,label.size = 5,split.by = "Sample")

table(redeem$Sample)
table(redeem$STD.CellType)
table(redeem$ClonalGroup)

table(redeem$STD.CellType,redeem$ClonalGroup)

cells = meta[meta$Sample == 'DN9_BMMC' & !is.na(meta$Sample),]
BMMC = subset(redeem,cells = rownames(cells))

# add velocyto
# young2  sample7 (BMMC)
ldat <- ReadVelocity("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/fastq_young2/young2_7human_cellranger/velocyto/young2_7human_cellranger.loom")
bm <- as.Seurat(x = ldat)

m = bm@meta.data
m$cellname = rownames(m)
m = separate(m,col="cellname",sep = ":",into=c("cellname2","barcode"),remove = F)
m = separate(m,col="cellname2",sep ="_possorted_genome_bam_",into=c("sample","supname"))
m$barcode_new = gsub("x","",m$barcode)
cellname_new = paste0(m$barcode_new,"-1")

colnames(bm@assays$spliced@counts) <- cellname_new
colnames(bm@assays$spliced@data) <- cellname_new

colnames(bm@assays$unspliced@counts) <- cellname_new
colnames(bm@assays$unspliced@data) <- cellname_new

colnames(bm@assays$ambiguous@counts) <- cellname_new
colnames(bm@assays$ambiguous@data) <- cellname_new

rownames(bm@meta.data) <- cellname_new
Idents(bm) <- cellname_new

cells_overlap = intersect(cellname_new,rownames(cells))
bm <- subset(bm,idents = cells_overlap)

# to h5ad
redeem = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/Young2.All.Seurat.RDS")
setwd("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/")
DefaultAssay(redeem) <- 'ATAC'
redeem[["SCT"]] <- NULL
redeem[["RNA"]] <- NULL
SaveH5Seurat(redeem, filename = "redeem.h5Seurat")
Convert("redeem.h5Seurat", dest = "h5ad")


# get neighbors
redeem = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/Young2.All.Seurat.RDS")


meta = redeem@meta.data
cells = meta[meta$Sample == 'DN9_HSPC' & !is.na(meta$Sample),]
# 
# # difference in anndata
# cc = rownames(cells)
# cc <- cc[cc != "CGGACCTAGTTATCCT-2"]
 
redeem = subset(redeem,cells = rownames(cells))

cells_new <- gsub("-2","",rownames(redeem@meta.data) )
redeem <- RenameCells(redeem, new.names = cells_new)

redeem@meta.data

cells <- read.csv("/media/liyaru/LYR/Diff_change/11_redeem/filtered_cells.txt",header = F)
redeem <- subset(redeem,cells = cells$V1)

redeem <- FindMultiModalNeighbors(redeem, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 50)

redeem@neighbors
nn_idx <- redeem@neighbors$weighted.nn@nn.idx
nn_dist <- redeem@neighbors$weighted.nn@nn.dist
nn_cells <- redeem@neighbors$weighted.nn@cell.names

setwd("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/")
# save neighborhood graph
write.table(nn_idx, "nn_idx.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, "nn_dist.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, "nn_cells.txt", sep = ',', row.names = F, col.names = F, quote = F)

