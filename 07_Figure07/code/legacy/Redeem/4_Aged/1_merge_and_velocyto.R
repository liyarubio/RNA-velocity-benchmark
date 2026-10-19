library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(data.table)

setwd("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/young2")

young2 = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/Young2.All.Seurat.RDS")
DimPlot(young2)
m = young2@meta.data
mm = m[m$Sample == 'DN9_HSPC',]
young2_HSPC = subset(young2,cells = rownames(mm))

DefaultAssay(young2_HSPC) <- 'RNA'
young2_HSPC <- RunUMAP(young2_HSPC,dims = 1:20)
DimPlot(young2_HSPC,reduction = "umap",label = T)

# check UMAP
# DefaultAssay(young2_HSPC) <- 'RNA'
# young2_HSPC <- FindVariableFeatures(young2_HSPC)
# young2_HSPC <- NormalizeData(young2_HSPC)
# young2_HSPC <- ScaleData(young2_HSPC)
# young2_HSPC <- RunPCA(young2_HSPC)
# young2_HSPC <- FindNeighbors(young2_HSPC, dims = 1:20)
# young2_HSPC <- RunUMAP(young2_HSPC,dims = 1:20)
# DimPlot(young2_HSPC,label = T)
# Idents(young2_HSPC)
# m = young2_HSPC@meta.data

aged2 = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/23290004/Old2.BMMC_HSPC.Seurat.RDS")
DimPlot(aged2)
m = aged2@meta.data
mm = m[m$Sample == 'Old2_HSPC',]
aged2_HSPC = subset(aged2,cells = rownames(mm))

DefaultAssay(aged2_HSPC) <- 'RNA'
aged2_HSPC <- RunUMAP(aged2_HSPC,dims = 1:20)
DimPlot(aged2_HSPC,reduction = "umap",label = T)

HSPC_list <- list(young2_HSPC,aged2_HSPC)
anchors <- FindIntegrationAnchors(object.list = HSPC_list,dims=1:20)
merge_inte <- IntegrateData(anchorset = anchors, dims = 1:20)

merge_inte

DefaultAssay(merge_inte) <- 'integrated'
merge_inte <- ScaleData(merge_inte)
merge_inte <- RunPCA(merge_inte)
merge_inte <- FindNeighbors(merge_inte, dims = 1:20)
merge_inte <- RunUMAP(merge_inte,dims = 1:20)

DimPlot(merge_inte,label = T,split.by = "Sample")
merge_inte$Celltype <- Idents(merge_inte)

# add velocyto
meta = merge_inte@meta.data
meta = meta[!is.na(meta$Sample),]
meta$cell = rownames(meta)
meta = separate(meta,col = "cell",into = c("barcode","sup"),remove = F,sep = "-")
sup = meta[,c("Sample","sup")] %>% unique()


files = c("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/young2/young2_8_cellranger/velocyto/young2_8_cellranger.loom", #young2 HSPC
          "/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/aged2/aged2_HSPC_cellranger/velocyto/aged2_HSPC_cellranger.loom" ) # aged2_HSPC

sup_names = c("-2_1","-2_2")
result <- list()

for (i in 1:length(files)){
  
  f = files[i]
  print(f)
  
  s = sup_names[i]
  
  ldat <- ReadVelocity(f)
  bm <- as.Seurat(x = ldat)
  
  m = bm@meta.data
  m$cellname = rownames(m)
  m = separate(m,col="cellname",sep = ":",into=c("cellname2","barcode"),remove = F)
  m = separate(m,col="cellname2",sep ="_possorted_genome_bam_",into=c("sample","supname"))
  m$barcode_new = gsub("x","",m$barcode)
  cellname_new = paste0(m$barcode_new,s)
  
  colnames(bm@assays$spliced@counts) <- cellname_new
  colnames(bm@assays$spliced@data) <- cellname_new
  
  colnames(bm@assays$unspliced@counts) <- cellname_new
  colnames(bm@assays$unspliced@data) <- cellname_new
  
  colnames(bm@assays$ambiguous@counts) <- cellname_new
  colnames(bm@assays$ambiguous@data) <- cellname_new
  
  rownames(bm@meta.data) <- cellname_new
  Idents(bm) <- cellname_new
  
  cells_overlap = intersect(cellname_new,rownames(meta))
  bm <- subset(bm,idents = cells_overlap)
  
  result <- c(result,bm)
  
}

velo <- merge(result[[1]],
              y=result[[2]],
              project="velocity")
velo

# intersect
DefaultAssay(merge_inte) <- 'RNA'
colname_intersect <- intersect(colnames(merge_inte),colnames(velo))
rowname_intersect <- intersect(rownames(merge_inte),rownames(velo))
velo2  <- subset(velo,cells=colname_intersect,features=rowname_intersect)
merge_inte <- subset(merge_inte,cells=colname_intersect,features=rowname_intersect)

merge_inte[["spliced"]] <- velo2[["spliced"]]
merge_inte[["unspliced"]] <- velo2[["unspliced"]]
merge_inte[["ambiguous"]] <- velo2[["ambiguous"]]

merge_inte@assays$RNA@meta.features[VariableFeatures(merge_inte), 'highly_variable_genes'] <- 1
merge_inte@assays$RNA@meta.features[is.na(merge_inte@assays$RNA@meta.features$highly_variable_genes), 'highly_variable_genes'] <- 0

source("/media/liyaru/LYR/radiation/R/50_addvelocyto.R")

DefaultAssay(merge_inte) <- "RNA"
saveRDS(merge_inte,"/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/young2_aged2_HSPC.rds")
seurat2anndata(obj = merge_inte,outFile = "/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/young2_aged2_HSPC.h5ad")


