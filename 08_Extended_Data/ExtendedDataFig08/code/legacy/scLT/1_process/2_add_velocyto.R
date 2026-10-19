library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(data.table)

LSK <- readRDS("/home/liyaru/DATA/8_transition/1_DATA/4_sup/FigShare/lsk_rna_doub_remv2.Rds")
#DefaultAssay(LSK) <- "integrated"
DimPlot(LSK,reduction = "cce",group.by = "cell_type")

 # read sup name
meta = LSK@meta.data
meta$cellname = rownames(meta)
meta = separate(meta,col = "cellname",sep="-",into = c("name1","name2","name3","name4"),remove = F)

supname = unique(meta[,c("name1","name3","name4","sample")])

filename = c("d2_5_possorted_genome_bam_QVXBK.loom",
  "d5_1_possorted_genome_bam_BSO8O.loom",
  "d5_2_possorted_genome_bam_AHZRZ.loom",
  "d5_3_possorted_genome_bam_NYHEL.loom",
  "d5_4_possorted_genome_bam_NBOY3.loom",
  "d5_5_possorted_genome_bam_F2EGD.loom",
  "d5_6_possorted_genome_bam_RY6HN.loom",
  "d5_7_possorted_genome_bam_C1ZG8.loom",
  "d5_8_possorted_genome_bam_QLDHR.loom")
supname$file = filename

# build velocyto object  in seurat
result <- list()
path = "/home/liyaru/DATA/8_transition/1_DATA/4_NBT_2023/GSE216606_LSK_scRNA/bam/velocyto/"
for (i in 1:nrow(supname)){
  
  #i = 1
  filepath = paste0(path,supname[i,'file'])
  ldat <- ReadVelocity(filepath)
  bm <- as.Seurat(x = ldat)
  
  sample = supname[i,'sample']
  sample_num = supname[i,'name3']
  
  m = bm@meta.data
  m$cellname = rownames(m)
  m = separate(m,col="cellname",sep = ":",into=c("cellname2","barcode"),remove = F)
  m = separate(m,col="cellname2",sep ="_possorted_genome_bam_",into=c("sample","supname"))
  cellname_new = paste0(sample,"-",m$barcode,"-",sample_num,"-rna")
  
  colnames(bm@assays$spliced@counts) <- cellname_new
  colnames(bm@assays$spliced@data) <- cellname_new

  colnames(bm@assays$unspliced@counts) <- cellname_new
  colnames(bm@assays$unspliced@data) <- cellname_new

  colnames(bm@assays$ambiguous@counts) <- cellname_new
  colnames(bm@assays$ambiguous@data) <- cellname_new

  rownames(bm@meta.data) <- cellname_new
  Idents(bm) <- cellname_new
  
  cells = intersect(cellname_new,rownames(meta))
  bm <- subset(bm,idents = cells)

  result <- c(result,bm)
}

saveRDS(result,"/home/liyaru/DATA/8_transition/1_DATA/4_NBT_2023/velocyto_list.rds")

velo <- merge(result[[1]],
              y=result[2:length(result)],
              project="velocity")

saveRDS(velo,"/home/liyaru/DATA/8_transition/1_DATA/4_NBT_2023/velocyto_merge.rds")


velo = readRDS("/home/liyaru/DATA/8_transition/1_DATA/4_NBT_2023/velocyto_merge.rds")

DefaultAssay(LSK) <- 'RNA'
colname_intersect <- intersect(colnames(LSK),colnames(velo))
rowname_intersect <- intersect(rownames(LSK),rownames(velo))
velo2  <- subset(velo,cells=colname_intersect,features=rowname_intersect)
LSK <- subset(LSK,cells=colname_intersect,features=rowname_intersect)

LSK[["spliced"]] <- velo2[["spliced"]]
LSK[["unspliced"]] <- velo2[["unspliced"]]
LSK[["ambiguous"]] <- velo2[["ambiguous"]]

LSK@assays$RNA@meta.features[VariableFeatures(LSK), 'highly_variable_genes'] <- 1
LSK@assays$RNA@meta.features[is.na(LSK@assays$RNA@meta.features$highly_variable_genes), 'highly_variable_genes'] <- 0

source("/media/liyaru/LYR/radiation/R/50_addvelocyto.R")

DefaultAssay(LSK) <- "RNA"
seurat2anndata(obj = LSK,outFile = "/home/liyaru/DATA/8_transition/1_DATA/4_NBT_2023/LSK_velocity_RNAassay.h5ad")


c = fread("/media/liyaru/LYR/Diff_change/7_lineage_tracing_multitag/CellTag-multi-2023-main/clone_tables/hsc.rna&atac.r1&2_master_v2.csv") %>% as.data.frame()
c = c[c$assay == "rna",]
c_t1 = c[c$day == "d2",]
length(unique(c_t1$clone.id))

c_t2 = c[c$day == "d5",]
length(unique(c_t2$clone.id))


















