library(Seurat)
library(CytoTRACE2)
library(reticulate)
library(anndata)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

use_condaenv('/home/liyr/miniconda3/envs/R4.3.3')

qc <- data.frame(matrix(nrow = 0,ncol = 7))
colnames(qc) <- c('uuid','n_gene','n_cell',
                  'n_gene_filter','n_cell_filter',
                  "nFeature_RNA_median_filter","nCount_RNA_median_filter")


files = list.files('/home/liyr/HuBMAP/Data/Expression_batch_download')

files_complete =list.files('/home/liyr/HuBMAP/Data/h5ad/')
files_complete = gsub('.h5ad','',files_complete)
files = setdiff(files,files_complete)


for (uuid in files){
  
  print(uuid)
  ad_file = paste0("/home/liyr/HuBMAP/Data/Expression_batch_download/",uuid,"/expr.h5ad")
  
  if (file.exists(ad_file)){

    ad <- read_h5ad(ad_file)
    # exp =  as.data.frame(ad$X) # cell X gene
    exp =  as.matrix(ad$X)
    
    if (ncol(ad$var)!=0){
      gene  = ad$var %>% as.data.frame()
      gene$ENSEMBL <- rownames(gene)
      gene <- gene[!is.na(gene$hugo_symbol) & (!duplicated(gene$hugo_symbol)),]
      exp <- exp[1:nrow(exp),gene$ENSEMBL]
      colnames(exp) <- gene$hugo_symbol

      gene_name = gene
      gene_name$SYMBOL = gene_name$hugo_symbol
      gene_name$gene = gene_name$ENSEMBL

    }else{
      gene = colnames(exp) %>% as.data.frame()
      colnames(gene) <- 'gene'
      gene = separate(gene,col="gene",sep= "[.]", into=c("ENSEMBL","Sup"),remove=FALSE)
      
      gene_symbol = bitr(gene$ENSEMBL,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb='org.Hs.eg.db')
      gene_name = merge(gene,gene_symbol,col='ENSEMBL',all.x=T)
      gene_name = gene_name[!duplicated(gene_name$ENSEMBL),]
      rownames(gene_name) = gene_name$gene
      gene_name <- gene_name[!duplicated(gene_name$SYMBOL) & (!is.na(gene_name$SYMBOL)),]
      exp <- exp[1:nrow(exp),gene_name$gene]
      colnames(exp) <- gene_name$SYMBOL
    }
    
    exp <- t(exp) %>% as.data.frame() # gene X cell
    
    seu <- CreateSeuratObject(exp, min.cells = 3, min.features = 200) # filter cells and genes
    # VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA")) # too few features, poor quality
    
    qc[nrow(qc) + 1,] = c(uuid,dim(exp),dim(seu),median(seu$nFeature_RNA),median(seu$nCount_RNA))
    fwrite(qc,'/home/liyr/HuBMAP/Data/QC.csv')
    
    # saveRDS(seu,
    #         paste0('/home/liyr/HuBMAP/RDS/',uuid,'.rds'))
    

    # Run cytotrace2
    cytotrace2_result <- cytotrace2(seu,species = "human",is_seurat=T,seed = 2024,ncores = 26)
    
    saveRDS(cytotrace2_result,
            paste0('/home/liyr/HuBMAP/Data/rds/',uuid,'.rds'))
    
    # Save filtered anndata
    rownames(gene_name) = gene_name$SYMBOL
    gene_name = gene_name[rownames(cytotrace2_result),] # same genes as filtered
    gene_name = gene_name[gene_name$gene %in% ad$var_names,]

    ad2 = ad[colnames(cytotrace2_result),gene_name$gene] # filter gene & cell

    ad2$var['SYMBOL'] = gene_name$SYMBOL
    ad2$var['ENSEMBL'] = gene_name$ENSEMBL
    ad2$var_names = gene_name$SYMBOL 
    
    ad2$obs[colnames(cytotrace2_result@meta.data)] = cytotrace2_result@meta.data[colnames(cytotrace2_result@meta.data)]
    write_h5ad(ad2,paste0('/home/liyr/HuBMAP/Data/h5ad/',uuid,'.h5ad'))

    # rm(cytotrace2_result)
    # rm(seu)
    # rm(exp)
    # rm(ad)
    # rm(ad2)
    
    print('---------------------------------------------')
    
  }
  
}