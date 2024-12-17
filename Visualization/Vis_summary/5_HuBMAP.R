library(ggsci)
library(Seurat)
library(tibble)
library(stringr)
library(pheatmap)
library(tidyr)

check_method_name <- function(metrics){
  metrics <- metrics %>%
    filter(Method != "Dynamo_deterministic") %>%
    mutate(Method = gsub("Velocyto", "velocyto", Method,ignore.case = T))  %>%
    mutate(Method = gsub("scvelo_dynamical|scvelo_dynamic" , "scVelo(dyn)", Method,ignore.case = T))  %>%
    mutate(Method = gsub("scvelo_stochastic|scvelo$" , "scVelo(sto)", Method,ignore.case = T))  %>%
    mutate(Method = gsub("veloAE" , "VeloAE", Method,ignore.case = T))  %>%
    # mutate(Method = gsub("Dynamo_deterministic", "dynamo", Method))  %>%
    mutate(Method = gsub("Dynamo_stochastic|Dynamo", "dynamo", Method,ignore.case = T))  %>%
    mutate(Method = gsub("veloVAE", "VeloVAE", Method,ignore.case = T))  %>%
    mutate(Method = gsub("UniTvelo|UniteVelo", "UniTVelo", Method,ignore.case = T))  %>%
    mutate(Method = gsub("DeepVelo_SA|DeepveloSA", "DeepVelo(VAE)", Method,ignore.case = T)) %>%
    mutate(Method = gsub("Celldancer", "cellDancer", Method,ignore.case = T)) %>%
    mutate(Method = gsub("Latentvelo", "LatentVelo", Method,ignore.case = T)) %>%
    mutate(Method = gsub("DeepVelo_GB|DeepveloGB", "DeepVelo(GCN)", Method,ignore.case = T))%>%
    mutate(Method = gsub("stt", "STT", Method,ignore.case = T)) %>%
    mutate(Method = gsub("multivelo", "MultiVelo", Method,ignore.case = T))%>%
    mutate(Method = gsub("phylovelo", "PhyloVelo", Method,ignore.case = T)) %>%
    mutate(Method = gsub("tfvelo", "TFvelo", Method,ignore.case = T))
  
  metrics <- metrics[!metrics$Method %in% c('dynamo_deter','silhouette'),]
  # metrics$Method <- factor(metrics$Method,levels = names(method_palette))
  # metrics <- metrics[order(metrics$Method),]

  return(metrics)
}

method_palette = c(
  "velocyto" = "#D62728FF",
  "scVelo(sto)" = '#cab2d6',"scVelo(dyn)" = '#6a3d9a',
  "VeloAE" = '#33a02c' , 
  "dynamo" = '#FFB547FF',
  "VeloVAE" = '#b2df8a',
  "MultiVelo"  ="#843C39FF",
  "UniTVelo" = '#bcbd22',
  "DeepVelo(VAE)" = '#1f78b4',
  "cellDancer" = '#e377c2',
  'PhyloVelo'="#E7CB94FF",
  "veloVI" = '#ff9896',
  "LatentVelo" = "#17BECFFF",
  "DeepVelo(GCN)" = '#a6cee3',
  'TFvelo' = "#7B4173FF",
  "STT" = '#ff7f00',
  'Best performance' = 'grey'
)


####-------all-------------------
HuBMAP <- fread('/media/liyaru/LYR/Benchmark/HuBMAP/new_hubmap_global_all.csv') %>% as.data.frame()

HuBMAP <- column_to_rownames(HuBMAP,'file_name')
HuBMAP <- HuBMAP %>% t() %>% as.data.frame()
# HuBMAP$HuBMAP_GCoh_mean <- apply(HuBMAP, 1, mean, na.rm = TRUE)
HuBMAP$Method <- rownames(HuBMAP)
HuBMAP <-  check_method_name(HuBMAP)
# HuBMAP <-  HuBMAP[!is.na(HuBMAP$Method),]
  
max_values <- apply(HuBMAP[,1:176], 2, function(x) max(x, na.rm = TRUE))
HuBMAP['Best performance',] <- c(max_values,'Best performance')

# HuBMAP[HuBMAP$Method == 'best_method','Method'] = 'Best performance'

HuBMAP_long <- pivot_longer(HuBMAP, cols = -Method, 
                            names_to = "datasets", values_to = "GCoh") %>% 
  as.data.frame()

HuBMAP_long$GCoh <- as.numeric(HuBMAP_long$GCoh)

# order
HuBMAP_long <- HuBMAP_long %>%
  group_by(Method) %>%
  mutate(median_value = median(GCoh, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(median_value)

HuBMAP_long$Method <- factor(HuBMAP_long$Method,levels = unique(HuBMAP_long$Method))

# HuBMAP_long[HuBMAP_long$Method == 'Best performance','Clsss'] = 'Best'
# HuBMAP_long[HuBMAP_long$Method != 'Best performance','Clsss'] = 'Method'

p <- ggplot(HuBMAP_long, aes(x = Method, y = GCoh,fill = Method)) +
  geom_boxplot() +
  # coord_flip()+
  scale_fill_manual(values=method_palette)+
  theme_minimal() +
  theme(
    # legend.position  = "bottom",
    axis.title       = element_text(size = 20),
    axis.text        = element_text(size = 14),
    panel.border     = element_rect(fill = NA),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(size = 10, colour = "white"),
    axis.text.x      = element_text(angle = 45, hjust = 1)
    # legend.title     = element_blank(),
    # legend.text      = element_text(size = 14)
  )
  # facet_wrap(~Clsss,nrow = 1,drop = TRUE) 
p


pdf('/media/liyaru/LYR/Benchmark/Fig/7.HuBMAP.pdf',width = 8,height = 5)
p
dev.off()

####-------best -------------------
HuBMAP <- fread('/media/liyaru/LYR/Benchmark/HuBMAP/new_hubmap_global_all.csv') %>% as.data.frame()

HuBMAP <- column_to_rownames(HuBMAP,'file_name')
HuBMAP <- HuBMAP %>% t() %>% as.data.frame()
# HuBMAP$HuBMAP_GCoh_mean <- apply(HuBMAP, 1, mean, na.rm = TRUE)
HuBMAP$Method <- rownames(HuBMAP)
HuBMAP <-  check_method_name(HuBMAP)
# HuBMAP <-  HuBMAP[!is.na(HuBMAP$Method),]

max_values <- apply(HuBMAP[,1:176], 2, function(x) max(x, na.rm = TRUE))
HuBMAP['Best performance',] <- c(max_values,'Best performance')

HuBMAP <- t(HuBMAP) %>% as.data.frame()

b <- data.frame(ID = rownames(HuBMAP)[1:176],
                best_method = HuBMAP$`Best performance`[1:176])

b <- b %>%
  mutate(ID = str_replace(ID, "\\.h5ad$", ""))

a = fread('/media/liyaru/LYR/Diff_change/0_DATA/hubmap-datasets-metadata-2024-08-28_11-51-13.tsv') %>% 
  as.data.frame()

a = a[,c('uuid','origin_samples_unique_mapped_organs')]
colnames(a)[1] = 'ID'

t = a[a$ID %in% b$ID,]

b = merge(a,b,by='ID',all.y=T)
colnames(b)[2] = 'organ'
b$organ <- gsub("'", "", b$organ)  
b$organ <- gsub("\\[|\\]", "", b$organ)  
b$organ <- gsub("\\([^)]*\\)", "", b$organ)

b$best_method <- as.numeric(b$best_method)

b <- b %>%
  group_by(organ) %>%
  mutate(median_value = median(best_method, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(median_value)

b$organ <- factor(b$organ,levels = unique(b$organ))

# pal= pal_npg("nrc")(10)
# scales::show_col(pal)
pal= pal_d3('category20b')(10)
# scales::show_col(pal)

ggplot(b, aes(x = organ, y = best_method,fill=organ))+ 
  geom_boxplot()  +
  # coord_flip()+
  scale_fill_manual(values = pal)+
  #scale_fill_aaas()+
  theme_minimal() +
  theme(
    # legend.position  = "bottom",
    axis.title       = element_text(size = 20),
    axis.text        = element_text(size = 14),
    panel.border     = element_rect(fill = NA),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(size = 10, colour = "white"),
    # legend.title     = element_blank(),
    # legend.text      = element_text(size = 14)
  )+
  RotatedAxis()+
  coord_cartesian(ylim = c(0.5, 1))

pdf('/media/liyaru/LYR/Benchmark/Fig/7.HuBMAP_best_by_tissue.pdf',width = 8,height = 5)
p
dev.off()

####--------heatmap---------------------
library(pheatmap)
library(viridis)
HuBMAP <- fread('/media/liyaru/LYR/Benchmark/HuBMAP/new_hubmap_global_all.csv') %>% as.data.frame()

HuBMAP <- column_to_rownames(HuBMAP,'file_name')
HuBMAP <- HuBMAP %>% t() %>% as.data.frame()
HuBMAP$Method <- rownames(HuBMAP)
HuBMAP <-  check_method_name(HuBMAP)
max_values <- apply(HuBMAP[,1:176], 2, function(x) max(x, na.rm = TRUE))
HuBMAP['Best performance',] <- c(max_values,'Best performance')

rownames(HuBMAP) = HuBMAP$Method
HuBMAP$Method <- factor(HuBMAP$Method,levels = names(method_palette))
HuBMAP <- HuBMAP[order(HuBMAP$Method),1:176]

HuBMAP[1:176] <- lapply(HuBMAP, as.numeric)

annotation_col = data.frame(
  Method = names(method_palette)
)
rownames(annotation_col) = names(method_palette)

# annotation row
a = fread('/media/liyaru/LYR/Diff_change/0_DATA/hubmap-datasets-metadata-2024-08-28_11-51-13.tsv') %>% 
  as.data.frame()
a = a[,c('uuid','origin_samples_unique_mapped_organs')]
colnames(a)[1] = 'ID'

ID <- gsub("\\.h5ad$", "", colnames(HuBMAP))
b = a[a$ID %in% ID,]
colnames(b)[2] = 'organ'
b$organ <- gsub("'", "", b$organ)  
b$organ <- gsub("\\[|\\]", "", b$organ)  
b$organ <- gsub("\\([^)]*\\)", "", b$organ)

annotation_row = data.frame(
  Organ = b$organ
)
rownames(annotation_row) = paste0(b$ID,'.h5ad') 

organ_pal =  c("#2873B3","#2EBEBE","#74B346","#167153","#F1CC2F","#7D4444","#A14462","#8264CC","#D55E00", "#CC79A7")
names(organ_pal) = unique(annotation_row$Organ)

ann_colors = list(
  Method = method_palette,
  Organ = organ_pal 
)

p <- pheatmap(t(HuBMAP[,1:176]),cluster_cols = F,
         show_rownames = F,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         color = viridis(100),
         annotation_row = annotation_row
         # color = brewer.pal(100, "Blues")
         # color = colorRampPalette(c("blue", "red"))(100)
         )

pdf('/media/liyaru/LYR/Benchmark/Fig/7.HuBMAP_heatmap.pdf',width = 6,height = 8)
p
dev.off()









