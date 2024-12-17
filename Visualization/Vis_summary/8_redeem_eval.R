library(ggplot2)
library(data.table)
library(tidyr)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

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
  "STT" = '#ff7f00'
)

palettes = list('Metrics' = 'Greens',
                'Score' = 'Greens')

####-----------function-------------------
check_method_name <- function(metrics){
  metrics <- metrics %>%
    filter(Method != "Dynamo_deterministic") %>%
    mutate(Method = gsub("Velocyto", "velocyto", Method,ignore.case = T))  %>%
    mutate(Method = gsub("scvelo_dynamical" , "scVelo(dyn)", Method,ignore.case = T))  %>%
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
  
  metrics$Method <- factor(metrics$Method,levels = names(method_palette))
  metrics <- metrics[!is.na(metrics$Method),]
  metrics <- metrics[order(metrics$Method),]
  
  return(metrics)
}

# add defaults for optional values
add_column_if_missing <- function(df, ...) {
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))
    
    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), default_val, df[[column_name]])
    } else {
      df[[column_name]] <- default_val
    }
  }
  df
}

####-----------heatmap-------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/Redeem/eval_csv/all_eval_CBDir.csv') %>%
  t() %>%
  as.data.frame()

colnames(metrics) <- metrics[1,]
metrics <- metrics[2:nrow(metrics),]
metrics$Method = rownames(metrics)
metrics <- check_method_name(metrics)
rownames(metrics) <- metrics$Method

metrics<- metrics[,1:11] %>% mutate_all(~ as.numeric(.))

breaks <- seq(-1, 1, length.out = 100)

pheatmap(metrics,breaks = breaks)

pheatmap(metrics,breaks = breaks,cluster_rows = F,cluster_cols = F)


####-----------bar plot -------------------

####------read data-----------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/Redeem/eval_csv/all_eval_CBDir.csv') %>%
  t() %>%
  as.data.frame()

colnames(metrics) <- metrics[1,]
metrics <- metrics[2:nrow(metrics),]
metrics$Method = rownames(metrics)
metrics <- check_method_name(metrics)
rownames(metrics) <- metrics$Method

cols = 1:11
metrics[,cols] <- metrics[,cols] %>% mutate_all(~ as.numeric(.))

metrics$Method <- factor(metrics$Method,levels = names(method_palette))

metrics <- metrics[order(metrics$Method),]

####-----------bar plot-----------------
library(ggplot2)
library(reshape2)
library(gridExtra)

df <- metrics %>% t() %>% as.data.frame()
df <- df[1:11,]
df <- df %>% mutate_all(~ as.numeric(.))
df$celltype <- rownames(df)

df$celltype <- gsub("\\(|\\)", "", df$celltype)
df$celltype <- gsub("'", "", df$celltype)
df$celltype <- gsub(",", " →", df$celltype)

color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

create_barplot <- function(col) {
  t = df[,c(col,'celltype')] %>% as.data.frame()
  t['value'] = t[col]
  
  t$celltype <- factor(t$celltype,levels = rev(t$celltype))
  
  p <- ggplot(t,aes(x=celltype,y=value,fill=value)) +
    geom_bar(stat = "identity",width = 0.8,color='black') +
    coord_flip()+
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 1))+
    labs(title =col) + 
    scale_fill_gradientn(colors = color_palette,limits = c(-1, 1))
  
  # if (col == names(df)[1]) {
  #   p + theme(axis.text.y = element_text())
  # } else {
  #   p
  # }
  
}

plots <- lapply(names(df)[1:16], create_barplot)

pdf('/media/liyaru/LYR/Benchmark/Fig/12.redeem_CBDir.pdf',width = 12,height = 5)
do.call(grid.arrange, c(plots, ncol = 8)) 
dev.off()

# legend and celltype names
col = 'velocyto'
t = df[,c(col,'celltype')] %>% as.data.frame()
t['value'] = t[col]

t$celltype <- factor(t$celltype,levels = rev(t$celltype))

p <- ggplot(t,aes(x=celltype,y=value,fill=value)) +
  geom_bar(stat = "identity",width = 0.8,color='black') +
  coord_flip()+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 1))+
  labs(title =col)+ 
  scale_fill_gradientn(colors = color_palette,limits = c(-1, 1))
p

pdf('/media/liyaru/LYR/Benchmark/Fig/12.redeem_CBDir_legend.pdf',width = 5,height = 5)
p
dev.off()