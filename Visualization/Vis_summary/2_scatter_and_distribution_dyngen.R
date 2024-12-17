library(data.table)
library(tidyr)
library(reshape2) 
library(patchwork)
library(ggpubr)

`%>%` <- magrittr::`%>%`


####-------function-----------------
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
  metrics <- metrics[order(metrics$Method),]
  
  return(metrics)
}



####----- color----------------------
setwd('/media/liyaru/LYR/Diff_change/2_benchmark/25_multiomic_benchmark/scib-reproducibility-main/visualization')

# pal_d3("category20")(20)
# pal_d3("category20b")(20)

method_palette = c(
  "velocyto" = "#D62728FF",
  "scVelo(sto)" = '#cab2d6',"scVelo(dyn)" = '#6a3d9a',
  "VeloAE" = '#33a02c' , 
  "dynamo" = '#FFB547FF',
  "VeloVAE" = '#b2df8a',
  "UniTVelo" = '#bcbd22',
  "DeepVelo(VAE)" = '#1f78b4',
  "cellDancer" = '#e377c2',
  "veloVI" = '#ff9896',
  "LatentVelo" = "#17BECFFF",
  "DeepVelo(GCN)" = '#a6cee3',
  "STT" = '#ff7f00',
  "MultiVelo"  ="#843C39FF",
  'PhyloVelo'="#E7CB94FF",
  'TFvelo' ="#7B4173FF"
)



####----- Cosine high vs low ------------
df1 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv')
df1_long <- melt(df1, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Cosine_V')
colnames(df1_long)[1] = 'Method'
df1_long$Method_Dataset <- paste0(df1_long$Method,".",df1_long$Dataset)


df2 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_pca2.csv')
df2_long <- melt(df2, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Global')
colnames(df2_long)[1] = 'Method'
df2_long$Method_Dataset <- paste0(df2_long$Method,".",df2_long$Dataset)

df_merge <- merge(df1_long[,3:4],df2_long,by='Method_Dataset',all.x = T,all.y=T)

df_merge <- df_merge %>%
  mutate(Method= gsub("^[0-9]+_", "", Method)) %>%
  filter(Method != "11_Dynamo_deterministic")

df_merge <- check_method_name(df_merge)

df_merge$Method <- factor(df_merge$Method,levels = names(method_palette))

table(df_merge$Method)

df_summ <- df_merge %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    CosineMean = mean(`Cosine_V`,na.rm=T),
    CosineSD   = sd(`Cosine_V`,na.rm=T),
    GlobalMean   = mean(`Global`,na.rm=T),
    GlobalSD     = sd(`Global`,na.rm=T),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(Method)

df_summ <- df_summ[df_summ$Method != 'VeloAE',]


summary_plot <- ggplot(df_summ) +
  aes(
    x      = CosineMean,
    y      = GlobalMean,
    colour = Method
  ) +
  geom_errorbar(
    aes(ymin = GlobalMean   - GlobalSD, ymax = GlobalMean   + GlobalSD)
  ) +
  geom_errorbarh(
    aes(xmin = CosineMean - CosineSD, xmax = CosineMean + CosineSD)
  ) +
  geom_point(size = 3, stroke = 1, fill = "white") +
  scale_color_manual(values = method_palette)+
  #coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-1, 1)) +
  labs(
    x = "High-dim Cosine Similarity",
    y = "Low-dim Cosine Similarity"
  ) +
  guides(
    colour = guide_legend(
      title          = "Method",
      title.position = "top",
      ncol           = 2,
      order          = 10
    )
  ) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    axis.title       = element_text(size = 20),
    axis.text        = element_text(size = 14),
    panel.border     = element_rect(fill = NA),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(size = 10, colour = "white"),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 14)
  )

summary_plot

# save fig
#now_str <- format(lubridate::now(), "%Y%m%d_%H%M%S")
for (extension in c("pdf")) { # , "tiff", "png"
  ggsave(
    as.character(glue::glue(
      "/media/liyaru/LYR/Benchmark/Fig/2.dyngen_scatter_Cosine_high_and_low.{extension}")),
    summary_plot,
    width = 210, height = 297, units = "mm"
  )
}

####----- Cosine high vs low split by dataset ------------
df1 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv')
df1_long <- melt(df1, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Cosine_V')
colnames(df1_long)[1] = 'Method'
df1_long$Method_Dataset <- paste0(df1_long$Method,".",df1_long$Dataset)


df2 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_pca2.csv')
df2_long <- melt(df2, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Global')
colnames(df2_long)[1] = 'Method'
df2_long$Method_Dataset <- paste0(df2_long$Method,".",df2_long$Dataset)

df_merge <- merge(df1_long[,3:4],df2_long,by='Method_Dataset',all.x = T,all.y=T)

df_merge <- df_merge %>%
  mutate(Method= gsub("^[0-9]+_", "", Method)) %>%
  filter(Method != "11_Dynamo_deterministic")

df_merge <- check_method_name(df_merge)

df_merge$Method <- factor(df_merge$Method,levels = names(method_palette))

df_merge$Dataset_type <- sub("_seed\\d+$", "", df_merge$Dataset)
Dataset_type <- unique(df_merge$Dataset_type)

df_all = df_merge

plots <- list() 

for (d in Dataset_type){
  
  print(d)
  df_merge = df_all[df_all$Dataset_type ==d,]
  
  df_summ <- df_merge %>%
    dplyr::group_by(Method) %>%
    dplyr::summarise(
      CosineMean = mean(`Cosine_V`,na.rm=T),
      CosineSD   = sd(`Cosine_V`,na.rm=T),
      GlobalMean   = mean(`Global`,na.rm=T),
      GlobalSD     = sd(`Global`,na.rm=T),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(Method)
  
  df_summ <- df_summ[df_summ$Method != 'VeloAE',]
  
  summary_plot <- ggplot(df_summ) +
    aes(
      x      = CosineMean,
      y      = GlobalMean,
      colour = Method
    ) +
    geom_errorbar(
      aes(ymin = GlobalMean   - GlobalSD, ymax = GlobalMean   + GlobalSD)
    ) +
    geom_errorbarh(
      aes(xmin = CosineMean - CosineSD, xmax = CosineMean + CosineSD)
    ) +
    geom_point(size = 2, stroke = 1, fill = "white") +
    scale_color_manual(values = method_palette)+
    #coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-1, 1)) +
    labs(
      x = "",
      y = "",
      title = d
    ) +
    guides(
      colour = guide_legend(
        title          = "Method",
        title.position = "top",
        ncol           = 2,
        order          = 10
      )
    ) +
    theme_minimal() +
    theme(
      # legend.position  = "bottom",
      legend.position  = "none",
      axis.title       = element_text(size = 8),
      axis.text        = element_text(size = 8),
      panel.border     = element_rect(fill = NA),
      strip.background = element_rect(fill = "black"),
      strip.text       = element_text(size = 8, colour = "white"),
      legend.title     = element_blank(),
      legend.text      = element_text(size = 8),
      plot.title       = element_text(hjust = 0.5, size = 12) 
    )
  
  plots[[d]] <- summary_plot 

}

combined_plot <- wrap_plots(plots, ncol = 5, nrow = 3)
combined_plot 

pdf("/media/liyaru/LYR/Benchmark/Fig/2.dyngen_scatter_Cosine_high_and_low_split.pdf",width = 15,height = 9)
combined_plot 
dev.off()


####--------distribution boxplot----------------
df1 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv')
df1_long <- melt(df1, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'High')
colnames(df1_long)[1] = 'Method'
df1_long$Method_Dataset <- paste0(df1_long$Method,".",df1_long$Dataset)


df2 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_pca2.csv')
df2_long <- melt(df2, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Low')
colnames(df2_long)[1] = 'Method'
df2_long$Method_Dataset <- paste0(df2_long$Method,".",df2_long$Dataset)

df_merge <- merge(df1_long[,3:4],df2_long,by='Method_Dataset',all.x = T,all.y=T)

df_merge <- df_merge %>%
  mutate(Method= gsub("^[0-9]+_", "", Method)) %>%
  filter(Method != "11_Dynamo_deterministic")

df_merge <- df_merge[,c('Method_Dataset','High','Low')]

df_merge_long <- melt(df_merge, id.vars = c("Method_Dataset"), 
                      variable.name = "Dim", value.name = "Cosine_sim")

# ggpaired(df_merge_long, 
#          x = "Dim", 
#          y = "Cosine_sim",
#          color = "Dim", 
#          line.color = "gray", 
#          line.size = 0.4,
#          palette = "jco") 


p <- ggplot(df_merge_long, aes(x = Dim, y = Cosine_sim, fill = Dim)) +
  geom_violin(size = 1) +
  geom_line(aes(group = Method_Dataset), color = "grey80", size = 0.1) +
  geom_point(size = 0.1) +
  ggsci::scale_color_npg() +
  scale_fill_manual(values = c("#08306B","#C6DBEF"))+
  theme_classic2(base_size = 20)+
  theme(
    legend.position = "none"
  )

pdf("/media/liyaru/LYR/Benchmark/Fig/2.dyngen_distribution.pdf",width = 4,height = 6)
p
dev.off()

####------ distribution boxplot split by method ---------------
df_merge_long <- separate(df_merge_long,
                                 col = 'Method_Dataset',sep='[.]',
                                 into=c('Method','Dataset'),remove = F)

df_merge_long <- df_merge_long %>%
  mutate(Method= gsub("^[0-9]+_", "", Method))

df_merge_long <- check_method_name(df_merge_long)

df_merge_long$Method <- factor(df_merge_long$Method,levels = names(method_palette))

p <- ggplot(df_merge_long, aes(x = Dim, y = Cosine_sim, fill = Dim)) +
  geom_violin(size = 1) +
  geom_line(aes(group = Method_Dataset), color = "grey80", size = 0.5) +
  geom_point(size = 0.5) +
  ggsci::scale_color_npg() +
  scale_fill_manual(values = c("#08306B","#C6DBEF"))+
  theme_classic2(base_size = 15)+
  facet_wrap(~Method,nrow = 2) 

pdf("/media/liyaru/LYR/Benchmark/Fig/2.dyngen_distribution_split.pdf",width = 16,height = 6)
p
dev.off()


####--------distribution boxplot add global ----------------
df1 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv')
df1_long <- melt(df1, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'High')
colnames(df1_long)[1] = 'Method'
df1_long$Method_Dataset <- paste0(df1_long$Method,".",df1_long$Dataset)


df2 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_pca2.csv')
df2_long <- melt(df2, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Low')
colnames(df2_long)[1] = 'Method'
df2_long$Method_Dataset <- paste0(df2_long$Method,".",df2_long$Dataset)


df3 <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/eval_global.csv')
df3_long <- melt(df3, id.vars = c("V1"),
                 variable.name = c('Dataset'),value.name = 'Global')
colnames(df3_long)[1] = 'Method'
df3_long$Method_Dataset <- paste0(df3_long$Method,".",df3_long$Dataset)

df_merge <- merge(df1_long[,3:4],df2_long,by='Method_Dataset',all.x = T,all.y=T)
df_merge <- merge(df_merge,df3_long[,3:4],by='Method_Dataset',all.x = T,all.y=T)

df_merge <- df_merge %>%
  mutate(Method= gsub("^[0-9]+_", "", Method)) %>%
  filter(Method != "11_Dynamo_deterministic")

df_merge <- df_merge[,c('Method_Dataset','High','Low','Global')]

df_merge_long <- melt(df_merge, id.vars = c("Method_Dataset"),
                      variable.name = "Dim", value.name = "Cosine_sim")

p <- ggplot(df_merge_long, aes(x = Dim, y = Cosine_sim, fill = Dim)) +
  geom_violin(size = 1) +
  geom_line(aes(group = Method_Dataset), color = "grey80", size = 0.1) +
  geom_point(size = 0.1) +
  ggsci::scale_color_npg() +
  # scale_fill_manual(values = c("#08306B","#C6DBEF"))+
  theme_classic2(base_size = 20)+
  theme(
    legend.position = "none"
  )

df_merge_long <- separate(df_merge_long,
                          col = 'Method_Dataset',sep='[.]',
                          into=c('Method','Dataset'),remove = F)

df_merge_long <- df_merge_long %>%
  mutate(Method= gsub("^[0-9]+_", "", Method))

df_merge_long <- check_method_name(df_merge_long)

df_merge_long$Method <- factor(df_merge_long$Method,levels = names(method_palette))

p <- ggplot(df_merge_long, aes(x = Dim, y = Cosine_sim, fill = Dim)) +
  geom_violin(size = 1) +
  geom_line(aes(group = Method_Dataset), color = "grey80", size = 0.5) +
  geom_point(size = 0.5) +
  ggsci::scale_color_npg() +
  # scale_fill_manual(values = c("#08306B","#C6DBEF"))+
  theme_classic2(base_size = 15)+
  facet_wrap(~Method,nrow = 2)


p <- ggplot(df_merge_long, aes(x = Dim, y = Cosine_sim, fill = Dim)) +
  geom_violin(size = 1) +
  geom_line(aes(group = Method_Dataset), color = "grey80", size = 0.5) +
  geom_point(size = 0.5) +
  ggsci::scale_color_npg() +
  # scale_fill_manual(values = c("#08306B","#C6DBEF"))+
  theme_classic2(base_size = 15)+
  facet_wrap(~Dataset,nrow = 4)

# pdf("/media/liyaru/LYR/Benchmark/Fig/2.dyngen_distribution_split.pdf",width = 16,height = 6)
# p
# dev.off()
