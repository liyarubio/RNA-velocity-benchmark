library(tidyr)
library(data.table)

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


####-----------cosine similarity-----------------------

metrics <- fread('/media/liyaru/LYR/Benchmark/simu_celldancer/celldancer_cosine_sim_all_reorder.csv') %>%
  as.data.frame()

metrics$Method <- factor(metrics$Method,levels = names(method_palette))

# metrics <- metrics[,c(1,22:31)] # multi forward
# metrics <- metrics[,c(1,2:11)] # mono
# metrics <- metrics[,c(1,12:21)] # boost
# metrics <- metrics[,c(1,32:41)] # multi backward

metrics <- metrics[,c(1,2:31)] 
metrics <- metrics[metrics$Method != 'LatentVelo',]


metrics_long <- pivot_longer(metrics, cols = -Method, 
                          names_to = "simulation", values_to = "Cosine") %>% 
  as.data.frame()
metrics_long$noise_level = sub(".*_(\\d+)$", "\\1", metrics_long$simulation)
metrics_long$noise_value = as.integer(metrics_long$noise_level) * 0.2
# metrics_long$group <- sub("_.*$", "", metrics_long$simulation)

df_summ <- metrics_long %>%
  dplyr::group_by(Method,noise_value) %>%
  dplyr::summarise(
    CosineMean = mean(`Cosine`,na.rm=T),
    CosineMedian = median(`Cosine`,na.rm=T),
    CosineSD   = sd(`Cosine`,na.rm=T),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(noise_value)

df_summ_wide <- pivot_wider(df_summ[,c("Method","noise_value","CosineMean")], 
                            names_from = noise_value, values_from = CosineMean)

fwrite(df_summ_wide,'/media/liyaru/LYR/Benchmark/simu_celldancer/simu_noise_mean.csv')

p <- ggplot(df_summ ,aes(noise_value,CosineMean,group= Method,color=Method))+
  geom_point(size=2)+
  geom_line()+
  scale_color_manual(values = method_palette)+
  #coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-1, 1)) +
  labs(    
    x = "Noise Level",
    y = "Cosine Similarity") +
  # geom_errorbar(aes(ymin = CosineMean - CosineSD, ymax = CosineMean + CosineSD),
  #               width = 0.15,cex=1.5)+
  guides(
    colour = guide_legend(
      title          = "Method",
      title.position = "top",
      ncol           = 2,
      order          = 10
    )
  ) +
  theme_minimal() +
  # theme_bw()+
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

pdf('/media/liyaru/LYR/Benchmark/Fig/6_noise_cosine.pdf')
p
dev.off()

####----------- Raito 1 : cal ratio first -----------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/simu_celldancer/celldancer_cosine_sim_all_reorder.csv') %>%
  as.data.frame()

metrics$Method <- factor(metrics$Method,levels = names(method_palette))

# metrics <- metrics[,c(1,22:31)] # multi forward
# metrics <- metrics[,c(1,2:11)] # mono
# metrics <- metrics[,c(1,12:21)] # boost
# metrics <- metrics[,c(1,32:41)] # multi backward

metrics <- metrics[,c(1,2:31)] 
metrics <- metrics[metrics$Method != 'LatentVelo',]

# ratio
metrics[,2:11] <- metrics[,2:11]/metrics[,2]
metrics[,12:21] <- metrics[,12:21]/metrics[,12]
metrics[,22:31] <- metrics[,22:31]/metrics[,22]

fwrite(metrics,'/media/liyaru/LYR/Benchmark/simu_celldancer/celldancer_cosine_sim_percentage.csv')

metrics_long <- pivot_longer(metrics, cols = -Method, 
                             names_to = "simulation", values_to = "Cosine") %>% 
  as.data.frame()
metrics_long$noise_level = sub(".*_(\\d+)$", "\\1", metrics_long$simulation)
metrics_long$noise_value = as.integer(metrics_long$noise_level) * 0.2
# metrics_long$group <- sub("_.*$", "", metrics_long$simulation)

df_summ <- metrics_long %>%
  dplyr::group_by(Method,noise_value) %>%
  dplyr::summarise(
    CosineMean = mean(`Cosine`,na.rm=T),
    CosineMedian = median(`Cosine`,na.rm=T),
    CosineSD   = sd(`Cosine`,na.rm=T),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(noise_value)

df_summ_wide <- pivot_wider(df_summ[,c("Method","noise_value","CosineMean")], 
                            names_from = noise_value, values_from = CosineMean)
fwrite(df_summ_wide,'/media/liyaru/LYR/Benchmark/simu_celldancer/simu_noise_mean_percentage.csv')

p <- ggplot(df_summ ,aes(noise_value,CosineMean*100,group= Method,color=Method))+
  geom_point(size=2)+
  geom_line()+
  scale_color_manual(values = method_palette)+
  #coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-1, 1)) +
  labs(    
    x = "Noise Level",
    y = "Percentage") +
  # geom_errorbar(aes(ymin = CosineMean - CosineSD, ymax = CosineMean + CosineSD),
  #               width = 0.15,cex=1.5)+
  guides(
    colour = guide_legend(
      title          = "Method",
      title.position = "top",
      ncol           = 2,
      order          = 10
    )
  ) +
  theme_minimal() +
  # theme_bw()+
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


pdf('/media/liyaru/LYR/Benchmark/Fig/6_noise_ratio.pdf')
p
dev.off()


# ####--------dotplot----------------
# metrics <- fread('/media/liyaru/LYR/Benchmark/simu_celldancer/celldancer_cosine_sim_all_reorder.csv') %>%
#   as.data.frame()
# 
# metrics$Method <- factor(metrics$Method,levels = names(method_palette))
# 
# # metrics <- metrics[,c(1,22:31)] # multi forward
# # metrics <- metrics[,c(1,2:11)] # mono
# # metrics <- metrics[,c(1,12:21)] # boost
# # metrics <- metrics[,c(1,32:41)] # multi backward
# 
# metrics <- metrics[,c(1,2:31)] 
# metrics <- metrics[metrics$Method != 'LatentVelo',]
# 
# metrics_long <- pivot_longer(metrics, cols = -Method, 
#                              names_to = "simulation", values_to = "Cosine") %>% 
#   as.data.frame()
# metrics_long$noise_level = sub(".*_(\\d+)$", "\\1", metrics_long$simulation)
# metrics_long$noise_value = as.character(as.integer(metrics_long$noise_level) * 0.2)  
# 
# brewer.pal(9, 'RdBu')
# brewer.pal()
# 
# p <- ggplot()+
#   geom_tile(data=metrics_long,aes(simulation,Method),fill="white",color="grey")+
#   geom_point(data=metrics_long,aes(simulation,Method,color=Cosine,size=Cosine,stroke =1),shape=15)+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0,size = 12,color = "black"),
#         axis.text.y = element_text(size = 10,color = "black"))+
#   scale_x_discrete(position = "top")+
#   scale_y_discrete(position = "right")+
#   labs(x=NULL,y=NULL)+
#   scale_colour_gradient2(low ="#2166AC", mid = "#F7F7F7",high ="#B2182B", na.value = "grey",limits = c(-1,1))
#   #scale_colour_viridis_c(limits = c(-1,1),oob = scales::squish, na.value = "grey50")
# p

###--------heatmap----------
library(pheatmap)
metrics <- fread('/media/liyaru/LYR/Benchmark/simu_celldancer/celldancer_cosine_sim_all_reorder.csv') %>%
  as.data.frame()

metrics$Method <- factor(metrics$Method,levels = names(method_palette))

metrics2 <- metrics[,2:31]
rownames(metrics2) <- metrics$Method

annotation_col = data.frame(
  data = colnames(metrics2)
)
annotation_col$noise_level = sub(".*_(\\d+)$", "\\1",annotation_col$data)
annotation_col$noise_level = as.integer(annotation_col$noise_level) * 0.2
annotation_col$data_type <- gsub("_[0-9]+$", "", annotation_col$data)
rownames(annotation_col) = annotation_col$data
annotation_col <- annotation_col[,2:3]

display_matrix <- metrics2 %>% mutate_all(~round(., 1))
display_matrix <- replace(display_matrix, is.na(metrics2), "")

data_type =  c("#DADAEB", "#9E9AC8","#54278F")
names(data_type) = unique(annotation_col$data_type)
grad_colors <- colorRampPalette(c("#F7FCF5" , "#00441B"))(100)

ann_colors = list(
  data_type = data_type,
  noise_level = grad_colors 
)

pdf('/media/liyaru/LYR/Benchmark/Fig/5.simu_cd_cosine_heatmap.pdf',width = 12,height = 4)
pheatmap(metrics2,cluster_cols = F,cluster_rows = F,
         annotation_col = annotation_col,
         display_numbers = display_matrix,
         show_colnames = F,
         annotation_colors = ann_colors)
dev.off()



####----------- Raito 2 : cal avg first -----------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/simu_celldancer/simu_noise_mean.csv',header = T) %>%
  as.data.frame()

metrics$Method <- factor(metrics$Method,levels = names(method_palette))

# ratio
metrics[,2:11] <- metrics[,2:11]/metrics[,2]

metrics_long <- pivot_longer(metrics, cols = -Method, 
                             names_to = "simulation", values_to = "Cosine") %>% 
  as.data.frame()

metrics_long$noise_value = as.numeric(metrics_long$simulation) 

p <- ggplot(metrics_long ,aes(noise_value,Cosine*100,group= Method,color=Method))+
  geom_point(size=2)+
  geom_line()+
  scale_color_manual(values = method_palette)+
  #coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-1, 1)) +
  labs(    
    x = "Noise Level",
    y = "Percentage") +
  # geom_errorbar(aes(ymin = CosineMean - CosineSD, ymax = CosineMean + CosineSD),
  #               width = 0.15,cex=1.5)+
  guides(
    colour = guide_legend(
      title          = "Method",
      title.position = "top",
      ncol           = 2,
      order          = 10
    )
  ) +
  theme_minimal() +
  # theme_bw()+
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


pdf('/media/liyaru/LYR/Benchmark/Fig/6_noise_ratio2.pdf')
p
dev.off()









