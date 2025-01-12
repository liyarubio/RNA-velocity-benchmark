library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

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
  'TFvelo' ="#7B4173FF",
  "STT" = '#ff7f00',
  "MultiVelo"  ="#843C39FF",
  'PhyloVelo'="#E7CB94FF"
)

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

####----------- time all hours  -----------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/time/time_new_0103.csv') %>%
  t() %>%
  as.data.frame()

colnames(metrics) <- metrics[1,]
metrics <- metrics[2:16,]
metrics$Method = rownames(metrics)

metrics <- check_method_name(metrics)
metrics[, 1:15] <- sapply(metrics[, 1:15], as.numeric)

metrics_long <- pivot_longer(metrics, cols = -Method, 
                          names_to = "Cell", values_to = "Time") %>% 
  as.data.frame()


metrics_long$Cell <- as.numeric(metrics_long$Cell) / 1000
metrics_long$Time <- as.numeric(metrics_long$Time / 60 / 60) # hours

p <- ggplot(metrics_long ,aes(Cell,Time,group= Method,color=Method))+
  geom_point(size=1)+
  geom_line()+
  scale_color_manual(values = method_palette)+
  # coord_cartesian(ylim = c(0, 200))+
  labs(    
    x = "Cell number",
    y = "Time") +
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
  )+
  scale_x_continuous(breaks = seq(0,60,10))
  # geom_hline(yintercept = 60, linetype = "dashed", color = "black") 
  # scale_y_continuous(breaks = seq(0,1800,60)) 
p

pdf('/media/liyaru/LYR/Benchmark/Fig/11.time_all.pdf')
p
dev.off()

p <- ggplot(metrics_long ,aes(Cell,Time,group= Method,color=Method))+
  geom_point(size=1)+
  geom_line()+
  scale_color_manual(values = method_palette)+
  coord_cartesian(ylim = c(0, 3))+
  labs(    
    x = "Cell number",
    y = "Time") +
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
  ) +
  scale_x_continuous(breaks = seq(0,60,10))
  #scale_y_continuous(breaks = c(0, 60, 120, 180)) 
p

pdf('/media/liyaru/LYR/Benchmark/Fig/11.time_0_200.pdf')
p
dev.off()


####----------- heatmap  -----------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/time/time_new_0103.csv') %>%
  t() %>%
  as.data.frame()

colnames(metrics) <- metrics[1,]
metrics <- metrics[2:16,]
metrics$Method = rownames(metrics)

metrics <- check_method_name(metrics)
metrics[, 1:15] <- sapply(metrics[, 1:15], as.numeric)
metrics_long <- pivot_longer(metrics, cols = -Method, 
                             names_to = "Cell", values_to = "Time") %>% 
  as.data.frame()

metrics_long$Cell <-as.numeric(metrics_long$Cell) / 1000 
metrics_long$Cell <- paste0(metrics_long$Cell,'K')
metrics_long$Cell <- factor(metrics_long$Cell,levels=paste0(rev(c(1:10,20,30,40,50,60)),'K'))


metrics_long$Time <- as.numeric(metrics_long$Time / 60 ) # Minutes

# metrics_long$Label = ceiling(metrics_long$Time)

metrics_long[!is.na(metrics_long$Time) & (metrics_long$Time < 60),'Label'] = 
  paste0(ceiling(metrics_long[!is.na(metrics_long$Time) & (metrics_long$Time < 60),'Time']),'m') 

metrics_long[!is.na(metrics_long$Time) & (metrics_long$Time > 60),'Label'] = 
  paste0(ceiling(metrics_long[!is.na(metrics_long$Time) & (metrics_long$Time > 60),'Time']/60),'h') 

library(RColorBrewer)
brewer.pal(9, 'Oranges')

p <- ggplot()+
  geom_tile(data=metrics_long,aes(Method,Cell),fill="white",color="grey")+
  geom_point(data=metrics_long,aes(Method,Cell,color=Time,size=10,alpha=0.7,stroke = 8),shape=15)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0,size = 12,color = "black"),
        axis.text.y = element_text(size = 10,color = "black"))+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "right")+
  labs(x=NULL,y=NULL)+
  scale_colour_gradient(low ="#FEE6CE", high = "#D94801",
                        limits = c(0, 60), oob = scales::squish, 
                        na.value = "grey") +
  # scale_colour_viridis_c(limits = c(0, 60),
  #                        oob = scales::squish, na.value = "grey50", guide = guide_legend(override.aes = list(color = NA)))+
  
  geom_text(data = metrics_long, aes(Method, Cell, label = Label), vjust = 0.5, size = 3)


pdf('/media/liyaru/LYR/Benchmark/Fig/11.time_new.pdf')
p
dev.off()



