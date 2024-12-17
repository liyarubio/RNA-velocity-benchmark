library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)
library(data.table)
library(dplyr)
library(scales)
library(ggimage)
library(cowplot)
library(tidyr)

#source("knit_table.R") # Please put knit_table.R in your working dir
setwd('/media/liyaru/LYR/Diff_change/2_benchmark/25_multiomic_benchmark/scib-reproducibility-main/visualization')

####------color-----------
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
  # "MultiVelo"  ="#843C39FF",
  # 'PhyloVelo'="#E7CB94FF",
  'TFvelo' ="#7B4173FF",
  "STT" = '#ff7f00'
)

palettes = list('Accuracy' = 'Oranges',
                'Accuracy2' = 'Purples',
                'Scalability' = 'Greens',
                'Stability' = 'Blues',
                'Overall' = 'Reds')

# palettes = list('Accuracy' = 'Reds',
#                 'Simu1' = 'Blues',
#                 'Simu2' = 'BuPu',
#                 'Real' = 'RdPu',
#                 'Multi' = 'OrRd',
#                 'Time' = 'Greens',
#                 'Stability' = 'YlOrBr')

####-----function---------------------
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

normalize <- function(x){
  # ignore na
  x_min <- min(x, na.rm = TRUE)-0.01
  x_max <- max(x, na.rm = TRUE)
  
  # ignore 0
  if (x_max == x_min) {
    return(rep(0, length(x))) 
  } else {
    return ((x - x_min) / (x_max - x_min))
  }
}

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
  
  metrics$Method <- factor(metrics$Method,levels = names(method_palette))
  metrics <- metrics[order(metrics$Method),]
  metrics <- metrics[!is.na(metrics$Method),]
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
  "STT" = '#ff7f00'
)

process_metrics <- function(metrics){
  colnames(metrics)[1] <- 'Method'
  #metrics <- separate(metrics,col='Method',into=c("no.",'Method'))
  
  metrics <- metrics %>%
    mutate(Method= gsub("^[0-9]+_", "", Method))
  
  metrics <- metrics %>% check_method_name()
  
  metrics$Method <- factor(metrics$Method,levels = names(method_palette))
  metrics <- metrics[order(metrics$Method)]
  
  metrics$Median <- apply(metrics[,2:43], 1, median, na.rm = TRUE)
  metrics$Mean <- apply(metrics[,2:43], 1, mean, na.rm = TRUE)
  
  return(metrics)
}

####------build data for plot-----------------------
####------dyngen simulation---------------------------

cos_v <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv') %>% 
  process_metrics()

####------simu noise / stability ----------------------------
simu_cd_all <- fread('/media/liyaru/LYR/Benchmark/simu_celldancer/celldancer_cosine_sim_all_reorder.csv') %>%
  as.data.frame()
simu_cd_all <- simu_cd_all[simu_cd_all$Method != 'LatentVelo',]

# simu_cd <- simu_cd_all[,c("Method","Mono-kinetic_1","Transcriptional_boost_1","Multi-forward_branching_1")]
# simu_cd$Stability <- apply(simu_noise[,2:4], 1, mean, na.rm = TRUE)


simu_cd_all[,2:11] <- simu_cd_all[,2:11]/simu_cd_all[,2]
simu_cd_all[,12:21] <- simu_cd_all[,12:21]/simu_cd_all[,12]
simu_cd_all[,22:31] <- simu_cd_all[,22:31]/simu_cd_all[,22]

simu_noise <- simu_cd_all[,c('Method',
                             "Mono-kinetic_2","Transcriptional_boost_2","Multi-forward_branching_2",
                             "Mono-kinetic_5", "Transcriptional_boost_5","Multi-forward_branching_5",
                             "Mono-kinetic_10","Transcriptional_boost_10","Multi-forward_branching_10")]

simu_noise$noise_0.4 <- apply(simu_noise[,c( "Mono-kinetic_2","Transcriptional_boost_2","Multi-forward_branching_2")], 1, mean, na.rm = TRUE)
simu_noise$noise_1 <- apply(simu_noise[,c(  "Mono-kinetic_5", "Transcriptional_boost_5","Multi-forward_branching_5")], 1, mean, na.rm = TRUE)
simu_noise$noise_2 <- apply(simu_noise[,c( "Mono-kinetic_10","Transcriptional_boost_10","Multi-forward_branching_10")], 1, mean, na.rm = TRUE)

simu_noise$Stability <- apply(simu_noise[,c('noise_0.4','noise_1','noise_2')], 1, mean, na.rm = TRUE)

# x_min <- min(simu_noise[,c('noise_0.4','noise_1','noise_2')], na.rm = TRUE)-0.1
# x_max <- max(simu_noise[,c('noise_0.4','noise_1','noise_2')], na.rm = TRUE)
# simu_noise[,c('noise_0.4','noise_1','noise_2')] <- (simu_noise[,c('noise_0.4','noise_1','noise_2')] - x_min) / (x_max - x_min)

# simu_noise$Stability  <- normalize(simu_noise$Stability )

####-----HuBAMP data-----------------------------
HuBMAP <- fread('/media/liyaru/LYR/Benchmark/HuBMAP/new_hubmap_global_all.csv') %>% as.data.frame()

HuBMAP <- column_to_rownames(HuBMAP,'file_name')
HuBMAP <- HuBMAP %>% t() %>% as.data.frame()
HuBMAP$Method <- rownames(HuBMAP)
HuBMAP <-  check_method_name(HuBMAP)
rownames(HuBMAP) = HuBMAP$Method
HuBMAP$Method <- factor(HuBMAP$Method,levels = names(method_palette))
HuBMAP <- HuBMAP[order(HuBMAP$Method),]
HuBMAP$HuBMAP_GCoh_mean <- apply(HuBMAP[,1:176], 1, function(x) mean(x, na.rm = TRUE))

####-------Datasets10-------------------
data10  <- fread('/media/liyaru/LYR/Benchmark/datasets10/global_1205.csv') %>% as.data.frame()

colnames(data10) <- c('Method','Pancreatic Endocrinogenesis',
                       'DentateGyrus (P0 and P5)',
                       'CD34+ Bone Marrow (BM)',
                       'Forebrain','DentateGyrus (P12 and P35)',
                       'Hindbrain','Organogenesis',
                       'Gastrulation Erythroid','Zebrafish',
                       'Gastrulation')

datasets <- c('DentateGyrus (P0 and P5)','DentateGyrus (P12 and P35)',
              'Pancreatic Endocrinogenesis','CD34+ Bone Marrow (BM)',
              'Gastrulation','Gastrulation Erythroid',
              'Forebrain','Hindbrain',
              'Organogenesis','Zebrafish')

data10 <- data10[,c('Method',datasets)]

data10[] <- lapply(data10, function(x) gsub('"', '', x))
data10[,2:11] <- lapply(data10[,2:11], as.numeric)

data10 <- check_method_name(data10)

data10$Method <- factor(data10$Method,levels = names(method_palette))

data10<- data10[order(data10$Method),]

data10 <- as.data.frame(data10)
rownames(data10) <- data10$Method

data10$data10_GCoh_mean <- apply(data10[,2:11], 1, mean,na.rm = TRUE)

####-------Redeem (lineage + scATAC)---------------
redeem <- fread('/media/liyaru/LYR/Benchmark/Redeem/eval_csv/all_eval_8.csv') %>% as.data.frame()
redeem <- redeem[,c(4,1)]
colnames(redeem) <- c("Method",'Redeem_GCoh')
redeem <-  check_method_name(redeem)


#####------Multi-tag--------------------------------
tag <- fread('/media/liyaru/LYR/Benchmark/multitag/Clone_info_low_dim_sim.csv')%>%
  t() %>% as.data.frame()

tag$multitag_accurcy <- rowMeans(tag)
tag$Method <- rownames(tag)

tag <- tag[2:15,c('Method','multitag_accurcy')]

tag <-  check_method_name(tag)


####------ metobalic labling data (scNT-seq)--------
scNT <- fread('/media/liyaru/LYR/Benchmark/scNT/scNT.csv')%>%
  t() %>% as.data.frame()

scNT$scNT_accurcy <- rowMeans(scNT)
scNT$Method <- rownames(scNT)

scNT <- scNT[2:15,c('Method','scNT_accurcy')]

scNT <-  check_method_name(scNT)

####--------time / scalability ----------------------------
metrics <- fread('/media/liyaru/LYR/Benchmark/time/time_new _1125.csv') %>%
  t() %>%
  as.data.frame()

colnames(metrics) <- metrics[1,]
metrics <- metrics[2:16,]
metrics$Method = rownames(metrics)

metrics <- check_method_name(metrics)

time <- metrics[,c('Method','1000','10000','50000')] %>% as.data.frame()

colnames(time) <- c('Method','1k cells','10k cells','50k cells')

# time[,2:4] <- -time[,2:4]
# time[,2:4] <- 1/((time[,2:4])/ 3600) # hours
time[,2:4] <- 1/(time[,2:4]) # seconds

time$time_mean <- apply(time[,2:3], 1,mean, na.rm = TRUE)

# time[,2:4] <- apply(time[,2:4],2,normalize)
x_min <- min(time[,2:4], na.rm = TRUE)-1
# x_min <- min(time[,2:4], na.rm = TRUE) - 0.1
# x_min <- 0
x_max <- max(time[,2:4], na.rm = TRUE)
time[,2:4] <- (time[,2:4] - x_min) / (x_max - x_min)

# time$time_mean <- normalize(time$time_mean)


####--------  merge data--------------

metrics <- cos_v[,c('Method','Mean')]
colnames(metrics)[2] <- c('Cosine_dyngen')

# add simulation celldancer
# metrics <- merge(metrics,simu_cd,by='Method',all.x = T,all.y = T)

# add HuBMAP
metrics <- merge(metrics,HuBMAP[,c('HuBMAP_GCoh_mean','Method')],by='Method',all.x = T,all.y = T)

# add datasets10
metrics <- merge(metrics,data10[,c('data10_GCoh_mean','Method')],by='Method',all.x = T,all.y = T)

# add redeem
metrics <- merge(metrics,redeem,by='Method',all.x = T)

# # add multitag
# metrics <- merge(metrics,multitag[,c('Method','Cor_with_lineage_mean')],by='Method',all.x = T,all.y = T)
# # add scNT
# metrics <- merge(metrics,scNT[,c('Method','Cor_with_metabolic_in_low_dim_mean')],by='Method',all.x = T,all.y = T)

# add time
metrics <- merge(metrics,time,by='Method',all.x = T,all.y = T)

# add noise
metrics <- merge(metrics,simu_noise[,c('Method',"noise_0.4","noise_1", "noise_2","Stability")],
                 by='Method',all.x = T,all.y = T)

metrics <- merge(metrics,tag,
                 by='Method',all.x = T,all.y = T)

metrics <- merge(metrics,scNT,
                 by='Method',all.x = T,all.y = T)
# normalize
# metrics[,2:ncol(metrics)] <- apply(metrics[,2:ncol(metrics)],2,normalize)
# metrics[, c('1k cells','10k cells','50k cells')] <- metrics[, c('1k cells','10k cells','50k cells')] / 6000
# metrics[,2:6] <- as.data.frame(scale(metrics[,2:6]))



metrics$Method <- factor(metrics$Method,levels = names(method_palette))
metrics <- metrics[order(metrics$Method),]
metrics = as.data.frame(metrics)


colnames(metrics) <- c(
"Method",
"Cosine similarity in simulation datasets1 (n=42)",
"GCoh in HuBMAP datasets (n=176)",
"GCoh in classic datasets (n=10)",
"GCoh in Redeem dataset",
"Computational efficiency (1k cells)",
"Computational efficiency (10k cells)",
"Computational efficiency (50k cells)",
"Scalability",
"Percentage of performance (noise level 0.4)",
"Percentage of performance (noise level 1.0)",
"Percentage of performance (noise level 2.0)",
"Stability",
"Cosine similarity in scLT dataset",
"Cosine similarity in scNT dataset")

acc = c("Cosine similarity in simulation datasets1 (n=42)",
        "Cosine similarity in scLT dataset",
        "Cosine similarity in scNT dataset",
        "GCoh in HuBMAP datasets (n=176)",
        "GCoh in classic datasets (n=10)",
        "GCoh in Redeem dataset")

x_min <- -1
x_max <- 1
metrics[,acc] <- (metrics[,acc] - x_min) / (x_max - x_min)

# metrics[,acc] <- apply(metrics[,acc],2,normalize) # normalize by col
metrics$Accuracy <- rowMeans(metrics[,acc], na.rm = TRUE)

# overall score
# metrics$Overall = 0.5*metrics$Accuracy + 0.25*metrics$Scalability + 0.25*metrics$Stability
# metrics[metrics$Method=='VeloAE','Overall'] = 0.75*metrics[metrics$Method=='VeloAE','Accuracy'] + 0.25*metrics[metrics$Method=='VeloAE','Scalability']
# metrics[metrics$Method=='TFvelo','Overall'] = 0.75*metrics[metrics$Method=='TFvelo','Accuracy'] + 0.25*metrics[metrics$Method=='TFvelo','Scalability']
# metrics[metrics$Method=='LatentVelo','Overall'] = 0.75*metrics[metrics$Method=='LatentVelo','Accuracy'] + 0.25*metrics[metrics$Method=='LatentVelo','Scalability']

Accuracy_rank_re = rank(metrics$Accuracy,na.last = 'keep')
Scalability_rank_re = rank(metrics$Scalability,na.last = 'keep')
Stability_rank_re = rank(metrics$Stability,na.last = 'keep')

# metrics$Overall = 0.6*metrics$Accuracy + 0.4*metrics$Scalability
metrics$`Overall Score` = (0.6*Accuracy_rank_re  + 0.4*rank(Scalability_rank_re)) / 14
# metrics$Overall = normalize(metrics$Overall)

# # show rank
metrics$`Accuracy rank index` = Accuracy_rank_re/ 14
metrics$`Scalability rank index`= Scalability_rank_re/14
metrics$`Stability rank index`= Stability_rank_re/ 14

new_order = c(
  "Method",
  'Overall Score',
  'Accuracy rank index',
  "Cosine similarity in simulation datasets1 (n=42)",
  "Cosine similarity in scLT dataset",
  "Cosine similarity in scNT dataset",
  "GCoh in HuBMAP datasets (n=176)",
  "GCoh in classic datasets (n=10)",
  "GCoh in Redeem dataset",
  "Scalability rank index",
  "Computational efficiency (1k cells)",
  "Computational efficiency (10k cells)",
  "Computational efficiency (50k cells)",
  "Stability rank index",
  "Percentage of performance (noise level 0.4)",
  "Percentage of performance (noise level 1.0)",
  "Percentage of performance (noise level 2.0)"
  )
metrics <- metrics[,new_order]

# # normalize by group
# x_min = -1
# x_max = 1
# metrics[,2:10] <- (metrics[,2:10] - x_min) / (x_max - x_min)


####----- basic info ------------
row_height <- 1.1
row_space <- .1
row_bigspace <- .5
col_width <- 1.1
col_space <- .2
col_bigspace <- .5
segment_data <- NULL

data =as.data.frame(metrics)  

####----- row info and pos--------
row_info <- as.data.frame(metrics$Method)
# colnames(row_info) <- 'Method'
colnames(row_info) <- 'id'
row_pos <- 
  row_info %>%
  mutate(
    group_i = row_number(),
    row_i = row_number(),
    colour_background = group_i %% 2 == 1,
    # do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
    do_spacing = FALSE,
    ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
    y = - (row_i * row_height + cumsum(ysep)),
    ymin = y - row_height / 2,
    ymax = y + row_height / 2
  )

####----- col info and pos--------
column_info <- data.frame(id = colnames(metrics))

column_info$group <- c('Text',
                       'Overall',
                       rep("Accuracy",7),
                       # rep("Accuracy2",3),
                       rep("Scalability",4),
                       rep('Stability',4))

column_info$geom <- c('text',
                      'bar','bar',
                      rep('circle',6),
                      'bar',
                      rep('circle',3),
                      'bar',
                      rep('circle',3)
)

column_info$width <- c(6,
                       2,2,
                       rep(1,6),
                       2,
                       rep(1,3),
                       2,
                       rep(1,3))

column_info$overlay <- FALSE

column_pos <-
  column_info %>%
  mutate(
    do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
    xsep = case_when(
      overlay ~ c(0, -head(width, -1)),
      do_spacing ~ col_bigspace,
      TRUE ~ col_space
    ),
    xwidth = case_when(
      overlay & width < 0 ~ width - xsep,
      overlay ~ -xsep,
      TRUE ~ width
    ),
    xmax = cumsum(xwidth + xsep),
    xmin = xmax - xwidth,
    x = xmin + xwidth / 2
  )

####-----circle data------------
# gather circle data
ind_circle <- which(column_info$geom == "circle")

dat_mat <- as.matrix(data[, ind_circle])
col_palette <- data.frame(metric = colnames(dat_mat), 
                          group = column_info[match(colnames(dat_mat), column_info$id), "group"])

# col_palette$name_palette <- 'RdPu'
col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])

circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                function(x) rep(x, nrow(dat_mat)))), 
                          x0 = unlist(lapply(column_pos$x[ind_circle], 
                                             function(x) rep(x, nrow(dat_mat)))), 
                          y0 = rep(row_pos$y, ncol(dat_mat)),
                          r = row_height/2*as.vector(sqrt(dat_mat))
)

# rescale by col
# for(l in unique(circle_data$label)){
#   ind_l <- which(circle_data$label == l)
#   circle_data[ind_l, "r"] <- rescale(circle_data[ind_l, "r"], to = c(0.05, 0.55), from = range(circle_data[ind_l, "r"], na.rm = T))
# }

colors <- NULL


for(i in 1:ncol(dat_mat)){
  palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
  # palette <- colorRampPalette(brewer.pal(9, col_palette$name_palette[[i]]))(nrow(data)-sum(is.na(dat_mat[,i])))
  colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
}

circle_data$colors <- colors

####-------bar/rect data-----------
ind_bar <- which(column_info$geom == "bar")
dat_mat <- as.matrix(data[, ind_bar])
# colnames(dat_mat) <- 'Accuracy'

col_palette <- data.frame(metric = colnames(dat_mat), 
                          group = column_info[match(colnames(dat_mat), column_info$id), "group"])


col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
# col_palette$name_palette <- 'BuPu'

rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                              function(x) rep(x, nrow(dat_mat)))),
                        method = rep(row_info$id, ncol(dat_mat)),
                        value = as.vector(dat_mat),
                        xmin = unlist(lapply(column_pos[ind_bar, "xmin"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        xmax = unlist(lapply(column_pos[ind_bar, "xmax"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        ymin = rep(row_pos$ymin, ncol(dat_mat)),
                        ymax = rep(row_pos$ymax, ncol(dat_mat)),
                        xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"], 
                                               function(x) rep(x, nrow(dat_mat))))
)

rect_data <- rect_data %>%
  add_column_if_missing(hjust = 0) %>%
  mutate(
    xmin = xmin + (1 - value) * xwidth * hjust,
    xmax = xmax - (1 - value) * xwidth * (1 - hjust)
  )

colors <- NULL

for(i in 1:ncol(dat_mat)){
  palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
  # palette <- colorRampPalette(brewer.pal(9, col_palette$name_palette[[i]]))(nrow(data)-sum(is.na(dat_mat[,i])))
  colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
}

rect_data$colors <- colors

####-------text data------------
ind_text <- which(column_info$geom == "text")
dat_mat <- as.matrix(data[, ind_text])
colnames(dat_mat) <- 'Method'

text_data <- data.frame(label_value = as.vector(dat_mat), 
                        group = rep(colnames(dat_mat), each = nrow(dat_mat)),
                        xmin = unlist(lapply(column_pos[ind_text, "xmin"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        xmax = unlist(lapply(column_pos[ind_text, "xmax"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        ymin = rep(row_pos$ymin, ncol(dat_mat)),
                        ymax = rep(row_pos$ymax, ncol(dat_mat)),
                        size = 4, fontface = "plain", stringsAsFactors = F)

#####-------- add top3 ranking for each bar column-------------

cols_bar <- unique(rect_data$label)
cols_bar <- as.character(cols_bar[!is.na(cols_bar)])
for(c in cols_bar){
  rect_tmp <- rect_data[rect_data$label == c,]
  rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(-rect_tmp$value, ties.method = "min")))
  # rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(rect_tmp$value, ties.method = "min")))
  rect_tmp <- rect_tmp[rect_tmp$label_value %in% c("1", "2", "3","4",'5','6'), c("label_value", "xmin", "xmax", "ymin", "ymax")]
  rect_tmp <- add_column(rect_tmp, "size" = 2.5, .after = "ymax")
  rect_tmp <- add_column(rect_tmp, "colors" = "black", .after = "size")
  rect_tmp <- add_column(rect_tmp, "fontface" = "plain", .after = "colors")
  rect_tmp <- add_column(rect_tmp, "group" = "top3", .after = "fontface")
  text_data <- bind_rows(text_data, rect_tmp)
}



####-----add column names-------------
# ADD COLUMN NAMES
df <- column_pos %>% filter(id != "Method") %>% filter(id != "Ranking")

if (nrow(df) > 0) {
  segment_data <- segment_data %>% bind_rows(
    df %>% transmute(x = x, xend = x, y = -.3, yend = -.1, size = .5)
  )
  text_data <-
    bind_rows(
      text_data,
      df %>% transmute(
        xmin = x, xmax = x, ymin = 0, ymax = -0.5,
        angle = 30, vjust = 0, hjust = 0,
        label_value = id, 
        size = 3
      )
    )
}


####----------row annotation---------------
# # GENERATE ROW ANNOTATION
# if (plot_row_annotation) {
#   row_annotation <-
#     row_pos %>% 
#     select(group, ymin, ymax) %>%
#     group_by(group) %>%
#     summarise(
#       ymin = min(ymin),
#       ymax = max(ymax),
#       y = (ymin + ymax) / 2
#     ) %>%
#     ungroup() %>%
#     mutate(xmin = -.5, xmax = 5) %>%
#     filter(!is.na(group), group != "")
#   
#   text_data <- text_data %>% bind_rows(
#     row_annotation %>%
#       transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .), 
#                 hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
#       mutate(ymax = ymin + row_height)
#   )
# }
# 
# ####----------image data---------------
# # gather image data
# ind_img <- which(column_info$geom == "image")
# if(length(ind_img) > 0){
#   dat_mat <- as.matrix(data[, ind_img])
#   
#   image_data <- data.frame(x = unlist(lapply(column_pos$x[ind_img], 
#                                              function(x) rep(x, nrow(dat_mat)))), 
#                            y = rep(row_pos$y, ncol(dat_mat)),
#                            image = mapvalues(dat_mat, from = c("graph", "embed", "gene"), 
#                                              to = c("./img/graph.png", "./img/embedding.png", "./img/matrix.png")),
#                            stringsAsFactors = FALSE
#   )
#   
# }
# 
# suppressWarnings({
#   minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend, 
#                    text_data$xmin, na.rm = TRUE)
#   maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend, 
#                    text_data$xmax, na.rm = TRUE)
#   minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,  
#                    text_data$ymin, na.rm = TRUE)
#   maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend, 
#                    text_data$ymax, na.rm = TRUE)
# })

####---- legend---------
# CREATE LEGEND for circle scores

minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend, 
                 text_data$xmin, na.rm = TRUE)
maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend, 
                 text_data$xmax, na.rm = TRUE)
minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,  
                 text_data$ymin, na.rm = TRUE)
maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend, 
                 text_data$ymax, na.rm = TRUE)

x_min_output <- minimum_x+0.5
x_min_scaling <- minimum_x + 5.5
x_min_ranking <- minimum_x + 10.5
x_min_score <-  minimum_x + 17

leg_max_y <- minimum_y - .5

# # Create legend for Output
# leg_min_x <- x_min_output
# output_title_data <- data.frame(xmin = leg_min_x, 
#                                 xmax = leg_min_x+ 2, 
#                                 ymin = leg_max_y - 1, 
#                                 ymax = leg_max_y, 
#                                 label_value = "Output", 
#                                 hjust = 0, vjust = 0, 
#                                 fontface = "bold",
#                                 size = 3)
# 
# output_img <- data.frame(x = leg_min_x+0.5,
#                          y = c(leg_max_y-2, leg_max_y-3.2,leg_max_y-4.4),
#                          image = c("./img/matrix.png", "./img/embedding.png", "./img/graph.png")
# )
# if(atac || atac_best){
#   output_text <- data.frame(xmin = leg_min_x+1.5, 
#                             xmax = leg_min_x+3, 
#                             ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
#                             ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
#                             label_value = c("feature", "embed", "graph"), 
#                             hjust = 0, vjust = 0, 
#                             fontface = "plain",
#                             size = 3)
# } else{
#   output_text <- data.frame(xmin = leg_min_x+1.5, 
#                             xmax = leg_min_x+3, 
#                             ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
#                             ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
#                             label_value = c("gene", "embed", "graph"), 
#                             hjust = 0, vjust = 0, 
#                             fontface = "plain",
#                             size = 3)
# }
# 
# text_data <- bind_rows(text_data, output_text, output_title_data)
# image_data <- bind_rows(image_data, output_img)
# 
# # Create legend for scaling
# if(!atac && !atac_best){
#   leg_min_x <- x_min_scaling
#   scaling_title_data <- data.frame(xmin = leg_min_x, 
#                                    xmax = leg_min_x+ 2, 
#                                    ymin = leg_max_y - 1, 
#                                    ymax = leg_max_y, 
#                                    label_value = "Scaling", 
#                                    hjust = 0, vjust = 0, 
#                                    fontface = "bold",
#                                    size = 3)
#   
#   scaling_text <- data.frame(xmin = c(leg_min_x, leg_min_x+1), 
#                              xmax = c(leg_min_x+0.5, leg_min_x+3), 
#                              ymin = c(rep(leg_max_y-2,2), rep(leg_max_y-3,2)), 
#                              ymax = c(rep(leg_max_y-1,2), rep(leg_max_y-2,2)), 
#                              label_value = c("+", ": scaled", "-", ": unscaled"), 
#                              hjust = 0, vjust = 0, 
#                              fontface = c("bold","plain", "bold", "plain"),
#                              size = c(5,3,5,3))
#   
#   text_data <- bind_rows(text_data, scaling_title_data, scaling_text)
# }

# CREATE LEGEND for ranking colors

####-----rank legend-----------
# CREATE LEGEND for ranking colors
leg_min_x <- x_min_ranking
# rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])
# rank_groups <- as.character(column_info[column_info$geom == "circle", "group"])
rank_groups <- as.character(column_info[column_info$geom %in% c("circle",'bar'), "group"])


rank_minimum_x <- list('Overall' = leg_min_x,
                       "Accuracy" =leg_min_x+1, 
                       "Accuracy2" =leg_min_x+2, 
                       "Scalability" = leg_min_x+3, 
                       "Stability" = leg_min_x+4)

leg_max_x <- leg_min_x+4


rank_title_data <- data.frame(xmin = leg_min_x, 
                              xmax = leg_min_x+ 2, 
                              ymin = leg_max_y - 1, 
                              ymax = leg_max_y, 
                              label_value = "Ranking", 
                              hjust = 0, vjust = 0, 
                              fontface = "bold")


for(rg in rank_groups){
  # rg = rank_groups[2]
  rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)
  # rank_palette <- colorRampPalette(brewer.pal(9, palettes[[rg]]))(5)
  
  rank_data <- data.frame(xmin = rank_minimum_x[[rg]],
                          xmax = rank_minimum_x[[rg]] + .8,
                          ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                          ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                          border = TRUE,
                          colors = rank_palette
  )
  rect_data <- bind_rows(rect_data, rank_data)
  
}

# create arrow for ranking
arrow_data <- data.frame(x = leg_max_x + 1.5, 
                         xend = leg_max_x +1.5, 
                         y = leg_max_y-4, 
                         yend = leg_max_y -1.5)


# add text next to the arrow
arrow_text <- data.frame(xmin = leg_max_x +2, 
                         xmax = leg_max_x +2.5, 
                         ymin = c(leg_max_y-2, leg_max_y-4), 
                         ymax = c(leg_max_y-1.5, leg_max_y-3.5 ), 
                         label_value = c("1", as.character(nrow(data))), 
                         hjust = 0, vjust = 0, size = 2.5)


text_data <- bind_rows(text_data, rank_title_data, arrow_text)

####---------circle legend--------
cir_minimum_x <- x_min_score

cir_legend_size <- 1
cir_legend_space <- .1

cir_legend_dat <-
  data.frame(
    value = seq(0, 1, by = .2),
    r = row_height/2*seq(0, 1, by = .2)
  )
cir_legend_dat$r <- rescale(cir_legend_dat$r, to = c(0.05, 0.55), from = range(cir_legend_dat$r, na.rm = T))

x0 <- vector("integer", nrow(cir_legend_dat))
for(i in 1:length(x0)){
  if(i == 1){
    x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
  }
  else {
    x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
  }
}

cir_legend_dat$x0 <- x0+0.5
cir_legend_min_y <- leg_max_y-4
cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r

cir_legend_dat$colors <- NULL
cir_maximum_x <- max(cir_legend_dat$x0)

cir_title_data <- data_frame(xmin = cir_minimum_x, 
                             xmax = cir_maximum_x, 
                             ymin = leg_max_y -1, 
                             ymax = leg_max_y,
                             label_value = "Score", 
                             hjust = 0, vjust = 0, fontface = "bold")

cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                             xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                             ymin = cir_legend_min_y,
                             ymax = cir_legend_min_y +3,
                             hjust = .5, vjust = 0, size = 2.5,
                             label_value = ifelse(cir_legend_dat$value %in% c(0, 1), 
                                                  paste0(cir_legend_dat$value*100, "%"), ""))

circle_data <- bind_rows(circle_data, cir_legend_dat)
text_data <- bind_rows(text_data, cir_title_data, cir_value_data)

minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))

####----GGPLOT ------------
g <-
  ggplot() +
  coord_equal(expand = FALSE) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_size_identity() +
  scale_linetype_identity() +
  cowplot::theme_nothing()

# PLOT ROW BACKGROUNDS
df <- row_pos %>% filter(colour_background)

if (nrow(df) > 0) {
  g <- g + 
    geom_rect(aes(xmin = min(column_pos$xmin)-.25, 
                  xmax = max(column_pos$xmax)+.25, 
                  ymin = ymin - (row_space / 2), 
                  ymax = ymax + (row_space / 2)), 
              df, fill = "#DDDDDD")
} 
g

####------ PLOT CIRCLES ------------
if (length(ind_circle) > 0) {
  g <- g + 
    ggforce::geom_circle(
      aes(x0 = x0, y0 = y0, fill= colors, r = r), 
      circle_data, 
      size=.25)
}
g

####--- PLOT RECTANGLES ------------
if (nrow(rect_data) > 0) {
  # add defaults for optional values
  rect_data <- rect_data %>%
    add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
    mutate(border_colour = ifelse(border, border_colour, NA))
  
  g <- g + 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin,ymax = ymax, 
                  fill = colors, colour = border_colour, alpha = alpha), 
              rect_data, 
              size = .25)
}
g

####------ PLOT TEXT------------------

text_data <- text_data %>%
  add_column_if_missing(
    hjust = .5,
    vjust = .5,
    size = 3,
    fontface = "plain",
    colors = "black",
    lineheight = 1,
    angle = 0
  ) %>%
  mutate(
    angle2 = angle / 360 * 2 * pi,
    cosa = cos(angle2) %>% round(2),
    sina = sin(angle2) %>% round(2),
    alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
    alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
    x = (1 - alphax) * xmin + alphax * xmax,
    y = (1 - alphay) * ymin + alphay * ymax
  ) %>%
  filter(label_value != "")

# Set fontface for legend bold
text_data[text_data$label_value == "Ranking", "fontface"] <- "bold"

# Set fontface for ranking numbers bold
# if(usability || atac_best){
#   text_data[1:nrow(data), "fontface"] <- "bold"
# }
# subset text_data to left-aligned rows

text_data_left <- text_data[which(text_data$group == "Method" | text_data$group == "top3"), ]
text_data <- text_data[-which(text_data$group == "Method" | text_data$group == "top3"), ]

g <- g + geom_text(aes(x = x, y = y, label = label_value, 
                       colour = colors, hjust = hjust, vjust = vjust, 
                       size = size, fontface = fontface, angle = angle), 
                   data = text_data)
g

text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] - 3

# if(usability || atac_best){
#   text_data_left[text_data_left$group == "top3", "x"] <- text_data_left[text_data_left$group == "top3", "xmin"] + .3
#   text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] + .5
# }

g <- g + geom_text(aes(x = x, y = y, 
                       label = label_value, colour = colors, hjust = "left", 
                       vjust = vjust, size = size, fontface = fontface, angle = angle), 
                   data = text_data_left)
g


g <- g + 
  geom_text(aes(x = x, y = y, label = label_value, colour = colors, 
                hjust = "left", vjust = vjust, size = size, 
                fontface = fontface, angle = angle), 
            data = text_data_left)
g


#### -------- plot segments-----------------
# PLOT SEGMENTS
if (nrow(segment_data) > 0) {
  # add defaults for optional values
  segment_data <- segment_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
  
  g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), segment_data)
}

# PLOT ARROW RANKING
if (nrow(arrow_data) > 0) {
  # add defaults for optional values
  arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
  
  g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
}

# # PLOT IMAGES
# if(length(ind_img) > 0){
#   for(r in 1:nrow(image_data)){
#     g <- g + cowplot::draw_image(image = image_data$image[r], x = image_data[r, "x"]-.5, y = image_data[r, "y"]-.5)
#   }
#   
# }


####------  add  size------------
# reserve a bit more room for text that wants to go outside the frame

minimum_x <- minimum_x - 2
maximum_x <- maximum_x + 5
minimum_y <- minimum_y - 2
maximum_y <- maximum_y + 10

g$width <- maximum_x - minimum_x
g$height <- maximum_y - minimum_y

g <- g + 
  expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))

g

ggsave("/media/liyaru/LYR/Benchmark/Fig/3.score_merge.pdf",
       g, 
       device = "pdf", width = 400, height = 300, units = "mm")
