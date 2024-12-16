setwd(dir = "/media/zx/HDD1/cooperation/liyr/result/dyngen/result_0921/eval/gene_level/")

name_para <- c("Global","Incluster","Direction")
name_ground <- c("velocity","velocity_pca","trans_matrix")
simuID_list = c('bifurcating_converging_seed1', 'bifurcating_converging_seed2','bifurcating_converging_seed3','bifurcating_cycle_seed1', 'bifurcating_cycle_seed2', 'bifurcating_cycle_seed3' ,'bifurcating_loop_seed1', 'bifurcating_loop_seed2', 'bifurcating_loop_seed3' ,'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1' ,'binary_tree_seed2', 'binary_tree_seed3', 'branching_seed1', 'branching_seed2', 'branching_seed3', 'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3', 'converging_seed1', 'converging_seed2','converging_seed3', 'cycle_seed1', 'cycle_seed2', 'cycle_seed3' ,'cycle_simple_seed1' ,'cycle_simple_seed2', 'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2', 'disconnected_seed3' ,'linear_seed1' ,'linear_seed2', 'linear_seed3', 'linear_simple_seed1' ,'linear_simple_seed2', 'linear_simple_seed3', 'trifurcating_seed1' ,'trifurcating_seed2' ,'trifurcating_seed3' )
simu_group = unique(x = sapply(X = simuID_list,FUN = function(x)substr(x = x,start = 1,stop = nchar(x)-6)))


########################### 
library(corrplot)

file_all <- c("eval_global.csv","eval_incluster.csv","eval_direction.csv","cosine_sim_velocity.T.csv","cosine_sim_velocity_pca2.T.csv","cosine_sim_trans_matrix.T.csv")
name_all <- c("Global","Incluster","Direction","velocity","velocity_pca","trans_matrix")

simuID_list = c('bifurcating_converging_seed1', 'bifurcating_converging_seed2','bifurcating_converging_seed3','bifurcating_cycle_seed1', 'bifurcating_cycle_seed2', 'bifurcating_cycle_seed3' ,'bifurcating_loop_seed1', 'bifurcating_loop_seed2', 'bifurcating_loop_seed3' ,'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1' ,'binary_tree_seed2', 'binary_tree_seed3', 'branching_seed1', 'branching_seed2', 'branching_seed3', 'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3', 'converging_seed1', 'converging_seed2','converging_seed3', 'cycle_seed1', 'cycle_seed2', 'cycle_seed3' ,'cycle_simple_seed1' ,'cycle_simple_seed2', 'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2', 'disconnected_seed3' ,'linear_seed1' ,'linear_seed2', 'linear_seed3', 'linear_simple_seed1' ,'linear_simple_seed2', 'linear_simple_seed3', 'trifurcating_seed1' ,'trifurcating_seed2' ,'trifurcating_seed3' )
file0 <- file_all[1]
data0 <- read.csv(file = paste0("genes/",file0),header = T,row.names = 1)
data0 <- data0[-11,]
data_use <- data.frame("Dataset"=rep(simuID_list, rep(13,42)),"Method" = rep(rownames(data0),42))
data_use$Group <- sapply(X =data_use$Dataset,FUN = function(x)substr(x = x,start = 1,stop = nchar(x)-6))
for (i in 1:6){
    name0 <- name_all[i]
    file0 <- file_all[i]
    data0 <- read.csv(file = paste0("genes/",file0),header = T,row.names = 1)
    data0 <- data0[-11,]
    data_use[[name0]] <- unlist(data0)
}

library("Hmisc")
library(corrplot)
library(cowplot)
library(ggplotify)
library(gridGraphics)
library(scater)

list_data <- list()
list_cor <- list()
list_pvalue <- list()
# sum_cor <- data.frame(rep())
for (i in unique(x = data_use$Group)){
    list_data[[i]] <- data_use[data_use$Group == i,]
    rownames(list_data[[i]]) <- paste0(list_data[[i]]$Dataset, "-", list_data[[i]]$Method)
    list_data[[i]] <- list_data[[i]][,4:9]
    list_cor[[i]] <- rcorr(as.matrix(list_data[[i]]))$r
    list_pvalue[[i]] <- rcorr(as.matrix(list_data[[i]]))$P
    list_pvalue[[i]][is.na(list_pvalue[[i]])] <- 0
}

list_pic <- list()
for (i in unique(x = data_use$Group)){
    print(i)
    corrplot(list_cor[[i]], type="upper",p.mat = list_pvalue[[i]],col = COL2('RdBu', 200)[200:1], insig = "blank", title = i,mar = c(1,1,1,1))
    grid.echo()
    P1 <- grid.grab()
    P1 <- editGrob(P1,gPath("square"), grep = TRUE,gp = gpar(col = NA,fill = NA))
    #P1 <- editGrob(P1,gPath("circle"), grep = TRUE,gp = gpar(col = NA,fill = NA))
    # P1 <- editGrob(P1,gPath("symbols-rect-1"), grep = TRUE,gp = gpar(fill = matrix.colors))
    list_pic[[i]] <- editGrob(P1,gPath("background"), grep = TRUE,gp = gpar(fill = NA))
}
pdf(file = "corplot_split_dataset.pdf",width = 14*9,height = 7)
multiplot(plotlist = list_pic,cols = 14)
dev.off()

#################################
data_use2 = data_use
# data_use2 = data_use[!data_use$Group %in% c("bifurcating_cycle","bifurcating_loop"),]
cor_data <- rcorr(as.matrix(data_use2[,4:9]))$r
cor_data

pvalue_data <- rcorr(as.matrix(data_use2[,4:9]))$P
pvalue_data[is.na(pvalue_data)] <- 0
pdf(file = "corplot_all_dataset.pdf",width = 9,height = 7)
corrplot(cor_data, type="upper",p.mat = pvalue_data,col = COL2('RdBu', 200)[200:1],sig.level = 0.5,  insig = "blank",mar = c(1,1,1,1))
dev.off()
pdf(file = "corplot_all_dataset_mixed.pdf",width = 9,height = 7)
corrplot.mixed(cor_data,upper.col = COL2('RdBu', 200)[200:1],lower.col = COL2('RdBu', 200)[200:1], insig = "blank",mar = c(1,1,1,1))
dev.off()



#########################################################

#########################33 cosine
library(openxlsx)
library(lsa)

list_res <- list()

for (i in simuID_list){
    
    bb <- read.xlsx(xlsxFile = "/media/zx/HDD1/cooperation/liyr/result/dyngen/result_0921/eval/gene_level/genes//cosine_sim_velocity_pca2-30.T.xlsx",sheet = i,rowNames = T)
    #colnames(bb) <- method_names
    bb <- bb[,-11]
    list_res[[i]] <- bb
}

list_res2 <- list()
for (i in 1:29){
    list_tmp <- list()
    for (j in simuID_list){
        tt <- unlist(list_res[[j]][i,])
        list_tmp[[j]] <- tt
    }
    list_tmp2 <- data.frame(list_tmp)
    list_res2[[i]] <- list_tmp2
}

data_velocity <- read.csv(file = "genes/cosine_sim_velocity.T.csv",header = T,row.names = 1)
data_velocity <- data_velocity[-c(11,14),]



xx <- list()

for (g in simu_group){
    simuID_list_tmp <- grep(pattern = g,x = simuID_list,value = T)
    tt <- c() 
    for (i in 1:29){
        xxx <- unlist(list_res2[[i]][,simuID_list_tmp])
        xxx <- xxx[!is.na(xxx)]
        yyy <- unlist(data_velocity[,simuID_list_tmp])
        yyy <- yyy[!is.na(yyy)]
        tt <- c(tt,cosine(x = xxx,y = yyy)[1,1])
    }
    xx[[g]] <- tt
}


simuID_list_tmp <- simuID_list
tt <- c() 
for (i in 1:29){
    xxx <- unlist(list_res2[[i]][,simuID_list_tmp])
    xxx <- xxx[!is.na(xxx)]
    yyy <- unlist(data_velocity[,simuID_list_tmp])
    yyy <- yyy[!is.na(yyy)]
    tt <- c(tt,cosine(x = xxx,y = yyy)[1,1])
}
xx[["total"]] <- tt

xx <- data.frame(xx)
rownames(xx) <- paste0("dim_",2:30)
data_plot <- data.frame("Dim" = rep(2:30,ncol(xx)),"Dataset" = rep(colnames(xx),rep(nrow(xx),ncol(xx))),"Cosine_sim" = unlist(xx))
data_plot$Dim <- factor(x = data_plot$Dim,levels = 2:30)
data_plot$Dataset <- factor(x = data_plot$Dataset,levels = c(simu_group,"total"))

library(scales)
#color_model <- rainbow(14)
color_model <- hue_pal()(14)
color_model <- c(color_model,"black")
pdf(file = "velocity_pca_2-30_split_datasets_ALL_cosine_similarity.pdf",width = 10,height = 6)
ggplot(data = data_plot,mapping = aes(x = Dim,y = Cosine_sim,group = Dataset,color = Dataset))+
    geom_point()+
    geom_line()+
    scale_color_manual(values = color_model)+
    ylab(label = "")+
    ggtitle(label = paste0("Cosine similarity"))+
    theme_bw()+theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))
dev.off()



###################################### cosine ground_truth
xx <- list()

for (g in simu_group){
    simuID_list_tmp <- grep(pattern = g,x = simuID_list,value = T)
    tt <- c() 
    for (i in 1:29){
        xxx <- unlist(list_res2[[i]][,simuID_list_tmp])
        xxx <- xxx[!is.na(xxx)]
        tt <- c(tt,median(xxx))
    }
    xx[[g]] <- tt
}


simuID_list_tmp <- simuID_list
tt <- c() 
for (i in 1:29){
    xxx <- unlist(list_res2[[i]][,simuID_list_tmp])
    xxx <- xxx[!is.na(xxx)]
    tt <- c(tt,median( xxx))
}
xx[["total"]] <- tt

xx <- data.frame(xx)
rownames(xx) <- paste0("dim_",2:30)
data_plot <- data.frame("Dim" = rep(2:30,ncol(xx)),"Dataset" = rep(colnames(xx),rep(nrow(xx),ncol(xx))),"Cosine_sim" = unlist(xx))
data_plot$Dim <- factor(x = data_plot$Dim,levels = 2:30)
data_plot$Dataset <- factor(x = data_plot$Dataset,levels = c(simu_group,"total"))

library(scales)
color_model <- rainbow(14)
color_model <- hue_pal()(14)
color_model <- c(color_model,"black")
pdf(file = "velocity_pca_2-30_split_datasets_ALL_cosine_similarity_with_ground_truth.pdf",width = 10,height = 6)
ggplot(data = data_plot,mapping = aes(x = Dim,y = Cosine_sim,group = Dataset,color = Dataset))+
    geom_point()+
    geom_line()+
    scale_color_manual(values = color_model)+
    ylab(label = "")+
    ggtitle(label = paste0("Cosine similarity"))+
    theme_bw()+theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))
dev.off()

