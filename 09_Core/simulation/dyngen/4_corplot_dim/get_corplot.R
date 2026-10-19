setwd(dir = "/media/zx/HDD1/cooperation/liyr/result/dyngen/result_0921/eval/res1/")
library("Hmisc")
library(corrplot)
library(cowplot)
library(ggplotify)
library(gridGraphics)
##################################################
file_all1 <- c("eval_global.csv","eval_incluster.csv","eval_direction.csv")
file_all2 <- c("eval_global_median.csv","eval_incluster_median.csv","eval_direction_median.csv")

file_all4 <- c("cosine_sim_velocity.csv","cosine_sim_velocity_pca2.csv","cosine_sim_trans_matrix.csv")
file_all5 <- c("cosine_sim_velocity.T.csv","cosine_sim_velocity_pca2.T.csv","cosine_sim_trans_matrix.T.csv")


name_all1 <- c("GCoh","ICCoh","CBDir")
name_all2 <- c("Cosine similarity on high-dim velocity","Cosine similarity on low-dim velocity","Cosine similarity on transition probability")


simuID_list = c('bifurcating_converging_seed1', 'bifurcating_converging_seed2','bifurcating_converging_seed3','bifurcating_cycle_seed1', 'bifurcating_cycle_seed2', 'bifurcating_cycle_seed3' ,'bifurcating_loop_seed1', 'bifurcating_loop_seed2', 'bifurcating_loop_seed3' ,'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1' ,'binary_tree_seed2', 'binary_tree_seed3', 'branching_seed1', 'branching_seed2', 'branching_seed3', 'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3', 'converging_seed1', 'converging_seed2','converging_seed3', 'cycle_seed1', 'cycle_seed2', 'cycle_seed3' ,'cycle_simple_seed1' ,'cycle_simple_seed2', 'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2', 'disconnected_seed3' ,'linear_seed1' ,'linear_seed2', 'linear_seed3', 'linear_simple_seed1' ,'linear_simple_seed2', 'linear_simple_seed3', 'trifurcating_seed1' ,'trifurcating_seed2' ,'trifurcating_seed3' )

######################################### read data cell level
data_use <- data.frame("Dataset"=rep(simuID_list, rep(13,42)),"Method" = rep(rownames(data0),42))
data_use$Group <- sapply(X =data_use$Dataset,FUN = function(x)substr(x = x,start = 1,stop = nchar(x)-6))


for (i in 1:3){
    name0 <- name_all1[i]
    file0 <- file_all1[i]
    data0 <- read.csv(file = paste0("median_all/",file0),header = T,row.names = 1)
    data0 <- data0[-11,]
    data_use[[name0]] <- unlist(data0)
}
for (i in 1:3){
    name0 <- name_all2[i]
    file0 <- file_all4[i]
    data0 <- read.csv(file = paste0("data/",file0),header = T,row.names = 1)
    data0 <- data0[-11,]
    data_use[[name0]] <- unlist(data0)
}

write.csv(x = data_use,file = "median_all/cell_level.csv")

#########################################3 read data gene level
data_use <- data.frame("Dataset"=rep(simuID_list, rep(13,42)),"Method" = rep(rownames(data0),42))
data_use$Group <- sapply(X =data_use$Dataset,FUN = function(x)substr(x = x,start = 1,stop = nchar(x)-6))


for (i in 1:3){
    name0 <- name_all1[i]
    file0 <- file_all1[i]
    data0 <- read.csv(file = paste0("median_all/",file0),header = T,row.names = 1)
    data0 <- data0[-11,]
    data_use[[name0]] <- unlist(data0)
}
for (i in 1:3){
    name0 <- name_all2[i]
    file0 <- file_all5[i]
    data0 <- read.csv(file = paste0("data/",file0),header = T,row.names = 1)
    data0 <- data0[-11,]
    data_use[[name0]] <- unlist(data0)
}

write.csv(x = data_use,file = "median_all/gene_level.csv")


##################### cor calculate
aaa = "median_all"
bbb = "cell"

data_use <- read.csv(file = paste0(aaa,"/",bbb,"_level.csv"),header = T,row.names = 1)
colnames(data_use)[4:9] <- c(name_all1,name_all2)
list_data <- list()
list_cor <- list()
list_pvalue <- list()
# sum_cor <- data.frame(rep())
for (i in unique(x = data_use$Group)){
    list_data[[i]] <- data_use[data_use$Group == i,]
    rownames(list_data[[i]]) <- paste0(list_data[[i]]$Dataset, "-", list_data[[i]]$Method)
    list_data[[i]] <- list_data[[i]][,4:9]
    #list_data[[i]] <- list_data[[i]][rowSums(is.na(list_data[[i]]))==0,]
    #list_cor[[i]] <- lsa::cosine(as.matrix(list_data[[i]]))
    #list_cor[[i]] <- rcorr(as.matrix(list_data[[i]]))$r
    list_cor[[i]] <- rcorr(as.matrix(list_data[[i]]),type = "pearson")$r
    list_pvalue[[i]] <- rcorr(as.matrix(list_data[[i]]))$P
    list_pvalue[[i]][list_pvalue[[i]] != 0] <- 0
    list_pvalue[[i]][is.na(list_pvalue[[i]])] <- 0
}

library(RColorBrewer)
library(scater)
##COL2('RdBu', 200)[200:1]
color_panel = colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100)
list_pic <- list()
for (i in unique(x = data_use$Group)){
    corrplot(list_cor[[i]], type="upper",p.mat = list_pvalue[[i]],col = color_panel,tl.col = "black",outline = T,diag=F,  insig = "blank",title = i,mar = c(1,1,1,1))
    grid.echo()
    P1 <- grid.grab()
    P1 <- editGrob(P1,gPath("square"), grep = TRUE,gp = gpar(col = NA,fill = NA))
    #P1 <- editGrob(P1,gPath("circle"), grep = TRUE,gp = gpar(col = NA,fill = NA))
    # P1 <- editGrob(P1,gPath("symbols-rect-1"), grep = TRUE,gp = gpar(fill = matrix.colors))
    list_pic[[i]] <- editGrob(P1,gPath("background"), grep = TRUE,gp = gpar(fill = NA))
}
pdf(file = paste0(aaa,"/corplot_split_group_",bbb,"_level.pdf"),width = 14*9,height = 7)
multiplot(plotlist = list_pic,cols = 14)
dev.off()

list_pic <- list()
for (i in unique(x = data_use$Group)){
    corrplot.mixed(list_cor[[i]],p.mat = list_pvalue[[i]],upper.col = color_panel,lower.col = color_panel, tl.col = "black",outline = T,diag= "u", insig = "blank",title = i,mar = c(1,1,1,1))
    grid.echo()
    P1 <- grid.grab()
    P1 <- editGrob(P1,gPath("square"), grep = TRUE,gp = gpar(col = NA,fill = NA))
    #P1 <- editGrob(P1,gPath("circle"), grep = TRUE,gp = gpar(col = NA,fill = NA))
    # P1 <- editGrob(P1,gPath("symbols-rect-1"), grep = TRUE,gp = gpar(fill = matrix.colors))
    list_pic[[i]] <- editGrob(P1,gPath("background"), grep = TRUE,gp = gpar(fill = NA))
}
pdf(file = paste0(aaa,"/corplot_mixed_split_group_",bbb,"_level.pdf"),width = 14*9,height = 7)
multiplot(plotlist = list_pic,cols = 14)
dev.off()



############
cor_data2 <- list_cor$bifurcating_converging
for (i in names(list_cor)[-1]){
    cor_data2 = cor_data2 + list_cor[[i]]
}
cor_data2 <- cor_data2/14

global_use <- c()

list_data <- list()
for (i in names(list_cor)) {
    for (j in name_all1){
        tmp <- data.frame(Pearson = unlist(list_cor[[i]][j,4:6]),"Dataset"=i,"Metric"=j,"Group"=name_all2)
        list_data[[paste0(i,"_",j)]] = tmp
    }
}
data_plot <- do.call(what = rbind,args = list_data)

write.csv(x = data_plot,file = paste0(aaa,"/","pearson_cor_",bbb,"_level_group.csv"))

pdf(file = paste0(aaa,"/","pearson_cor_",bbb,"_level_group.pdf"),width = 10,height = 6)
ggplot(data = data_plot,aes(x = Metric,y = Pearson))+
    geom_boxplot(mapping = aes(fill = Metric))+
    #stat_summary(fun=mean, geom='point')+
    geom_point(aes(color = Dataset))+
    geom_line(mapping = aes(group = Dataset,color = Dataset))+
    facet_grid(~Group)+
    theme_bw()
dev.off()

################################# cor all
data_use2 = data_use
#data_use2 = data_use[!data_use$Group %in% c("bifurcating_cycle","bifurcating_loop"),]
cor_data <- rcorr(as.matrix(data_use2[,4:9]))$r
cor_data

pvalue_data <- rcorr(as.matrix(data_use2[,4:9]))$P
pvalue_data[is.na(pvalue_data)] <- 0
pdf(file = paste0(aaa,"/","corplot_all_dataset_",bbb,".pdf"),width = 9,height = 7)
corrplot(cor_data, type="upper",p.mat = pvalue_data,col = color_panel,sig.level = 0.5, tl.col = "black",outline = T,diag= F, insig = "blank",mar = c(1,1,1,1))
dev.off()


pdf(file = paste0(aaa,"/","corplot_all_dataset_mixed_",bbb,".pdf"),width = 9,height = 7)
corrplot.mixed(cor_data,upper.col = color_panel,lower.col = color_panel, insig = "blank", tl.col = "black",outline = T,diag= "n",mar = c(1,1,1,1))
dev.off()

