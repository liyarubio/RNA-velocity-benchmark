library(redeemR)
library(ggtree)
library(tidytree)
library(data.table)

setwd("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM")

####-------- data process by redeemR -----------
# Dir="/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone/Youn2.BMMC.Consensus.final"
Dir="figshare/clone/Youn2.HPC.Consensus.final"
VariantsGTSummary<-redeemR.read(path=Dir,thr='S',Processed=F)

# t = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone/Youn2.BMMC.Consensus.final/VariantsGTSummary.RDS")
# t2 = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/GSE219015_RAW/GSM8133910_Extended_Young_BMMC.mito.depth.RDS")

head(VariantsGTSummary,n=3)

Example_redeemR<-Create_redeemR(VariantsGTSummary) # This is using default parameter
# The full statement is Example_redeemR<-Create_redeemR(VariantsGTSummary,qualifiedCellCut=10,VAFcut=1,Cellcut=2,maxctscut=2)

head(Example_redeemR@GTsummary.filtered,n=2)
head(Example_redeemR@V.fitered,n=2)

plot_depth(Example_redeemR)

MutationProfile.bulk(Example_redeemR@UniqueV)

Example_redeemR<-Make_matrix(Example_redeemR,onlyhetero=T)
Example_redeemR@Cts.Mtx[1:3,1:3]

Example_redeemR<-suppressMessages(SeuratLSIClustering(Example_redeemR,lsidim=2:50))

data(CellPCT)
V.weight<-data.frame(weight=1-CellPCT$muRate)
V.weight$Variants<-paste("Variants",gsub("_","",CellPCT$Variant),sep="")

Example_redeemR<-AddDist(Example_redeemR,weightDF=V.weight, LSIdist=T) # if LSIdist=T, SeuratLSIClustering needs to be run before. 

Example_redeemR.mnn.igraph<-FromDist2Graph(Example_redeemR@DistObjects@w_jaccard)  # The input can be other distances as well

plot(Example_redeemR.mnn.igraph,vertex.size = 5,vertex.label = NA, edge.label = NA) 

Example_redeemR.nn<-MakeNN(Example_redeemR@DistObjects@w_jaccard)
Example_redeemR.knn.matrix<-NN2M(Example_redeemR.nn)
Example_redeemR.knn.igraph<-igraph::graph_from_adjacency_matrix(Example_redeemR.knn.matrix,diag = F,mode = "undirected")
Example_redeemR<-Make_tree(Example_redeemR,d = "w_jaccard",algorithm = "nj")

treeplot<-ggtree(Example_redeemR@TREE@treedata,layout="circular", branch.length='none') 
print(treeplot)

# saveRDS(Example_redeemR,"/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone/Youn2.BMMC.Consensus.final/young2_BMMC.rds")
saveRDS(Example_redeemR,"/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone/young2_HPC.rds")

# y2_BMMC = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone/Youn2.BMMC.Consensus.final/young2_BMMC.rds")
y2_BMMC = readRDS("figshare/clone/young2_HPC.rds")

treeplot<-ggtree(y2_BMMC@TREE@treedata,layout="circular", branch.length='none') 
print(treeplot)

####-----redeem tree-----------
# update the `@CellMeta` by adding VN(variant number) of mitochondrial mutations and TreePos(cell position on the tree). 

library(Matrix)
library(dplyr)
CellTreePos<-subset(treeplot$data,isTip==TRUE) %>% .[order(.$y),] %>% .[,c("label","node","y")] %>% rename(Cell=label,node=node,TreePos=y) # This extract the position of cells on the tree. This order information will be stored in the factor of cells
VN.summary<-data.frame(VN=rowSums(y2_BMMC@Cts.Mtx.bi)) %>% tibble::rownames_to_column("Cell") %>% merge(.,CellTreePos) # This summarize the mtDNA mutations per cell
if(!"TreePos" %in% colnames(y2_BMMC@CellMeta)){
  y2_BMMC@CellMeta<-merge(y2_BMMC@CellMeta,VN.summary)
}    

head(y2_BMMC@CellMeta,n=5)

path = "figshare/clone/Youn2.HPC.Consensus.final/QualifiedTotalCts"
q <-read.table(path)

y2_BMMC<-Add_DepthMatrix(y2_BMMC,q)
y2_BMMC<-Add_AssignVariant(y2_BMMC,n.cores=8)

Allnodes=MakeAllNodes(y2_BMMC,prob.cut=0.1)
y2_BMMC<-Add_tree_cut(y2_BMMC,MinCell = 2,N = 1)
head(y2_BMMC@CellMeta)

saveRDS(y2_BMMC,"figshare/clone/young2_HPC_tree_cut.rds")

redeemR = y2_BMMC
phy=redeemR@TREE@phylo   

redeemR = readRDS("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone/young2_HPC_tree_cut.rds")

redeemR@AssignedVariant$Variant.assign.report %>% head(n=50)  # For example

head(redeemR@CellMeta)

m = redeemR@CellMeta

t = table(m$Clade_merge) %>% as.data.frame()
tt = table(m$Clone_merge) %>% as.data.frame()

m$Cell = Translate_simple_ATAC2RNA(m$Cell)

setwd("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone")
fwrite(m,"young2_hspc_clone.csv")

####---------- tree -----------
tree = redeemR@TREE
p = tree@phylo
t = tree@treedata

x <- as_tibble(t)
x

x$label = Translate_simple_ATAC2RNA(x$label)

p = as.phylo(x)

setwd("/home/liyaru/DATA/8_transition/1_DATA/6_Nature_ReDeeM/figshare/clone")
library(tidytree)
library(treeio)
write.tree(p,file = "young2_HPC_tree.nwk")



# same as
# xx <- as_tibble(p)
# xx
child(x,1)
parent(x,10)

child(x,"AAACAAGCAGGTTCTT")
parent(x,"AAACAAGCAGGTTCTT")

xx = x[!is.na(x$label),]
tt = table(xx$parent) %>% as.data.frame()


####----------- Tree test----------------
library(ape)
library(tidytree)
set.seed(2024)

tree1 <- rtree(4)
tree1

x1 <- as_tibble(tree1)
x1
as.phylo(x1)

d <- tibble(label = paste0('t', 1:4),
            trait = rnorm(4))

y <- full_join(x, d, by = 'label')
y

as.treedata(y)

y %>% as.treedata %>% as_tibble

child(x1, 5)
parent(x1, 5)











