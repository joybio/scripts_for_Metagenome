################################################################################
##############################物种分类及主成分##################################
################################################################################
#载入相应的包
#如果未安装，请按下行代码安装
#install.package(c("ggplot2","reshape2")
setwd("D:\\2021.3.metagenome/Alpha_beta/")

library(ggplot2)
library(reshape2)
#设置工作路径
#载入数据，本数据是鱼类物种组成的绝对值数据，在excel中保存为csv格式．

taxa<-read.table("species.merged_metaphlan2.txt",row.names = 1,header = T)
taxa <- taxa[,-1]
###############################################################################
#物种组成柱形图

#将文件中的NA值改为０
taxa[is.na(taxa)] <- 0
#将数据中的种类根据数量的多少进行排序
taxa<-taxa[order(rowSums(taxa),decreasing = T),]
#Ｎ值代表选择数量排前10的物种，将剩下的物种合并成其他
N<-10
taxa_list<-rownames(taxa)[1:N]
new_x<-rbind(taxa[row.names(taxa) %in% taxa_list,],
             others=rowSums(taxa[!rownames(taxa) %in% taxa_list,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
#作图
dir.create("abundunce")
pdf("abundunce/metaphlan.abundunce.bar.plot.pdf",width = 10, height = 8)

ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
dev.off()

###############################################################################
#聚类分析
#清理旧数据
#设置工作路径
#载入工作包
library(ggplot2)
library(ggdendro)
library(vegan)
#读取数据
taxa<-read.table("species.merged_metaphlan2.txt",row.names = 1,header = T)
taxa[is.na(taxa)]<-0
taxa <- taxa[,-1]
#采用Bray Curtis方法，如需要更换其他方法，可在method参数中调整
beta_bray<-vegdist(t(taxa),method="bray")
#建树
hc<-hclust(beta_bray)
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type = "rectangle")
#绘图
pdf("normal.tree.pdf",width = 30, height = 16)
ggplot(dend_data$segments) + 
  theme_dendro()+
  scale_x_discrete(expand = c(0,1))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            size = 5,check_overlap = T,angle=45,vjust = 3,
            nudge_y = -0.02)
dev.off()
#α多样性分析

rm(list = ls())
###############################################################################
#PCA分析（主成分分析）
###############################################################################

library(ggbiplot)

#读取数据
taxa<-read.csv("species.abundunce.xls",header = T,row.names = 1)
meta<-read.csv("map.csv",header = T,sep = ",",stringsAsFactors = F,row.names = 1)
taxa<-t(taxa)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]

#计算PCA值
pca<-prcomp(taxa,scale. = T)
#作图
dir.create("Beta")

pdf("Beta/PCA.pdf")
ggbiplot(pca, obs.scale = 2, var.scale = 1,var.axes = F,
         groups = meta$group, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()
###############################################################################

library(ade4)
library(ggplot2)
library(RColorBrewer)
data(deug)
#https://www.jianshu.com/p/a40715fa9b04
#PCA分析
taxa<-read.csv("species.abundunce.xls",header = T,row.names = 1)
meta<-read.csv("species.annotation.xls",header = T,sep = "\t",row.names = 1,stringsAsFactors = F)
group <- read.csv("map.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = F)
#PCA分析
taxa<-t(taxa)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]

pca<- dudi.pca(taxa, scal = T, scan = FALSE)

#坐标轴解释量（前两轴）
pca_eig <- (pca$eig)[1:2] / sum(pca$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pca$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCA1', 'PCA2')

#以group为分组
sample_site$level<-factor(group$group)

library(ggplot2)
pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=level)) +
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
  scale_color_manual(values = brewer.pal(6,"Set2")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCA1: ', round(100 * pca_eig[1], 2), '%'), y = paste('PCA2: ', round(100 * pca_eig[2], 2), '%')) 
ggsave("Beta/PCA-2.pdf")

################################################################################
##############################PCoA##################################
################################################################################

library(ade4)
library(ggplot2)
library(RColorBrewer)
library(vegan)
taxa<-read.csv("species.abundunce.xls",header = T,row.names = 1)
meta<-read.csv("species.annotation.xls",header = T,sep = "\t",row.names = 1,stringsAsFactors = F)
group <- read.csv("map.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = F)
#PCA分析
taxa<-t(taxa)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]
tab.dist<-vegdist(taxa,method='euclidean')
pcoa<- dudi.pco(tab.dist, scan = T,nf=3)
####此处需要根据柱形图选择nf数目
#坐标轴解释量（前两轴）
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pcoa$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

#以group作为分组
sample_site$level<-factor(group$group)

library(ggplot2)
pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2,color=level)) +
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
  scale_color_manual(values = brewer.pal(6,"Set2")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) 
ggsave("Beta/PCoA.pdf")

pcoa_plot


################################################################################
##############################Alpha##################################
################################################################################
#http://blog.sciencenet.cn/home.php?mod=space&uid=3406804&do=blog&id=1184055
#读取 species 丰度表
rm(list=ls())
species <- read.csv('species.abundunce.xls', header = T,row.names = 1,check.names = F)
#species <- species[,-1]
# <- apply(species, 1, as.integer)
#species <- t(species)
species <- apply(species,2,function(x)x*50000/(sum(x)))
species<-species[order(rowSums(species),decreasing = T),]
species <- head(species,50)
species <- apply(species,2,round)
species <- as.data.frame(t(species))
library(vegan)    #用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(picante)    #用于计算 PD_whole_tree，若不计算它就无需加载。事实上，picante 包加载时默认同时加载 vegan

##定义函数
#计算多种 Alpha 多样性指数，结果返回至向量
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}
dir.create("Alpha")
#统计 species 丰度表中各样本的 Shannon 指数，对数底数使用 e
shannon_index <- alpha_index(species, method = 'shannon', base = exp(1))
write.csv(shannon_index,"Alpha/shannon_index.csv")
richness_index <- alpha_index(species, method = 'richness', base = exp(1))
write.csv(richness_index,"Alpha/richness_index.csv")
chao1_index <- alpha_index(species, method = 'chao1', base = exp(1))
ace_index <- alpha_index(species, method = 'ace', base = exp(1))
simpson_index <- alpha_index(species, method = 'simpson', base = exp(1))
write.csv(simpson_index,"Alpha/simpson_index.csv")
pielou_index <- alpha_index(species, method = 'pielou', base = exp(1))
gc_index <- alpha_index(species, method = 'gc', base = exp(1))



#以 1000 条序列为抽样步长，依次对 species 表稀释抽样，直到最大序列深度；并统计各抽样梯度下的 species 丰度表中各样本的 Shannon 指数，对数底数使用 e
#shannon_curves <- alpha_curves(species, step = 1000, method = 'shannon', base = exp(1))

#接下来以最简单的Richness指数（物种丰富度指数）为例，展示Alpha多样性稀释曲线的绘制方法。
#多计算几次以获取均值 ± 标准差，然后再展示出也是一个不错的选择
#重复抽样 5 次
#############################################################################
plot_richness <- data.frame()
for (n in 1:5) {
  richness_curves <- alpha_curves(species, step = 2000, method = 'richness')
  
  for (i in names(richness_curves)) {
    richness_curves_i <- (richness_curves[[i]])
    richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
    plot_richness <- rbind(plot_richness, richness_curves_i)
  }
}
library(doBy)
library(ggplot2)
#计算均值 ± 标准差（doBy 包中的 summaryBy() 函数）
plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA
#ggplot2 作图

ggplot(plot_richness_stat, aes(rare, alpha.mean, color = sample)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = alpha.mean - alpha.sd, ymax = alpha.mean + alpha.sd), width = 500) +
  labs(x = 'number of normolised reads', y = 'Richness', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(species)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 1000000, 200000), labels = as.character(seq(0, 1000000,200000)))
ggsave("Alpha/richness.pdf",width = 20,height = 20)
dev.off()
#############################################################################



#Shannon指数稀释曲线等
##若简单的“geom_line()”样式波动幅度过大，不平滑等，可以尝试拟合曲线的样式
#获得作图数据。前面多生成一个点，使得 Shannon 拟合曲线更加平滑（你把 shannon_curves1 注释掉就知道我说的啥了）
shannon_curves1 <- alpha_curves(species, step = 200, rare = 200, method = 'shannon')
shannon_curves2 <- alpha_curves(species, step = 2000, method = 'shannon')
shannon_curves <- c(shannon_curves1, shannon_curves2)

plot_shannon <- data.frame()
for (i in 1:length(shannon_curves)) {
  shannon_curves_i <- shannon_curves[[i]]
  shannon_curves_i <- data.frame(rare = names(shannon_curves_i), alpha = shannon_curves_i, sample = names(shannon_curves)[i], stringsAsFactors = FALSE)
  plot_shannon <- rbind(plot_shannon, shannon_curves_i)
}

rownames(plot_shannon) <- NULL
plot_shannon$rare <- as.numeric(plot_shannon$rare)
plot_shannon$alpha <- as.numeric(plot_shannon$alpha)
plot_shannon <- plot_shannon[order(plot_shannon$sample, plot_shannon$rare), ]

#ggplot2 作图（使用到 ggalt 包的 geom_xspline() 绘制平滑拟合线）
library(ggalt)    #若未加载时先加载
ggplot(plot_shannon, aes(rare, alpha, color = sample)) +
  geom_xspline() +
  labs(x = 'Number of normolised reads', y = 'Shannon', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(species)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 100000, 20000), labels = as.character(seq(0, 10000, 2000)))
ggsave("Alpha/Shannon Index.pdf",width = 20,height = 20)
dev.off()

#Rank-abundance曲线在R中的绘制方法
#虽然自己排序也很简单，但对于我这样的懒人来讲，还是导个包统计省事点……
library(BiodiversityR)

#统计（BiodiversityR 包 rankabundance() 实现 species 排序）
species_relative <- species/rowSums(species)
rank_dat <- data.frame()
for (i in rownames(species_relative)) {
  rank_dat_i <- data.frame(rankabundance(subset(species_relative, rownames(species_relative) == i), digits = 6))[1:2]
  rank_dat_i$sample <- i
  rank_dat <- rbind(rank_dat, rank_dat_i)
}
rank_dat <- subset(rank_dat, rank_dat$abundance != 0)

#ggplot2 作图
dev.off()
ggplot(rank_dat, aes(rank, log(abundance, 10), color = sample)) +
  geom_line() +
  labs(x = 'speciess rank', y = 'Relative adundance (%)', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  scale_y_continuous(breaks = 0:-5, labels = c('100', '10', '1', '0.1', '0.01', '0.001'), limits = c(-5, 0))
ggsave("Alpha/Rank-abundance.pdf",width = 20,height = 20)




################################################################################
##########################beta##########################################
################################################################################
setwd("D:\\2021.3.metagenome/Alpha_beta/")

library("phyloseq")
library("ggplot2")
library("plyr")

taxa<-read.csv("species.abundunce.xls",header = T,row.names = 1,check.names = F)
meta<-read.csv("species.annotation.xls",header = T,sep = "\t",row.names = 1,stringsAsFactors = F)

map <- read.csv("map.csv",header = T,sep = ",",stringsAsFactors = F)

#将数据制成phyloseq的格式
#otu_table的输入文件必须是data.frame
TAXA<-otu_table(taxa,taxa_are_rows = T)
TAXA[1:5,1:5]
#tax_table的输入文件必须是matrix
meta <- data.frame(meta)
taxmeta <- as.matrix(meta)
dim(taxmeta)
rownames(taxmeta) <- rownames(taxa)
colnames(taxmeta) <- colnames(meta)
SAD = tax_table(taxmeta)
SAD[1:5,1:5]
################################################################################
?tax_table
################################################################################
#phy_tree <- data.frame(data=NULL)
#phy_tree <- phy_tree(phy_tree)
#sample_data的输入文件必须是data.frame
map_matrix <- as.data.frame(map)
row.names(map_matrix) <- map_matrix[,1]
sample_data<-sample_data(map_matrix)
intersect(row.names(sample_data),colnames(taxa))
################################################################################
?sample_data
data(soilrep)
head(sample_data(soilrep))
################################################################################
phy<-phyloseq(TAXA,SAD,sample_data)
pdf("abundunce/Domain.abundunce.pdf")
plot_bar(phy, fill = "Domain")
dev.off()
pdf("abundunce/Phylum.abundunce.pdf")
plot_bar(phy, fill = "Phylum")
dev.off()
pdf("abundunce/Order.abundunce.pdf",width = 30,height = 10)
plot_bar(phy, fill = "Order")
dev.off()
pdf("abundunce/Family.abundunce.pdf",width = 60,height = 10)
plot_bar(phy, fill = "Family")
dev.off()
pdf("abundunce/Genus.abundunce.pdf",width = 100,height = 10)
plot_bar(phy, fill = "Genus")
dev.off()
pdf("abundunce/Species.abundunce.pdf",width = 400,height = 10)
plot_bar(phy, fill = "Species")
dev.off()
?plot_bar


pdf("Alpha/Alpha Diversity Measure.pdf")
plot_richness(phy,x=nsamples(phy),color="group",measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
dev.off()
?plot_richness
?phyloseq
?nsamples
##################################################
#library("ape")
#random_tree = rtree(ntaxa(phy), rooted=TRUE, tip.label=taxa_names(phy))
#plot(random_tree)
##################################################



dist_methods <- unlist(distanceMethodList)
print(dist_methods)
# These require tree
dist_methods[(1:3)]
# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
# This is the user-defined method:
dist_methods["designdist"]
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods = dist_methods[-which(dist_methods=="raup")]

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
#######################################################################
#test
##iDist <- distance(phy, method="bray")
#iMDS  <- ordinate(phy, "MDS", distance=iDist)
#p <- NULL
#p <- plot_ordination(phy, iMDS, color="group", shape="group")
#p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))


#######################################################################
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(phy, method=i)
  # Calculate ordination
  iMDS  <- ordinate(phy, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(phy, iMDS, color="group", shape="group")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
#Shade according to sequencing technology
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=group, shape=group))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for metagenome dataset")
ggsave("Beta/MDS on various distance metrics for metagenome dataset.pdf",width = 20,height = 20)
dev.off()
#Jensen-Shannon Divergence
pdf("Beta/jensen-shannon diversity.pdf")
print(plist[["jsd"]])
dev.off()
##############################################
#Jaccard
pdf("Beta/jaccard diversity.pdf")
print(plist[["jaccard"]])
dev.off()

##############################################
#Bray-Curtis
pdf("Beta/Bray-Curtis diversity.pdf")
print(plist[["bray"]])
dev.off()
#Gower
pdf("Beta/Grower diversity.pdf")
print(plist[["gower"]])
dev.off()
#morisita
pdf("Beta/morisita diversity.pdf")
print(plist[["morisita"]])
dev.off()
##############################################
# betadiver1 w
pdf("Beta/betadiver1 diversity.pdf")
print(plist[["w"]])
dev.off()
#cao
pdf("Beta/cao diversity.pdf")
print(plist[["cao"]])
dev.off()
#chao
pdf("Beta/chao diversity.pdf")
print(plist[["chao"]])
dev.off()


