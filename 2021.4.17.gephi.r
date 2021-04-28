setwd("D:\\项目/2021.3.metagenome/Alpha_beta/")

library(Hmisc)
library(psych)

species<-read.csv('species.abundunce.xls',sep=',',row.names=1,encoding="UTF-8",check.names = F)
#genus <- otu[,-1]
#如果是绝对丰度表可以可以依据如下方法转化为相对丰度表（蓝色部分，可选内容）
#species<-species[order(rowSums(species),decreasing = T),]
#species <- head(species,50)
species <- t(species)
species<-species/rowSums(species)
genus<-t(species)
write.table(genus,file='species.txt',sep='\t')
#otu<-t(otu)
#otu<-otu/rowSums(otu)
#otu<-t(otu）
#write.table(otu,file='otu1.txt',sep='\t')
#接下来就是对相对丰度表的一些数据的筛选（蓝色部分，可选内容）
#过滤一些低丰度或低频的类群
genus <- genus[which(rowSums(genus) >= 0.005), ]    #只保留相对丰度总和高于 0.005 的属

#例如只保留在 5 个及以上样本中出现的属
genus1 <- genus
genus1[genus1>0] <- 1
genus <- genus[which(rowSums(genus1) >= 5), ]
#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus), type = 'spearman')
#genus_corr = corr.test(genus,use="pairwise",method="spearman",adjust="fdr",alpha=0.05)
occor_r = genus_corr$r # 取相关性矩阵R值
occor_p = genus_corr$P # 取相关性矩阵p值
# p 值校正，这里使用 BH 法校正 p 值
p <- p.adjust(occor_p, method = 'BH') 
#将相关系数低于0.6和p值大于0.05的值赋值为0 
occor_r[occor_p>0.05|abs(occor_r)<0.6] = 0
diag(occor_r) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.csv(occor_r,file='genus_ccor.csv')#输出csv格式文件

bian <- read.csv('Gephi.bian.csv')
bian[which(bian$Weight > 0),'pn'] <- 'p'
bian[which(bian$Weight < 0),'pn'] <- 'n'
head(bian)
write.csv(bian,file='Gephi.bian.pn.csv',row.names = F)
##！！！！！注意要把第一列删掉