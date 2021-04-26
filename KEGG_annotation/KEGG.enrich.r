#Bioconductor 安装 clusterProfiler
#BiocManager::install('clusterProfiler')
#加载 clusterProfiler
library(clusterProfiler)
library(ggplot2)
rm(list = ls())
#### KEGG
#“kegg_anno.txt”为预先获得的 KEGG 注释文件，包含基因和 KEGG pathway 的对应关系
#我用蛋白序列做的 kegg 注释，后续又对应的蛋白-基因关系，不用担心同一途径中基因名称出现重复的问题，enricher() 函数可自动去重
argv <- commandArgs(T)
if (length(argv) != 2 ) {
        cat('\nuseage:\nRscript KEGG_enrich.r kegg.ko.list /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/8.Gene_abundunce/test.diff_gene_deseq2.csv\nexample:\nRscript GO_enrich.r /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/6.Function_annotation/eggNOG/kegg.ko.list /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/8.Gene_abundunce/test.diff_gene_deseq2.csv\n\n')
        q('no')
}

kegg_anno <- read.delim(argv[1],header = T, colClasses = 'character')
#colnames(kegg_anno) <- c("gene_id","pathway_id","pathway_description")
head(kegg_anno)

#目标基因列表
gene_select <- read.delim(argv[2], header = FALSE,sep=",",stringsAsFactors = FALSE)[,1]
head(gene_select)

#KEGG 富集分析
#默认以所有注释到 KEGG 的基因为背景集，也可通过 universe 参数指定其中的一个子集作为背景集
#默认以 p<0.05 为标准，Benjamini 方法校正 p 值，q 值阈值 0.2
#默认输出 top500 富集结果
kegg_rich <- enricher(gene = gene_select,
                      TERM2GENE = kegg_anno[c('pathway_id', 'gene_id')], 
                      TERM2NAME = kegg_anno[c('pathway_id', 'pathway_description')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH', 
                      qvalueCutoff = 0.2, 
                      maxGSSize = 500)

#输出默认结果，即根据上述 p 值等阈值筛选后的
write.table(kegg_rich, 'kegg_rich.significant.txt', sep = '\t', row.names = FALSE, quote = FALSE)

kegg_rich.df <- kegg_rich@result
write.table(kegg_rich.df, 'kegg_rich.txt', sep = '\t', row.names = FALSE, quote = FALSE)

ALL.KEGG.df.fifteen <- head(kegg_rich.df,n=10)
ALL.KEGG.df.fifteen$padjust = -log10(ALL.KEGG.df.fifteen$pvalue)
B=as.data.frame(ALL.KEGG.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="KEGG.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #璁剧疆x杞村埢搴︽爣绛剧殑瀛椾綋鏄剧ず鍊炬枩瑙掑害涓?15搴︼紝骞跺悜涓嬭皟鏁?1(hjust = 1)锛屽瓧浣撶皣涓篢imes澶у皬涓?20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #璁剧疆y杞村埢搴︽爣绛剧殑瀛椾綋绨囷紝瀛椾綋澶у皬锛屽瓧浣撴牱寮忎负plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #璁剧疆y杞存爣棰樼殑瀛椾綋灞炴€?
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #鍘婚櫎榛樿濉厖鐨勭伆鑹诧紝骞跺皢x=0杞村拰y=0杞村姞绮楁樉绀?(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

#clusterProfiler 包里的一些默认作图方法，例如
pdf("kegg.dot.pdf",width=10,height=10)
dotplot(kegg_rich)  #富集气泡图
dev.off()

pdf("kegg.cnet.pdf",width=10,height=10)
cnetplot(kegg_rich) #网络图展示富集功能和基因的包含关系
dev.off()
pdf("kegg.ema.pdf",width=10,height=10)
emapplot(kegg_rich) #网络图展示各富集功能之间共有基因关系
dev.off()
pdf("kegg.heat.pdf",width=10,height=10)
heatplot(kegg_rich) #热图展示富集功能和基因的包含关系
dev.off()
#如果想输出所有富集结果（不考虑 p 值阈值等）
#可以在上述 enricher() 函数中，将 p、q 等阈值设置为 1

pdf("kegg.pdf")
dotplot(kegg_rich)  #富集气泡图
cnetplot(kegg_rich) #网络图展示富集功能和基因的包含关系
emapplot(kegg_rich) #网络图展示各富集功能之间共有基因关系
heatplot(kegg_rich) #热图展示富集功能和基因的包含关系
dev.off()

#或者在 enrichResult 类对象中直接提取需要的结果
#names(attributes(kegg_rich))  #查看对象“kegg_rich”包含元素
#result <- kegg_rich@result    #其中的 result，即为所有基因的富集结果

#write.table(result, 'kegg_rich.all.txt', sep = '\t', row.names = FALSE, quote = FALSE)



