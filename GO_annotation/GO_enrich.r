#加载 clusterProfiler
library(clusterProfiler)
library(ggplot2)
#### GO
#预先获得的 GO 注释文件，基因和 GO 的对应关系
argv <- commandArgs(T)
if (length(argv) != 2 ) {
        cat('\nuseage:\nRscript GO_enrich.r GO_forstat /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/8.Gene_abundunce/test.diff_gene_deseq2.csv\nexample:\nRscript GO_enrich.r /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/6.Function_annotation/GO/GO_forstat /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/8.Gene_abundunce/test.diff_gene_deseq2.csv\n\n')
        q('no')
}

go_anno <- read.csv(argv[1], header = FALSE,sep="\t", stringsAsFactors = FALSE)
names(go_anno) <- c('gene_id', 'ID')

#添加 GO 注释详情
#所有 GO 的注释描述均可在 GO 官网自定义获取：http://geneontology.org/
go_class_all <- read.delim('/media/ruizhi/database/GO/godb.csv', header = T,sep=",", stringsAsFactors = FALSE)
go_class <- go_class_all[,c("GOID","TERM","ONTOLOGY")]
names(go_class) <- c('ID', 'Description', 'Ontology')
go_anno <- merge(go_anno, go_class, by = 'ID', all.x = TRUE)

#目标基因列表

gene_select <- read.delim(argv[2], header = FALSE,sep=",",stringsAsFactors = FALSE)[,1]

#GO 富集分析
#默认以所有注释到 GO 的基因为背景集，也可通过 universe 参数输入背景集
#默认以 p<0.05 为标准，Benjamini 方法校正 p 值，q 值阈值 0.2
#默认输出 top500 富集结果
go_rich <- enricher(gene = gene_select,
    TERM2GENE = go_anno[c('ID', 'gene_id')], 
    TERM2NAME = go_anno[c('ID', 'Description')], 
	pvalueCutoff = 0.05, 
    pAdjustMethod = 'BH', 
    qvalueCutoff = 0.2, 
    maxGSSize = 500)

#输出默认结果，即根据上述 p 值等阈值筛选后的
write.table(go_rich, 'go_tmp.txt', sep = '\t', row.names = FALSE, quote = FALSE)
tmp <- read.delim('go_tmp.txt')
tmp <- merge(tmp, go_class[c('ID', 'Ontology')], by = 'ID')
tmp <- tmp[c(10, 1:9)]
tmp <- tmp[order(tmp$pvalue), ]
write.table(tmp, 'go_rich.significant.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#同样地，如果想输出所有富集结果（不考虑 p 值阈值等），将 p、q 等值设置为 1 即可
#或者直接在 enrichResult 类对象中直接提取需要的结果
#names(attributes(go_rich))  #查看对象“go_rich”包含元素
#result <- go_rich@result    #其中的 result，即为所有基因的富集结果

#result <- merge(result, go_class[c('ID', 'Ontology')], by = 'ID')
#result <- result[c(10, 1:9)]
#result <- result[order(result$pvalue), ]
#write.table(result, 'go_rich.all.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#一些常见作图方法
#上文提到的 clusterProfiler 包自带的柱形图、气泡图、网络图等默认方法略
#这些参考上文 KEGG
pdf("go_all.dot.pdf",width=12,height=8)
dotplot(go_rich)
dev.off()

pdf("go_all.bar.pdf",width=12,height=8)
barplot(go_rich)
dev.off()


ALL.BP.df <- go_rich@result[as.vector(subset(tmp, Ontology == 'BP')$ID), ]
ALL.CC.df <- go_rich@result[as.vector(subset(tmp, Ontology == 'CC')$ID), ]
ALL.MF.df <- go_rich@result[as.vector(subset(tmp, Ontology == 'MF')$ID), ]
ALL.BP.df$Class <- rep("biological_process")
ALL.CC.df$Class <- rep("cellular_component")
ALL.MF.df$Class <- rep("molecular_function")
ALL.BP.df <- ALL.BP.df[order(-ALL.BP.df$Count),]
ALL.CC.df <- ALL.CC.df[order(-ALL.CC.df$Count),]
ALL.MF.df <- ALL.MF.df[order(-ALL.MF.df$Count),]
ALL.BP.df <- head(ALL.BP.df,n=10)
ALL.CC.df <- head(ALL.CC.df,n=10)
ALL.MF.df <- head(ALL.MF.df,n=10)
GO <- rbind(ALL.BP.df,ALL.CC.df,ALL.MF.df)
CPCOLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
#CPCOLS<-c("#CC6666", "#9999CC", "#66CC99")
#CPCOLS<-c("#999999", "#E69F00", "#56B4E9")
dorder = factor(as.character(GO$Description),levels = rev(as.character(GO$Description)))
p <- ggplot(GO,aes(x=Description,y=Count,fill=Class)) 
p + geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder)) +  
  coord_flip() + 
  scale_y_log10(breaks = c(1,10,100,1000)) +
  scale_fill_manual(values = CPCOLS)
#scale_fill_discrete(name="Ontology") + 
theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
  xlab("Term") 
ggsave("GO_annotation_all.pdf",width = 10,height = 10)
dev.off()

#此外，对于 GO 的包含关系，即 GO 的有向无环图，在这种纯自定义富集分析的情况下，仍然可以使用 topGO 展示
#此时需要在作图前手动修改 enrichResult 类对象中的元素，主要为 result 统计表和 ontology 的 GO 分类，便于 topGO 识别
#但需注意的是一次只能展示一类（BP、CC 或 MF）GO，如下以展示富集的“biological process”（即 BP）为例

#BiocManager::install('topGO')
#BiocManager::install('Rgraphviz')


pdf("go_topGO_BP.pdf",width=12,height=8)
library(topGO)

go_rich_BP <- go_rich
go_rich_BP@result <- go_rich_BP@result[as.vector(subset(tmp, Ontology == 'BP')$ID), ]
go_rich_BP@ontology <- 'BP'
plotGOgraph(go_rich_BP)
dev.off()

pdf("go_topGO_CC.pdf",width=12,height=8)
go_rich_CC <- go_rich
go_rich_CC@result <- go_rich_CC@result[as.vector(subset(tmp, Ontology == 'CC')$ID), ]
go_rich_CC@ontology <- 'CC'
plotGOgraph(go_rich_CC)
dev.off()

pdf("go_topGO_MF.pdf",width=12,height=8)
go_rich_MF <- go_rich
go_rich_MF@result <- go_rich_MF@result[as.vector(subset(tmp, Ontology == 'MF')$ID), ]
go_rich_MF@ontology <- 'MF'
plotGOgraph(go_rich_MF)
dev.off()
