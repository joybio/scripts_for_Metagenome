rm(list=ls())
setwd("D://项目/2021.3.metagenome/diff/")
library(DESeq2)
reads_matrix <- read.csv("total.NumReads",header = T,sep = "\t",row.names=1,check.names=F)
#reads_matrix <- reads_matrix[,-1]
reads_matrix_int <- apply(reads_matrix,2,as.integer)
row.names(reads_matrix_int) <-row.names(reads_matrix)
head(reads_matrix_int)
head(reads_matrix)
#ncol(reads_matrix_int)
#head(reads_matrix_int)
#reads_matrix_int <- head(reads_matrix_int,5000)
#condition <- factor(c(rep("control_c",2),rep("treat_d",2)), levels = c("control_c","treat_d"))
group <- read.csv("group.txt",header = F,sep = "\t")
compare <- read.csv("compare.txt",header = F,sep = "\t")
head(group)
head(compare)
i=1

for (i in 1:nrow(compare)){
  name <- as.character(compare[i,]) #R3.5里必须unlist。例如：name <- as.character(unlist(compare[i,]))
  compare_name <- paste(compare[i,1],compare[i,2],sep="_VS_")
  dir.create(compare_name)
  compare_name_file <- paste0(compare_name,"/",compare_name,".csv")
  group_for_compare <- subset(group,group$V3 %in% name)
  group_for_compare_name <- as.character(group_for_compare$V1)
  compare_name_num <- reads_matrix_int[,group_for_compare_name]
  write.csv(compare_name_num,compare_name_file)
  colData <- data.frame(group_for_compare[,3])
  colnames(colData) <- c("condition")
  row.names(colData) <- group_for_compare$V1
  reads_matrix_sample <- compare_name_num
  head(reads_matrix_sample)
  #colData <- data.frame(row.names=colnames(reads_matrix_sample), condition)
  dds <- DESeqDataSetFromMatrix(reads_matrix_sample, colData, design= ~ condition)
  dds <- DESeq(dds)
  res= results(dds)
  compare_name_pdf <- paste0(compare_name,"/",compare_name,".pdf")
  pdf(compare_name_pdf)
  plotMA(res)
  dev.off()
  res = res[order(res$pvalue),]
  compare_name_result <- paste0(compare_name,"/",compare_name,".result.csv")
  write.csv(res,file=compare_name_result)
  table(res$padj<0.05)
  compare_name_p_result <- paste0(compare_name,"/",compare_name,".diff.result.csv")
  diff_gene_deseq2 <-subset(res, pvalue < 0.05 & abs(log2FoldChange) >= 0.585)
  write.csv(diff_gene_deseq2,compare_name_p_result)
}

