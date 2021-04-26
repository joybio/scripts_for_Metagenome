rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
argv <- commandArgs(T)
if (length(argv) != 1 ) {
        cat('\nuseage:\nRscript barplot.r stast\nexample:\nRscript barplot.r GO.secondary.annotation.csv\n\n')
        q('no')
}
data <- read.csv(argv[1],header = T,sep=",",stringsAsFactors = F)
#data <- read.csv(file = "GO.secondary.annotation.csv",header = T,sep=",",stringsAsFactors = F)
#data <- data[order(data[,"ONTOLOGY"],-data[,"GENE_NUM"]),]
data <- dplyr::arrange(data,data[,"ONTOLOGY"],-data[,"GENE_NUM"])
data_BP <- subset(data,data$ONTOLOGY == "BP")
if(nrow(data_BP) >=10){
  data_BP <- head(data_BP,n=10)
}
data_CC <- subset(data,data$ONTOLOGY == "CC")
if(nrow(data_CC) >=10){
  data_CC <- head(data_CC,n=10)
}
data_MF <- subset(data,data$ONTOLOGY == "MF")
if(nrow(data_MF) >=10){
  data_MF <- head(data_MF,n=10)
}
data <- rbind(data_BP,data_CC,data_MF)
CPCOLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
data2 <- tidyr::unite(data,"Description",GO_ID,TERM,sep=":") 
dorder = factor(as.character(data2$Description),levels = rev(as.character(data2$Description))) #
p <- ggplot(data2,aes(x=Description,y=GENE_NUM,fill=ONTOLOGY))
p + geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder)) +  
  coord_flip() + 
  scale_y_log10(breaks = c(1,10,100,1000)) + 
  scale_fill_manual(values = CPCOLS)
#scale_fill_discrete(name="Ontology") +
theme(panel.background = element_rect(fill = "transparent",colour = NA)) +
  xlab("Term")
ggsave("GO_Secondary_annotation.pdf",width = 10,height = 10)
dev.off()



