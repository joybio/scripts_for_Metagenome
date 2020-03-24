library("ggplot2")
data <- read.csv("../stast",header = T,sep="\t")
p <- ggplot(data,aes(x=mark,y=number,fill=description)) + geom_bar(stat='identity',width=0.6) + scale_y_continuous(expand = c(0,0),limits=c(0,400)) + guides(fill = guide_legend( ncol = 1, byrow = TRUE))
#,position=position_dodge(0.7),width=0.5
#guides(fill = guide_legend( ncol = 1, byrow = TRUE)) one coloum
p + xlab("Function class") + ylab("Number of matched genes") + theme(legend.text = element_text(color = "black", size = 8, family="Times"))
ggsave("COG function classification.pdf",width=10,height=10)

