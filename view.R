args<-commandArgs(T)
data<-read.table(args[1],header=F)
library(ggplot2)
library(dplyr)
library(forcats)
path<-args[2]
path1<-gsub("^\\s+|\\s+$", "",path)
p<-data %>%mutate(V1 = fct_reorder(V1, V2)) %>%ggplot( aes(x=V1,y=V2))+geom_boxplot(width=0.4,linetype="dashed")+stat_boxplot(aes(ymin=..lower..,ymax=..upper..))+stat_boxplot(geom = "errorbar",width=0.15,aes(ymax=..ymin..)) +theme_bw()+theme(legend.title = element_text(size=12),legend.text = element_text(colour = 'black',size = 10),,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"),axis.text = element_text(color="black", size=10),axis.title = element_text(color="black", size=12))+xlab("Genotype") + ylab("Gene expression level")
pdf(path1,width = 7,height = 6)
print(p)
dev.off()

