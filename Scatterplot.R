
#First condition NCIL-1 VS OEIL-1
group1 <- read.table("./Scatterplot/分组信息1.txt",sep="\t",row.names=1,header=TRUE,as.is = TRUE)
group2 <- group1[,'group']
tTag1 <- as.data.frame(res1_total)
g1 <- unique(group2)[1]
g2 <- unique(group2)[2]
diff1 <- tTag1[((abs(tTag1$log2FoldChange) >= 0.58)),]

g1_exp = NCILVSOEIL[rownames(tTag1),rownames(group1)[which(group1$group==g1)]]
g2_exp = NCILVSOEIL[rownames(tTag1),rownames(group1)[which(group1$group==g2)]]

g1_mean = apply(g1_exp,1,mean)
g2_mean = apply(g2_exp,1,mean)

type1=rep('No',length(g1_mean))

type1[which(tTag1$log2FoldChange > 0.58)] = "Up"

type1[which(tTag1$log2FoldChange < -0.58)] = "Down"

datam1 = data.frame(g1_mean,g2_mean,logFC=tTag1$log2FoldChange,pvalue=tTag1$pvalue,type1,stringsAsFactors=FALSE)

ggplot(datam1,aes(log2(g1_mean),log2(g2_mean),colour=type1))+
  
  geom_point(stat="identity",size=1)+theme(legend.title=element_blank())+scale_color_manual(values =c("Down"='DarkOrange',"No"='grey',"Up"='MediumAquamarine'))+
  
  labs(x=paste(g1,' Log2(TPM)'),y=paste(g2,' Log2(TPM)'),title=paste(g2,' VS ',g1,sep=""))+
  
  coord_cartesian(ylim=c(-10,10),xlim=c(-10,10))+geom_segment(aes(x = -10, y = -10, xend = 10, yend = 10),size=1,colour="#999999",linetype="dotted")+theme(plot.title = element_text(hjust = 0.5),title=element_text(face="bold",size=15,colour="black"),axis.title=element_text(face="bold",size=13,colour="black"),axis.text.x=element_text(face="bold",size=12,colour="black"),axis.text.y=element_text(face="bold",size=12,colour="black"),legend.text=element_text(face="bold",size=13,colour="black"))


#Second condition NC VS NCIL-1
group3 <- read.table("./Scatterplot/分组信息2.txt",sep="\t",row.names=1,header=TRUE,as.is = TRUE)
group4 <- group3[,'group']
tTag2 <- as.data.frame(res2)
g3 <- unique(group4)[1]
g4 <- unique(group4)[2]

g3_exp = NCVSNCIL[rownames(tTag2),rownames(group3)[which(group3$group==g3)]]
g4_exp = NCVSNCIL[rownames(tTag2),rownames(group3)[which(group3$group==g4)]]

g3_mean = apply(g3_exp,1,mean)
g4_mean = apply(g4_exp,1,mean)

type2=rep('No',length(g3_mean))

type2[which(tTag2$log2FoldChange > 0.58 & tTag2$pvalue<0.05)] = "Up"

type2[which(tTag2$log2FoldChange < -0.58& tTag2$pvalue<0.05)] = "Down"

datam2 = data.frame(g3_mean,g4_mean,logFC=tTag2$log2FoldChange,pvalue=tTag2$pvalue,type2,stringsAsFactors=FALSE)

ggplot(datam2,aes(log2(g3_mean),log2(g4_mean),colour=type2))+
  
  geom_point(stat="identity",size=1)+theme(legend.title=element_blank())+scale_color_manual(values =c("Down"='MediumAquamarine',"No"='grey',"Up"='DarkOrange'))+
  
  labs(x=paste(g3,' Log2(TPM)'),y=paste(g4,' Log2(TPM)'),title=paste(g4,' VS ',g3,sep=""))+
  
  coord_cartesian(ylim=c(-10,10),xlim=c(-10,10))+geom_segment(aes(x = -10, y = -10, xend = 10, yend = 10),size=1,colour="#999999",linetype="dotted")+theme(plot.title = element_text(hjust = 0.5),title=element_text(face="bold",size=15,colour="black"),axis.title=element_text(face="bold",size=13,colour="black"),axis.text.x=element_text(face="bold",size=12,colour="black"),axis.text.y=element_text(face="bold",size=12,colour="black"),legend.text=element_text(face="bold",size=13,colour="black"))

#Third condition NC VS NCIL-1
group5 <- read.table("./Scatterplot/分组信息3.txt",sep="\t",row.names=1,header=TRUE,as.is = TRUE)
group6 <- group5[,'group']
tTag3 <- as.data.frame(res3)
g5 <- unique(group6)[1]
g6 <- unique(group6)[2]

g5_exp = NCVSOEIL[rownames(tTag3),rownames(group5)[which(group5$group==g5)]]
g6_exp = NCVSOEIL[rownames(tTag3),rownames(group5)[which(group5$group==g6)]]

g5_mean = apply(g5_exp,1,mean)
g6_mean = apply(g6_exp,1,mean)

type3=rep('No',length(g5_mean))

type3[which(tTag3$log2FoldChange > 0.58 & tTag3$pvalue<0.05)] = "Up"

type3[which(tTag3$log2FoldChange < -0.58& tTag3$pvalue<0.05)] = "Down"

datam3 = data.frame(g5_mean,g6_mean,logFC=tTag3$log2FoldChange,pvalue=tTag3$pvalue,type3,stringsAsFactors=FALSE)

ggplot(datam3,aes(log2(g5_mean),log2(g6_mean),colour=type3))+
  
  geom_point(stat="identity",size=1)+theme(legend.title=element_blank())+scale_color_manual(values =c("Down"='MediumAquamarine',"No"='grey',"Up"='DarkOrange'))+
  
  labs(x=paste(g5,' Log2(TPM)'),y=paste(g6,' Log2(TPM)'),title=paste(g6,' VS ',g5,sep=""))+
  
  coord_cartesian(ylim=c(-10,10),xlim=c(-10,10))+geom_segment(aes(x = -10, y = -10, xend = 10, yend = 10),size=1,colour="#999999",linetype="dotted")+theme(plot.title = element_text(hjust = 0.5),title=element_text(face="bold",size=15,colour="black"),axis.title=element_text(face="bold",size=13,colour="black"),axis.text.x=element_text(face="bold",size=12,colour="black"),axis.text.y=element_text(face="bold",size=12,colour="black"),legend.text=element_text(face="bold",size=13,colour="black"))


