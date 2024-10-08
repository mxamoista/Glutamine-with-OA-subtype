
library(DESeq2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("DESeq2")
# Construct matrix 
NegativeControl <- tpmtotal[,(colnames(tpmtotal) %in% c("NC1","NC2","NC3")) ]  # 
OE <- tpmtotal[,(colnames(tpmtotal) %in% c("OE1","OE2","OE3")) ]
NCIL1B <- tpmtotal[,(colnames(tpmtotal) %in% c("NC.IL.1b1","NC.IL.1b2","NC.IL.1b3")) ]
OEIL1B <- tpmtotal[,(colnames(tpmtotal) %in% c("OE.IL.1b1","OE.IL.1b2","OE.IL.1b3")) ]

#First condition NCIL-1 VS OEIL-1
NCILVSOEIL<- cbind(NCIL1B,OEIL1B)

condition1 <- factor(c(rep("Negative Control+IL-1b",3),rep("GLSOE+IL-1b",3)))
colData1 <- data.frame(row.names=colnames(NCILVSOEIL), condition1)
matrix.round1<-round(NCILVSOEIL)
dds1 <- DESeqDataSetFromMatrix(countData = matrix.round1, colData = colData1, design = ~ condition1)
res1 <- DESeq(dds1, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
res1 <- results(res1,contrast = c('condition1', 'GLSOE+IL-1b', 'Negative Control+IL-1b'))
summary(res1) 
res1 <- data.frame(res1, stringsAsFactors = FALSE, check.names = FALSE)# res格式转化：用data.frame转化为表格形式
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]# 依次按照pvalue值log2FoldChange值进行排序
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 0.58 ),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -0.58 ),]    # 表达量显著下降的基因
write.table(res1,file="res1.txt",quote = F,sep="\t")
res1_total <- rbind(res1_up,res1_down)
#根据res1 做火山图
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
res1$Pvalue2<--log10(as.matrix(res1$pvalue))
res1$Geneid<-rownames(res1)
library(dplyr)
res1=distinct(res1,Geneid,.keep_all = T)

res1$threshold = factor(ifelse( abs(res1$log2FoldChange) >= 0.58, ifelse(res1$log2FoldChange>= 0.58,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ggplot(res1,aes(x=log2FoldChange,y=Pvalue2,color=threshold))+
  geom_point(size=2)+
  scale_color_manual(values=c("yellow","blue","grey"))+#ȷ????????ɫ
  geom_text_repel(
    data =res1[c('Mmp3','Adamts5','Gls','Got2'),],
    aes(label =Geneid),fontface="bold", color="black", 
    box.padding=unit(2, "lines"), size=4.5,
    segment.colour = "grey50")+ theme_classic(base_size = 30)+#???ӹ?ע?ĵ??Ļ?????
  theme_bw()+#?޸?ͼƬ????
  theme(
    legend.title = element_blank()#????ʾͼ??????
  )+
  ylab('-LOG(P-value)')+#?޸?y??????
  xlab('log2FoldChange')#?޸?x??????



#Second condition NC VS NCIL-1
NCVSNCIL<- cbind(NegativeControl,NCIL1B)
condition2 <- factor(c(rep("Negative Control",3),rep("Negative Control+IL-1b",3)))
colData2 <- data.frame(row.names=colnames(NCVSNCIL), condition2)
matrix.round2<-round(NCVSNCIL)
dds2 <- DESeqDataSetFromMatrix(countData = matrix.round2, colData = colData2, design = ~ condition2)
res2 <- DESeq(dds2, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res2 <- results(res2,contrast = c('condition2', 'Negative Control+IL-1b', 'Negative Control'))
summary(res2) 
res2 <- data.frame(res2, stringsAsFactors = FALSE, check.names = FALSE)# res格式转化：用data.frame转化为表格形式
res2 <- res2[order(res2$pvalue, res2$log2FoldChange, decreasing = c(FALSE, TRUE)), ]# 依次按照pvalue值log2FoldChange值进行排序
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res2_up<- res2[which(res2$log2FoldChange >= 0.58 & res2$pvalue<0.05),]      # 表达量显著上升的基因
res2_down<- res2[which(res2$log2FoldChange <= -0.58 & res2$pvalue<0.05),]    # 表达量显著下降的基因
res2$Pvalue2<--log10(as.matrix(res2$pvalue))
res2$Geneid<-rownames(res2)
library(dplyr)
res2=distinct(res2,Geneid,.keep_all = T)
res2$threshold = factor(ifelse(res2$'pvalue'<0.05 & abs(res2$log2FoldChange) >= 0.58, ifelse(res2$log2FoldChange>= 0.58,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ggplot(res2,aes(x=log2FoldChange,y=Pvalue2,color=threshold))+
  geom_point(size=2)+
  scale_color_manual(values=c("yellow","blue","grey"))+#ȷ????????ɫ
  geom_text_repel(
    data =res2[c('Col2a1','Sox9','Slc1a5'),],
    aes(label =Geneid),fontface="bold", color="black", size=4.5,
    box.padding=unit(2, "lines"),
    segment.colour = "grey50")+ theme_classic(base_size = 30)+#???ӹ?ע?ĵ??Ļ?????
  theme_bw()+#?޸?ͼƬ????
  theme(
    legend.title = element_blank()#????ʾͼ??????
  )+
  ylab('-LOG(P-value)')+#?޸?y??????
  xlab('log2FoldChange')#?޸?x??????


#Third condition NC VS OEIL-1
NCVSOEIL<- cbind(NegativeControl,OEIL1B)
condition3 <- factor(c(rep("Negative Control",3),rep("GLSOE+IL-1b",3)))
colData3 <- data.frame(row.names=colnames(NCVSOEIL), condition3)
matrix.round3<-round(NCVSOEIL)
dds3 <- DESeqDataSetFromMatrix(countData = matrix.round3, colData = colData3, design = ~ condition3)
res3 <- DESeq(dds3, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
res3 <- results(res3,contrast = c('condition3', 'GLSOE+IL-1b', 'Negative Control'))
summary(res3) 
res3 <- data.frame(res3, stringsAsFactors = FALSE, check.names = FALSE)# res格式转化：用data.frame转化为表格形式
res3 <- res3[order(res3$pvalue, res3$log2FoldChange, decreasing = c(FALSE, TRUE)), ]# 依次按照pvalue值log2FoldChange值进行排序
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res3_up<- res3[which(res3$log2FoldChange >= 0.58 & res3$pvalue<0.05),]      # 表达量显著上升的基因
res3_down<- res3[which(res3$log2FoldChange <= -0.58 & res3$pvalue<0.05),]    # 表达量显著下降的基因
res3$Pvalue2<--log10(as.matrix(res3$pvalue))
res3$Geneid<-rownames(res3)
library(dplyr)
res3=distinct(res3,Geneid,.keep_all = T)
res3$threshold = factor(ifelse(res3$'pvalue'<0.05 & abs(res3$log2FoldChange) >= 0.58, ifelse(res3$log2FoldChange>= 0.58,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ggplot(res3,aes(x=log2FoldChange,y=Pvalue2,color=threshold))+
  geom_point(size=2)+
  scale_color_manual(values=c("yellow","blue","grey"))+#ȷ????????ɫ
  geom_text_repel(
    data =res3[c('Col2a1','Sox9','Slc1a5'),],
    aes(label =Geneid),fontface="bold", color="black", size=4.5,
    box.padding=unit(2, "lines"),
    segment.colour = "grey50")+ theme_classic(base_size = 30)+#???ӹ?ע?ĵ??Ļ?????
  theme_bw()+#?޸?ͼƬ????
  theme(
    legend.title = element_blank()#????ʾͼ??????
  )+
  ylab('-LOG(P-value)')+#?޸?y??????
  xlab('log2FoldChange')#?޸?x??????

#Fourth condition NC VS OE
NCVSOE <- cbind(NegativeControl,OE)
condition4 <- factor(c(rep("Negative Control",3),rep("GLSOE",3)))
colData4<- data.frame(row.names=colnames(NCVSOE), condition4)
matrix.round4<-round(NCVSOE)
dds4 <- DESeqDataSetFromMatrix(countData = matrix.round4, colData = colData4, design = ~ condition4)
res4 <- DESeq(dds4, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
res4 <- results(res4,contrast = c('condition4', 'GLSOE', 'Negative Control'))
summary(res4) 
res4 <- data.frame(res4, stringsAsFactors = FALSE, check.names = FALSE)# res格式转化：用data.frame转化为表格形式
res4 <- res4[order(res4$pvalue, res4$log2FoldChange, decreasing = c(FALSE, TRUE)), ]# 依次按照pvalue值log2FoldChange值进行排序
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res4_up<- res4[which(res4$log2FoldChange >= 0.58 & res4$pvalue<0.05),]      # 表达量显著上升的基因
res4_down<- res4[which(res4$log2FoldChange <= -0.58 & res4$pvalue<0.05),]    # 表达量显著下降的基因
res4$Pvalue2<--log10(as.matrix(res4$pvalue))
res4$Geneid<-rownames(res4)
library(dplyr)
res4=distinct(res4,Geneid,.keep_all = T)
res4$threshold = factor(ifelse(abs(res4$log2FoldChange) >= 0.58, ifelse(res4$log2FoldChange>= 0.58,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ggplot(res4,aes(x=log2FoldChange,y=Pvalue2,color=threshold))+
  geom_point(size=2)+
  scale_color_manual(values=c("yellow","blue","grey"))+#ȷ????????ɫ
  geom_text_repel(
    data =res4[c('Mmp13','Atp6v0e2'),],
    aes(label =Geneid),fontface="bold", color="black", size=4.5,
    box.padding=unit(2, "lines"),
    segment.colour = "grey50")+ theme_classic(base_size = 30)+#???ӹ?ע?ĵ??Ļ?????
  theme_bw()+#?޸?ͼƬ????
  theme(
    legend.title = element_blank()#????ʾͼ??????
  )+
  ylab('-LOG(P-value)')+#?޸?y??????
  xlab('log2FoldChange')#?޸?x??????

