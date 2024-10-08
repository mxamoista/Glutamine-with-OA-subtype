library(ComplexHeatmap)
Total<- read.csv(file ="Total.csv", row.names=1) 
kbCar <- Total$Length / 1000
countCar <- Total[,1:12]
rpkCar <-countCar / kbCar 
tpmtotal <- t(t(rpkCar)/colSums(rpkCar) * 1000000)
write.table(tpmtotal,file="tpmtotal.txt",quote = F,sep = "\t") 
tpmCar02 <- tpmtotal[(rownames(tpmtotal) %in% c("Gls","Mmp3","Adamts5","Col2a1","SOX9","Slc1a5","Got2","Atp6v0e2")), ]  # select genes

annotationGLSOEIL_col = data.frame(
  group = c(rep("Negative Control",3),rep("GLSOE",3),rep("Negative Control+IL-1b",3),rep("GLSOE+IL-1b",3)))
color.palette =c(colorRampPalette(c("#87CEFA","#B0E0E6","#efedf5","#fff5eb","#F08080","#FF0000"))(100))

col.C1=rgb(31, 149, 203, max = 255)
col.C2=rgb(255, 233, 0, max = 255)
col.C4=rgb(236, 95, 92, max = 255)
col.C3=rgb(236,171,85, max = 255)


annotation_col = data.frame(group=as.factor(annotationGLSOEIL_col$group)) 

ann_colors <- list(
  group = c("Negative Control"="#8A2BE2", "GLSOE"="#FFC0CB","Negative Control+IL-1b"="#00BFFF","GLSOE+IL-1b"="#008000"))

ha_column = HeatmapAnnotation(df = annotation_col, col = ann_colors,gap = unit(0.5, "mm"),annotation_height = 0.35,annotation_name_gp = gpar(col = "black")) 

#Dset <- Dset[order(apply(Dset,1,var),decreasing = T),][1:900,]


ht1 = Heatmap(t(scale(t(tpmCar02))), cluster_columns= T, cluster_rows= T, show_row_names = T ,show_row_dend=F,show_column_names=F,
              #              ,column_dend_gp=,column_dend_reorder=
              name=" ",col=color.palette, top_annotation = ha_column)
