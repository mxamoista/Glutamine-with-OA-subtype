library(ComplexHeatmap)
CIBERSORT1 <- read.table("CIBERSORT-Results.txt",header = TRUE , row.names =1)
CIBERSORT1<- CIBERSORT1[,(colnames(CIBERSORT1) %in% c("ProC","HomC","HTC","RegC","preHTC","FC","EC"))]# select clusters


ann_colorsDeconvolution <- list(group = c("Negative Control"="#8A2BE2", "GLSOE"="#FFC0CB","Negative Control+IL-1b"="#00BFFF","GLSOE+IL-1b"="#008000"))

color.palette3 =c(colorRampPalette(c("#efedf5","#ffffff","#fff5eb","#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704"	))(100))

annotationGLSOEIL_col = data.frame(
  group = c(rep("Negative Control",3),rep("GLSOE",3),rep("Negative Control+IL-1b",3),rep("GLSOE+IL-1b",3)))

ha_columnDeconvolution = rowAnnotation(df = annotationGLSOEIL_col,col=ann_colorsDeconvolution,gap = unit(0.3, "mm"),annotation_height = 0.35,annotation_name_gp = gpar(col = "black")) 

ht2 = Heatmap(t((t(CIBERSORT1))), cluster_columns= F,cluster_rows = F, show_row_names = F ,show_row_dend=T,show_column_names=T,
              #              ,column_dend_gp=,column_dend_reorder=
              name="Normalized abundance",col=color.palette3, left_annotation = ha_columnDeconvolution)