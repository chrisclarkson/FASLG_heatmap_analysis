data=read.table('~/Downloads/data_RNA_Seq_v2_mRNA_median_Zscores.txt',header=T,fill = T,stringsAsFactors = F)
genes=as.vector(data[,1])
data=data[,-c(1,2)]
quantile(data[which(genes=='ENG'),], prob = seq(0, 1, length = 11))
quantile(data[which(genes=='FASLG'),], prob = seq(0, 1, length = 11))
data_up=data[,-which(data[which(genes=='FASLG'),]< 0.34621 & data[which(genes=='ENG'),]< 0.73)]
data_down=data[,-which(data[which(genes=='FASLG'),]> -0.46436 & data[which(genes=='ENG'),]> -0.7837)]

conditions=c('FALSG/ENG down','FALSG/ENG up')
data=cbind(data_down,data_up)
rownames(data)=genes
colnames(data)=c(rep('FALSG/ENG down',ncol(data_down)),rep('FALSG/ENG up',ncol(data_up)))
data=na.omit(data)
library(ComplexHeatmap)
library(circlize)

heat=function(cell){
  genes_OI=read.table(paste0(cell,'.txt'),header = F)
  data_OI=data.frame()
  for(i in genes_OI$V1){
    data_OI=rbind(data_OI,data[which(rownames(data)==i),])
  }
  sum(is.na(data_OI))
  data_OI=apply(data_OI,2,as.numeric)
  rownames(data_OI)=genes_OI$V1
  #data_OI=data_OI[,-apply(data_OI,2,function(x) any(is.na(x)))]
  require("RColorBrewer")
  col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "black", "darkred"))
  tiff(paste0(cell,'_TCGA_patients_ENG_FASLGd.tiff'),width=1000,height=1000)
    print(Heatmap(as.matrix(data_OI),
    col=col_fun,
    cluster_columns = F,
    cluster_rows = F,
    show_column_names=F,
    cluster_column_slices = T,
    show_row_names=T,
    show_row_dend=F,
    column_gap = unit(2, "mm"),
    column_split=colnames(data_OI), 
    heatmap_legend_param = list(title = "Z-scores"),
    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3),
        labels = c("FASLG/ENG low", "FASLG/ENG high"), 
        labels_gp = gpar(col = "white", fontsize = 15)))
   )
  )
  dev.off()
}
    
heat('mesenchymal')
heat('chemokine_receptors')
heat('Macrophage_M2_markers')
