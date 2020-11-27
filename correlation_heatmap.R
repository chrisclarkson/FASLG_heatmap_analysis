data=read.table('~/Downloads/data_RNA_Seq_v2_mRNA_median_Zscores.txt',header=T,fill = T,stringsAsFactors = F)
genes=as.vector(data[,1])
data=data[,-c(1,2)]
data=na.omit(data)

data_OI=data.frame()
genes_OI=c('NGFR','ENG','CD33','CD86','CD163','CLEC7A','CSF1','HIF1A','IL33','MRC1','TLR1','TLR8','CCL2','IL6','IL1A','IL1B')


for(i in genes_OI){
	data_OI=rbind(data_OI,data[which(genes==i),])
}
sum(is.na(data_OI))

data_OI=apply(data_OI,2,as.numeric)
rownames(data_OI)=genes_OI
# library(psych)
# correl=psych::cor.test(t(data_OI))
# library(Hmisc)
# correl <- Hmisc::rcorr(as.matrix(t(data_OI)), type="pearson")$P
correl=cor(t(data_OI))
colnames(correl)=genes_OI
rownames(correl)=genes_OI
library(ComplexHeatmap)
tiff(paste0('correlations.tiff'))
    print(Heatmap(as.matrix(correl),
    cluster_columns = T,
    cluster_rows = T,
    show_column_names=T,
    cluster_column_slices = F,
    show_row_names=T,
    show_row_dend=T,
    heatmap_legend_param = list(title = "Pvalues"),
  ))
  dev.off()

