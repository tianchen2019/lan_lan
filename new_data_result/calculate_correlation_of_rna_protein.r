### input data
qcdata<-read.table(file="/fshare2/Rotation/tianchen/young_hoxb5_cell/new_data/readcountsperENSMUSGene-afterQC-v1.tsv",header = T,sep=",",check.names = F)
rownames(qcdata)<-qcdata$gene
qc_matrix<-as.matrix(qcdata[,-1])
sorting_info<-read.table(file="/fshare2/Rotation/tianchen/young_hoxb5_cell/new_data/mergedsorting.csv",header = T,check.names = F,sep=",")
rownames(sorting_info)<-sorting_info$cell

#### calculate correlation of RNA and protein
gene_expression<-qc_matrix[,intersect(rownames(sorting_info),colnames(qc_matrix))]
gene_expression<-as.matrix(gene_expression)
rownames(gene_expression)<-toupper(rownames(gene_expression))
sorting_info2<-sorting_info[intersect(rownames(sorting_info),colnames(qc_matrix)),]
sorting_info2<-as.data.frame(sorting_info2)
colnames(sorting_info2)<-toupper(colnames(sorting_info2))
colnames_gene<-sub("^.*_+?(.*)+_.*$","\\1",colnames(sorting_info2))
genes<-intersect(toupper(rownames(qc_matrix)),colnames_gene)
need_col<-colnames(sorting_info2)[colnames_gene %in% genes]

plot_list_pos=list()
cor_val<-list()
cor_pval<-list()
for (i in 1:length(need_col)){
  cells<-rownames(sorting_info2)[!is.na(sorting_info2[,need_col[i]])]
  gene<-sub("^.*_+?(.*)+_.*$","\\1",need_col[i])
  cor_val[i]=c(cor.test(gene_expression[gene,cells],sorting_info2[cells,need_col[i]])$estimate)
  cor_pval[i]=c(cor.test(gene_expression[gene,cells],sorting_info2[cells,need_col[i]])$p.value)
  names(cor_val)[i]<-need_col[i]
  names(cor_pval)[i]<-need_col[i]
  plot(gene_expression[gene,cells],sorting_info2[cells,need_col[i]],pch=16,
       ylab=paste(need_col[i],"protein",sep=" of "),xlab=paste(need_col[i],"RNA",sep=" of "))
  #p1<-plot(gene_expression["Il7r",],sorting_info2$`pHSC_IL7R_Alexa Fluor 700-A Mean`,pch=16,
       #ylab=paste(gene,"protein level",sep=" "),xlab=paste(gene,"RNA level",sep=" "))
  #plot_list_pos[[i]]=p1
}

plot_list_pos=list()
for(j in 1:length(plot_neg_gene[1:20])){
  p1<-FeaturePlot(object = list.integrated, dims = c(1,2),features = plot_neg_gene[j],reduction ="tsne",
                  cols = c("lightgrey", "red"),pt.size=1)
  plot_list_pos[[j]]=p1
}

for (j in 1:length(plot_list_pos)){
  file_name1=paste(j,"_neg",".png",sep="")
  png(file_name1)
  print(plot_list_pos[[j]])
  dev.off()
}
