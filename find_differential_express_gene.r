data1<-newhoxb@raw.data[,grep("^Hoxb5P4",colnames(newhoxb@raw.data))]
data2<-newhscfil@raw.data[,grep("^Hoxb5p4",colnames(newhscfil@raw.data))]
data_analysis<-cbind(data1[intersect(rownames(data1),rownames(data2)),],data2[intersect(rownames(data1),rownames(data2)),])
myFun <- function(x){
  x = as.numeric(x)
  v1 = x[1:dim(data1)[2]]
  v2 = x[(dim(data1)[2]+1):(dim(data_analysis)[2])]
  out <- wilcox.test(v1,v2)
  out <- out$p.value
}

p_value <- apply(data_analysis,1,myFun)
p_value[is.nan(p_value)] <- 1
FDR <- p.adjust(p_value,method = "fdr")
###method=FDR
avgdata1 <- data.frame("avgdata1"=apply(data_analysis[,1:dim(data1)[2]],1,mean))
avgdata2 <- data.frame("avgdata2"=apply(data_analysis[,(dim(data1)[2]+1):(dim(data_analysis)[2])],1,mean))
log2fc <-  data.frame("log2fc"=log2((avgdata1$avgdata1)/(avgdata2$avgdata2)))
results1 <- cbind(avgdata1,avgdata2,log2fc,p_value,FDR)
rownames(results1)<-rownames(data_analysis)

### 三组之间比较的code
gene_expression<-read.table(file="genes.fpkm_tracking_bulk_rnaseq_mhsc",header=T)
rownames(gene_expression)<-gene_expression$tracking_id
high<-gene_expression[,c("HSC_MH2_B6J_FPKM","HSC_MH2_B6J_FPKM","HSC_FH2_B6J_FPKM","HSC_FH1_B6J_FPKM")]
low<-gene_expression[,c("HSC_ML1_B6J_FPKM","HSC_ML2_B6J_FPKM","HSC_FL1_B6J_FPKM","HSC_FL2_B6J_FPKM")]
neg<-gene_expression[,c("HSC_MN1_B6J_FPKM","HSC_MN2_B6J_FPKM","HSC_FN1_B6J_FPKM","HSC_FN2_B6J_FPKM")]
data1<-high
data2<-low
data3<-neg
gene1<- intersect(rownames(data1),rownames(data2))
gene_final<-intersect(gene1,rownames(data3))
data_analysis<-cbind(data1[gene_final,],data2[gene_final,],data3[gene_final,])

myFun1 <- function(x){
  x = as.numeric(x)
  v1 = x[1:dim(data1)[2]]
  v2 = x[(dim(data1)[2]+1):(dim(data1)[2]+dim(data2)[2])]
  out <- wilcox.test(v1,v2)
  out <- out$p.value
  }
myFun2 <- function(x){
  x = as.numeric(x)
  v1 = x[1:dim(data1)[2]]
  v3 = x[(dim(data1)[2]+dim(data2)[2]+1):(dim(data_analysis)[2])]
  out <- wilcox.test(v1,v3)
  out <- out$p.value
}
myFun3 <- function(x){
  x = as.numeric(x)
  v2 = x[(dim(data1)[2]+1):(dim(data1)[2]+dim(data2)[2])]
  v3 = x[(dim(data1)[2]+dim(data2)[2]+1):(dim(data_analysis)[2])]
  out <- wilcox.test(v2,v3)
  out <- out$p.value
}
p_value1 <- apply(data_analysis,1,myFun1)
p_value1[is.nan(p_value1)] <- 1
FDR1 <- p.adjust(p_value1,method = "fdr")
avgdata1 <- data.frame("avgdata1"=apply(data_analysis[,1:dim(data1)[2]],1,mean))
avgdata2 <- data.frame("avgdata2"=apply(data_analysis[,(dim(data1)[2]+1):(dim(data1)[2]+dim(data2)[2])],1,mean))
avgdata3 <- data.frame("avgdata3"=apply(data_analysis[,(dim(data1)[2]+dim(data2)[2]+1):(dim(data_analysis)[2])],1,mean))
log2fc1 <-  data.frame("log2fc1"=log2((avgdata1$avgdata1+1)/(avgdata2$avgdata2+1)))

p_value2 <- apply(data_analysis,1,myFun2)
p_value2[is.nan(p_value2)] <- 1
FDR2 <- p.adjust(p_value2,method = "fdr")
###method=FDR
log2fc2 <-  data.frame("log2fc2"=log2((avgdata1$avgdata1+1)/(avgdata3$avgdata3+1)))

p_value3 <- apply(data_analysis,1,myFun3)
p_value3[is.nan(p_value3)] <- 1
FDR3 <- p.adjust(p_value3,method = "fdr")
###method=FDR
log2fc3 <-  data.frame("log2fc3"=log2((avgdata2$avgdata2+1)/(avgdata3$avgdata3+1)))

results <- cbind(avgdata1,avgdata2,avgdata2,log2fc1,p_value1,FDR1,log2fc2,p_value2,FDR2,log2fc3,p_value3,FDR3)
rownames(results)<-rownames(data_analysis)

#### select gene
data1_sig_gene1<-rownames(results)[results$p_value1<0.05 & results$log2fc1>0]
data1_sig_gene2<-rownames(results)[results$p_value2<0.05 & results$log2fc2>0]
data1_sig_gene<-unique(c(data1_sig_gene1,data1_sig_gene2))

data2_sig_gene1<-rownames(results)[results$p_value1<0.05 & results$log2fc1<0]
data2_sig_gene2<-rownames(results)[results$p_value3<0.05 & results$log2fc3>0]
data2_sig_gene<-unique(c(data2_sig_gene1,data2_sig_gene2))

data3_sig_gene1<-rownames(results)[results$p_value2<0.05 & results$log2fc1<0]
data3_sig_gene2<-rownames(results)[results$p_value3<0.05 & results$log2fc3<0]
data3_sig_gene<-unique(c(data3_sig_gene1,data3_sig_gene2))

list_3_class_2<-list("high"=setdiff(data1_sig_gene,c(data2_sig_gene,data3_sig_gene)),
                   "low"=setdiff(data2_sig_gene,c(data1_sig_gene,data3_sig_gene)),
                   "neg"=setdiff(data3_sig_gene,c(data1_sig_gene,data2_sig_gene)))

#### select fc
high_fc_mean<-apply(results[list_3_class_2$high,c(4,7)],1,mean)
high_fc_mean<-data.frame("gene"=list_3_class_2$high,"mean_fc"=high_fc_mean)
rownames(high_fc_mean)<-list_3_class_2$high

low_fc_mean<-apply(abs(results[list_3_class_2$low,c(4,10)]),1,mean)
low_fc_mean<-data.frame("gene"=list_3_class_2$low,"mean_fc"=low_fc_mean)
rownames(low_fc_mean)<-list_3_class_2$low

neg_fc_mean<-apply(abs(results[list_3_class_2$neg,c(7,10)]),1,mean)
neg_fc_mean<-data.frame("gene"=list_3_class_2$neg,"mean_fc"=neg_fc_mean)
rownames(neg_fc_mean)<-list_3_class_2$neg

####rank
high_gene_rank<-arrange(high_fc_mean,desc(mean_fc))
low_gene_rank<-arrange(low_fc_mean,desc(mean_fc))
neg_gene_rank<-arrange(neg_fc_mean,desc(mean_fc))

list_3_class_rank<-list("high"=high_gene_rank$gene,
                        "low"=low_gene_rank$gene,
                        "neg"=neg_gene_rank$gene)

####high+low vs neg sig gene
myFun4 <- function(x){
  x = as.numeric(x)
  v1_v2 = x[1:(dim(data1)[2]+dim(data2)[2])]
  v3 = x[(dim(data1)[2]+dim(data2)[2]+1):(dim(data_analysis)[2])]
  out <- wilcox.test(v1_v2,v3)
  out <- out$p.value
}
p_value4 <- apply(data_analysis,1,myFun1)
p_value4[is.nan(p_value4)] <- 1
FDR4 <- p.adjust(p_value4,method = "fdr")
avgdata1_2 <- data.frame("avgdata1_2"=apply(data_analysis[,1:(dim(data1)[2]+dim(data2)[2])],1,mean))
avgdata3 <- data.frame("avgdata3"=apply(data_analysis[,(dim(data1)[2]+dim(data2)[2]+1):(dim(data_analysis)[2])],1,mean))
log2fc4 <-  data.frame("log2fc"=log2((avgdata1_2$avgdata1_2+1)/(avgdata3$avgdata3+1)))
results2<-cbind(avgdata1_2,avgdata3,log2fc4,p_value4,FDR4)
rownames(results2)<-rownames(data_analysis)
data1_2_sig_gene<-rownames(results2)[p_value4<0.05 & log2fc4>0]
data3_sig_gene<-rownames(results2)[p_value4<0.05 & log2fc4<0]
list_2_class<-list("high_low"=data1_2_sig_gene,"neg"=data3_sig_gene)
high_low_sig_gene_mt<-data.frame("gene"=list_2_class$high_low,"fc"=results2[list_2_class$high_low,3])
high_low_gene_rank<-arrange(high_low_sig_gene_mt,desc(fc))

neg_sig_gene_mt<-data.frame("gene"=list_2_class$neg,"fc"=abs(results2[list_2_class$neg,3]))
neg_2_class_gene_rank<-arrange(neg_sig_gene_mt,desc(fc))
list_2_class_rank<-list("high_low"=high_low_gene_rank$gene,
                        "neg"=neg_2_class_gene_rank$gene)

##### gene featureplot
load("/Share2/home/lanxun11/tianchen_rotation/young_hoxb_analysis/result_move_batch/all_cell_move_batch4.rdata")
plot_high_gene<-intersect(list_3_class_rank$high,rownames(list.integrated@assays$RNA))
plot_low_gene<-intersect(list_3_class_rank$low,rownames(list.integrated@assays$RNA))
plot_neg_gene<-intersect(list_3_class_rank$neg,rownames(list.integrated@assays$RNA))

plot_list_pos=list()
for(j in 1:length(plot_high_gene[1:20])){
	p1<-FeaturePlot(object = list.integrated, dims = c(1,2),features = plot_high_gene[j],reduction ="tsne",
    cols = c("lightgrey", "red"),pt.size=1)
    plot_list_pos[[j]]=p1
}

for (j in 1:length(plot_list_pos)){
  file_name1=paste(j,"_high",".png",sep="")
  png(file_name1)
  print(plot_list_pos[[j]])
  dev.off()
}

###low_plot
plot_list_pos=list()
for(j in 1:length(plot_low_gene[1:20])){
	p1<-FeaturePlot(object = list.integrated, dims = c(1,2),features = plot_low_gene[j],reduction ="tsne",
    cols = c("lightgrey", "red"),pt.size=1)
    plot_list_pos[[j]]=p1
}

for (j in 1:length(plot_list_pos)){
  file_name1=paste(j,"_low",".png",sep="")
  png(file_name1)
  print(plot_list_pos[[j]])
  dev.off()
}
#### neg
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



