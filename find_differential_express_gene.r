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
avggdata1 <- data.frame("avgdata1"=apply(data_analysis[,1:dim(data1)[2]],1,mean))
avgdata2 <- data.frame("avgdata2"=apply(data_analysis[,(dim(data1)[2]+1):(dim(data_analysis)[2])],1,mean))
log2fc <-  data.frame("log2fc"=log2((avgdata1$avgdata1)/(avgdata2$avgdata2)))
results1 <- cbind(avgdata1,avgdata2,log2fc,p_value,FDR)
rownames(results1)<-rownames(data_analysis)
