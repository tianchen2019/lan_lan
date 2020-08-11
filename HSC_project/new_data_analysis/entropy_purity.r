#### read_meta data
load("/fshare2/Rotation/tianchen/young_hoxb5_cell/result/list_integrated.rdata")
meta_data<-read.csv(file="/fshare2/Rotation/tianchen/juptyer/RNOTE/meta.data3.csv",header = T,sep = ",",check.names = F)
rownames(meta_data)<-meta_data$cell
meta_data<-meta_data[,2:23]
umap<-read.csv(file="/fshare2/Rotation/tianchen/juptyer/RNOTE/umap.csv",sep=",",header = T)
rownames(umap)<-rownames(meta_data)
umap<-umap[,2:3]
tsne<-read.csv(file="/fshare2/Rotation/tianchen/juptyer/RNOTE/tsne.csv",sep=",",header = T)
rownames(tsne)<-rownames(meta_data)
tsne<-tsne[,2:3]


data_plot<-data.frame("cell"=rownames(umap),
                      "dataset"=as.factor(sub("(^.*)+?[0-9].*$","\\1",rownames(umap))),
                      "cluster"=as.factor(meta_data$leiden),
                      "umap1"=umap[,1],
                      "umap2"=umap[,2],
                      "hoxb"=rep("other",length(rownames(umap))),
                      "young_hoxb"=rep("other",length(rownames(umap))),
                      "fraction"=rep("other",length(rownames(umap))),
                      "young_fraction"=rep("other",length(rownames(umap))),
                      "mcherry_value"=rep("NA",length(rownames(umap))),
                      "value"=rep(1,length(rownames(umap))),
                      "tsne1"=tsne[,1],
                      "tsne2"=tsne[,2])

load("/fshare2/Rotation/tianchen/young_hoxb5_cell/new_data/RSEMgRDS-afterQC-withercc-readscount-filplate_integratedata.rdata")
sorting_info<-totaldata$prot
qc_matrix<-as.data.frame(t(totaldata$exps))
#rownames(qc_matrix)<-gsub("_","-",rownames(qc_matrix))

ryo_f1<-intersect(rownames(sorting_info)[!is.na(sorting_info$`Fraction1 Events`) & 
                                             sorting_info$`Fraction1 Events`=="1" ],colnames(qc_matrix))
ryo_f2<-intersect(rownames(sorting_info)[!is.na(sorting_info$`Fraction2 Events`) & 
                                             sorting_info$`Fraction2 Events`=="1" ],colnames(qc_matrix))
ryo_f3<-intersect(rownames(sorting_info)[!is.na(sorting_info$`Fraction3 Events`) & 
                                             sorting_info$`Fraction3 Events`=="1" ],colnames(qc_matrix))
hoxb5<-qc_matrix[,grep("^Hoxb5",colnames(qc_matrix))]
young_pHSC_cell<-intersect(rownames(sorting_info)[!is.na(sorting_info$LTHSC_Events) & 
                                                    sorting_info$LTHSC_Events=="0" ],colnames(hoxb5))
young_LTHSC_cell<-intersect(rownames(sorting_info)[!is.na(sorting_info$LTHSC_Events) & 
                                                     sorting_info$LTHSC_Events=="1" ],colnames(hoxb5))
ryoyoung_cell<-colnames(qc_matrix)[c(grep("^RYP2",colnames(qc_matrix)),
                                     grep("^RYP5",colnames(qc_matrix)))]
ryoaged_cell<-colnames(qc_matrix)[c(grep("^RYP1",colnames(qc_matrix)),
                                    grep("^RYP3",colnames(qc_matrix)),
                                    grep("^RYP4",colnames(qc_matrix)))]
data_plot<-as.matrix(data_plot)
rownames(data_plot)<-rownames(umap)
data_plot[ryo_f1,8]<-rep("f1",length(ryo_f1))             
data_plot[ryo_f2,8]<-rep("f2",length(ryo_f2))
data_plot[ryo_f3,8]<-rep("f3",length(ryo_f3))             
#data_plot[all_pHSC_cell,6]<-rep("all_HOXB5N",length(all_pHSC_cell))             
#data_plot[all_LTHSC_cell,6]<-rep("all_LTHSC",length(all_LTHSC_cell))

data_plot[ryoyoung_cell,2]<-rep("ryoyoung",length(ryoyoung_cell))
data_plot[ryoaged_cell,2]<-rep("ryoaged",length(ryoaged_cell))


###
data_plot[young_pHSC_cell,7]<-rep("young_HOXB5N",length(young_pHSC_cell))             
data_plot[young_LTHSC_cell,7]<-rep("young_LTHSC",length(young_LTHSC_cell))

data_plot[intersect(ryoyoung_cell,ryo_f1),9]<-rep("young_f1",length(intersect(ryoyoung_cell,ryo_f1)))  
data_plot[intersect(ryoyoung_cell,ryo_f2),9]<-rep("young_f2",length(intersect(ryoyoung_cell,ryo_f2)))  
data_plot[intersect(ryoyoung_cell,ryo_f3),9]<-rep("young_f3",length(intersect(ryoyoung_cell,ryo_f3)))  

data_point_plot<-as.data.frame(data_plot,stringsAsFactors=F)
data_point_plot[,4]<-as.numeric(data_point_plot[,4])
data_point_plot[,5]<-as.numeric(data_point_plot[,5])
data_point_plot[,12]<-as.numeric(data_point_plot[,12])
data_point_plot[,13]<-as.numeric(data_point_plot[,13])
ggplot(data_point_plot) +
  geom_point(mapping = aes(x = umap1, y = umap2, color = young_hoxb )) +
  scale_color_manual(values = c("grey","#0099cc","#ff6c5f"))+               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot(data_point_plot) +
  geom_point(mapping = aes(x = umap1, y = umap2, color = young_fraction ),size=1,alpha=0.7) +
  scale_color_manual(values = c("grey","#ff6c5f","green","#0099cc"))+               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot(data_point_plot) +
  geom_point(mapping = aes(x = tsne1, y = umap2, color = young_hoxb )) +
  scale_color_manual(values = c("grey","#0099cc","#ff6c5f"))+               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot(data_point_plot) +
  geom_point(mapping = aes(x = tsne1, y = umap2, color = young_fraction ),size=1,alpha=0.7) +
  scale_color_manual(values = c("grey","#ff6c5f","green","#0099cc"))+               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.line = element_line(colour = "black"))



entropy_purity<-function(a,b,c,d){ #a:range,b:tsne_mt,c:x_axis(2) or y_axis(3),d: class list in cluster
  save_entropy<-0
  save_purity<-0
  for (i in 1:(length(a)-1)){
    cell<-b[(b[,c]>=a[i] & b[,c]<a[i+1]),1]
    entropy<-0
    p_list<-0
    for (j in 1:length(d)){
      tmp<-(sum(cell==d[j])/length(cell))*log2((sum(cell==d[j])/length(cell)))*(-1)
      if(is.na(tmp)){
        tmp<-0 
      }else if(tmp==-Inf){
        tmp<-0 
      }else if(tmp!= -Inf){
        tmp<-tmp
      }
      entropy<-entropy+tmp
      p_list[j]<-sum(cell==d[j])/length(cell)
    }
    purity<-max(p_list)
    if(is.na(purity)){
      purity<-0 
    }else if(!is.na(purity)){
      purity<-purity
    }
    #save_entropy<-list()
    #save_purity<-list()
    #save_entropy[[i]]<-entropy
    #save_purity[[i]]<-purity
    save_entropy[i]<-entropy
    save_purity[i]<-purity
  }
  result_data<-data.frame("entropy"=save_entropy,"purity"=save_purity)
  rownames(result_data)<-a
  #return(save_entropy)
  #return(save_purity)
  return(result_data)
  #return(cell)
}

#### entropy_purity
y_range<-seq(-50,10)
tsne<-as.matrix(list.integrated@reductions$tsne@cell.embeddings)
tsne1<-data.frame("cell"=rownames(tsne),"tsne1"=tsne[,1],"tsne2"=tsne[,2],stringsAsFactors = F)
rownames(tsne1)<-tsne1$cell
tsne1[,3]<-as.numeric(tsne1[,3])
tsne1[,2]<-as.numeric(tsne1[,2])


young_tsne<-tsne1[c(intersect(ryoyoung_cell,ryo_f1),intersect(ryoyoung_cell,ryo_f2),
                    intersect(ryoyoung_cell,ryo_f3),young_LTHSC_cell,young_pHSC_cell),]
young_tsne[c(intersect(ryoyoung_cell,ryo_f1),young_LTHSC_cell),1]<-rep("LTHSC",length(c(intersect(ryoyoung_cell,ryo_f1),
                                                                                        young_LTHSC_cell)))
young_tsne[c(intersect(ryoyoung_cell,ryo_f2),young_pHSC_cell),1]<-rep("pHSC",length(c(intersect(ryoyoung_cell,ryo_f2),young_pHSC_cell)))
young_tsne[intersect(ryoyoung_cell,ryo_f3),1]<-rep("F3",length(intersect(ryoyoung_cell,ryo_f3)))

d<-c("LTHSC","pHSC","F3")

result_entropy_purity<-entropy_purity(y_range,young_tsne,3,d)


      

entropy_purity<-function(a,b,c,d,e,f){ #a:y_range,b:tsne_mt,c: y_axis(3),d: class list in cluster,e: x_range,f:x_axis(2) 
  save_entropy<-matrix(0,length(a)-1,(length(e)-1))
  save_purity<-matrix(1,length(a)-1,(length(e)-1))
  result_mt<-list()
  for (i in 1:(length(a)-1)){
    for (k in 1:(length(e)-1)){
      cell<-b[((b[,c]>=a[i] & b[,c]< a[i+1]) & (b[,f]>=e[k] & b[,f]< e[k+1])),1]
      entropy<-0
      p_list<-0
      for (j in 1:length(d)){
        tmp<-(sum(cell==d[j])/length(cell))*log2((sum(cell==d[j])/length(cell)))*(-1)
        if(is.na(tmp)){
          tmp<-0 
        }else if(tmp==-Inf){
          tmp<-0 
        }else if(tmp!= -Inf){
          tmp<-tmp
        }
       entropy<-entropy+tmp
       p_list[j]<-sum(cell==d[j])/length(cell)
      }
     purity<-max(p_list)
     if(is.na(purity)){
       purity<-0 
     }else if(!is.na(purity)){
       purity<-purity
     }
    save_entropy[i,k]<- entropy
    save_purity[i,k]<- purity
    }
  }
  
  result_mt$entropy<-save_entropy
  result_mt$purity<-save_purity
  return(result_mt)
}
x_range<-seq(-56,72)
result_entropy_purity_mt<-entropy_purity(y_range,young_tsne,3,d,x_range,2)
save
library(pheatmap)
p<-pheatmap(log2(result_entropy_purity_mt$entropy+1),display_numbers = TRUE)

save(result_entropy_purity_mt,file="result_entropy_purity_mt.rdata")
pheatmap(log2(result_entropy_purity_mt+1),#gaps_col=c(6,15,24,30,39,48),
         #gaps_row = c(41,56,66,93,136,150,165,186,189,198,206,222,235,240,244),
         breaks = bk,cluster_rows=F,cluster_cols=F,
         border_color="grey46",scale = 'row',
         color = colorRampPalette(c("blue", "white","red"))(15),
         show_rownames = T,show_colnames =F,annotation_col = annotation_col,
         annotation_row = annotation_row_zong,
         annotation_colors=annotation_color,
         fontsize=4,fontsize_row = 1)




