# 设置文件路径
setwd("/fshare2/Rotation/tianchen/young_hoxb5_cell")

#### cell sorting info : 得到细胞的mcherry pos/neg细胞
cell_info<-read.table(file="scanorma_meta_data.csv",header = T,sep=",")
hoxb_cell<-cell_info[grep("^HOXB",cell_info$X),]
hoxb_cell2<-hoxb_cell[!is.na(hoxb_cell$Hoxb),]
pos_cell<-hoxb_cell2[hoxb_cell2$Hoxb=="pos",]
rownames(pos_cell)<-pos_cell$X
neg_cell<-hoxb_cell2[hoxb_cell2$Hoxb=="neg",]
rownames(neg_cell)<-neg_cell$X

### 其他方法 得到mcherry neg/pos的细胞
pos_cell<-hoxb_cell %>% filter(Hoxb=="pos")   
neg_cell<-hoxb_cell %>% filter(Hoxb=="neg")   


#### rna-seq data
load("/fshare2/Rotation/tianchen/sc_rna/hong_data/new/origin-RSEM-tRDS.rdata")
rownames(y)<-sub("EF","DQ",rownames(y))
noGMP<-y[,setdiff(c(1:ncol(y)),grep("GMP",colnames(y)))]
origin_data<-as.matrix(noGMP)
rownames(origin_data)<-sub("_","-",rownames(origin_data))
#####move ERCC
panduan<-grepl("^DQ",rownames(origin_data))
panduan<-!panduan
row_new<-rownames(origin_data)[panduan]
origin_data2<-origin_data[row_new,]
save(origin_data2,file="old_data_noGMP_no_ercc.rdata")

####### 找出mcherry pos/neg 细胞的差异基因
hsc<-origin_data2
hoxb_matrix<-hsc[,intersect(c(rownames(pos_cell),rownames(neg_cell)),colnames(hsc))]

# 去除P2的HOXb.P2细胞进行差异分析。
hsc_matrix_no_P2<-hoxb_matrix[,-c(grep("HOXB5.P2",colnames(hsc_hoxb_cell_mt)))]
hsc_matrix_no_P2<-hsc_matrix_no_P2[,intersect(c(rownames(pos_cell),rownames(neg_cell)),colnames(hsc_matrix_no_P2))]

########### 以下是data后的操作
#### filter:
#### choose data
test_data<-hsc_matrix_no_P2
filter_hoxb <- apply(
  test_data,
  1,
  function(x) length(x[x > 0]) >= 2
)

hoxb_filter <- test_data[filter_hoxb,]
####  hong's way:wilcox
####
myFun <- function(x){
  x = as.numeric(x)
  v1 = x[1:(length(intersect(rownames(pos_cell),colnames(hoxb_filter))))]
  v2 = x[(length(intersect(rownames(pos_cell),colnames(hoxb_filter)))+1):(dim(hoxb_filter)[2])]
  out <- wilcox.test(v1,v2)
  out <- out$p.value
}

p_value <- apply(hoxb_filter,1,myFun)
p_value[is.nan(p_value)] <- 1
FDR <- p.adjust(p_value,method = "fdr")
###method=FDR
avggroup12 <- data.frame("avggroup12"=apply(hoxb_filter[,1:(length(intersect(pos_cell$X,colnames(hoxb_filter))))],1,mean))
avggroup3 <- data.frame("avggroup3"=apply(hoxb_filter[,(length(intersect(pos_cell$X,colnames(hoxb_filter)))+1):(dim(hoxb_filter)[2])],1,mean))
log2fc <-  data.frame("log2fc"=log2((avggroup12$avggroup12)/(avggroup3$avggroup3)))
results1 <- cbind(avggroup12,avggroup3,log2fc,p_value,FDR)
rownames(results1)<-rownames(hoxb_filter)
#### 储存文件
write.table(results1,file = "qu_P2_hoxb_young_pos_neg.tsv",row.names = T, sep="\t", quote=F)
###sig_gene :转录本
sig_gene<-rownames(results1)[results1$p_value<0.05]
### 选择的是基因
sig_gene2<-unique(sub("^ENS.*_+(.*)-.*$","\\1",sig_gene))

#### young_hoxb5_cell use pos_vs_different_gene to plot pca
#####修改数据
zong_hoxb_matrix<-hsc[,grep("^HOXB5.*",colnames(hsc))]
#### 新数据
zong_hoxb_matrix<-as.matrix(newhoxb@raw.data)

####选择输入数据
seurat_matrix<-zong_hoxb_matrix
clu_hoxb<-CreateSeuratObject(counts = seurat_matrix, project = "all_cell_cluster", min.cells = 5, min.features = 500)
clu_hoxb[["percent.mt"]] <- PercentageFeatureSet(clu_hoxb, pattern = "-mt-")
VlnPlot(clu_hoxb, features = c("percent.mt"))
VlnPlot(clu_hoxb, features = c("nFeature_RNA"))
VlnPlot(clu_hoxb, features = "nCount_RNA")
clu_hoxb <- subset(clu_hoxb, subset = nFeature_RNA > 500  & nCount_RNA >200  & nCount_RNA < 5000000  )

clu1 <- NormalizeData(clu_hoxb, normalization.method = "LogNormalize", scale.factor = 10000)
clu1<- FindVariableFeatures(clu1,selection.method = "vst", nfeatures = 3000)
#### change_pca_gene

###scale
all.genes <- rownames(clu1)
#length(rownames(clu2))
#clu_hgv<- ScaleData(clu2,vars.to.regress = c("nCount_RNA","percent.ercc" ))

clu1<-ScaleData(clu1,features = all.genes)

###PCA
clu1 <- RunPCA(clu1, features = VariableFeatures(object = clu1))
#select_feature<- VariableFeatures(object = clu1)
DimPlot(clu1, reduction = "pca",pt.size = 1.5)


clu1 <- JackStraw(clu1, num.replicate = 100)
clu1 <- ScoreJackStraw(clu1, dims = 1:20)

clu1 <- FindNeighbors(clu1, dims = 1:20)
clu1 <- FindClusters(clu1, resolution = 1)

clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 10)

#clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 30)
TSNEPlot(object = clu1,pt.size = 1.5)
TSNEPlot(object = clu1,pt.size = 1.5,group.by="orig.ident")

#### use sig gene to cluster
clu2 <- RunPCA(clu1, features =sig_gene2)
#select_feature<- VariableFeatures(object = clu1)
DimPlot(clu2, reduction = "pca",pt.size = 1.5)
DimPlot(clu2, reduction = "pca",pt.size = 1.5,group.by = "orig.ident")


clu2 <- JackStraw(clu2, num.replicate = 100)
clu2 <- ScoreJackStraw(clu2, dims = 1:20)

clu2 <- FindNeighbors(clu2, dims = 1:20)
clu2 <- FindClusters(clu2, resolution = 1)
clu2 <- RunTSNE(clu2, dims.use = 1:10, perplexity = 10)
#### RUN UMAP
clu2<-RunUMAP(clu2,features = sig_gene2)
DimPlot(clu2, reduction = "umap",pt.size = 1.5,group.by = "orig.ident")

#clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 30)
#TSNEPlot(object = clu2,pt.size = 1.5)
TSNEPlot(object = clu2,pt.size = 1.5,group.by="orig.ident",group.by = "orig.ident")
#### 改变分辨率分cluster
clu3 <- FindClusters(clu3, resolution = 0.3)
clu3 <- RunTSNE(clu3, dims.use = 1:10, perplexity = 10)
TSNEPlot(object = clu3,pt.size = 1.5)

#### find batch effect: move batch effect
data1<-CreateSeuratObject(counts=newhoxb@raw.data[,c(grep("^Hoxb5P1.HOXB5.P1",colnames(newhoxb@raw.data)),
                                                     grep("^Hoxb5p2.H5P2",colnames(newhoxb@raw.data)),
                                                     grep("^Hoxb5P4.HOXB5.P4",colnames(newhoxb@raw.data)))],
                          project = "P1", min.cells = 0, min.features =0)
data3<-CreateSeuratObject(counts=newhoxb@raw.data[,grep("^Hoxb5P3.HOXB5.P3",colnames(newhoxb@raw.data))],
                          project = "P3", min.cells = 0, min.features =0)

two_data.list<-list("P1_P3_P4"=data1,"P3"=data3)
for (i in 1:length(two_data.list)) {
  two_data.list[[i]] <- NormalizeData(two_data.list[[i]], verbose = FALSE)
  two_data.list[[i]] <- FindVariableFeatures(two_data.list[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

reference.list <- two_data.list
two_data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3000)
two_data.integrated <- IntegrateData(anchorset = two_data.anchors, dims = 1:30,features.to.integrate=rownames(data1@assays$RNA))

DefaultAssay(two_data.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
two_data.integrated <- ScaleData(two_data.integrated, verbose = FALSE)
two_data.integrated <- RunPCA(two_data.integrated, npcs = 30, verbose = FALSE)
DimPlot(two_data.integrated, reduction = "pca", group.by = "tech")
DimPlot(two_data.integrated, reduction = "pca",group.by = "orig.ident")
two_data.integrated <- JackStraw(two_data.integrated, num.replicate = 100)
two_data.integrated <- ScoreJackStraw(two_data.integrated, dims = 1:20)

two_data.integrated <- FindNeighbors(two_data.integrated, dims = 1:20)
two_data.integrated<- FindClusters(two_data.integrated, resolution = 0.3)
two_data.integrated <- RunTSNE(two_data.integrated, dims.use = 1:10, perplexity = 10)
TSNEPlot(object = two_data.integrated,pt.size = 1.5, group.by = "orig.ident")
two_data.integrated<-RunUMAP(two_data.integrated,npcs = 30, verbose = FALSE)
DimPlot(two_data.integrated, reduction = "umap",pt.size = 1.5,group.by="orig.ident")

# two_data.integrated<-RunUMAP(two_data.integrated,npcs = 30, verbose = FALSE)
# DimPlot(two_data.integrated, reduction = "umap",pt.size = 1.5,group.by="orig.ident")


###plot pos neg
color<-as.matrix(data.frame("color"=rep("grey",dim(clu2@assays$RNA)[2])))
new_cell_name<-sub("Hoxb5P1.","",colnames(clu2@assays$RNA))
new_cell_name<-sub("Hoxb5P2.","",new_cell_name)
new_cell_name<-sub("Hoxb5P3.","",new_cell_name)
new_cell_name<-sub("Hoxb5P4.","",new_cell_name)

rownames(color)<-new_cell_name
color[intersect(rownames(pos_cell),new_cell_name),1]<-rep("red",length(intersect(rownames(pos_cell),new_cell_name)))
color[intersect(rownames(neg_cell),new_cell_name),1]<-rep("green",length(intersect(rownames(neg_cell),new_cell_name)))
plot(clu2@reductions$tsne@cell.embeddings[,1],clu2@reductions$tsne@cell.embeddings[,2],
     pch=16,col=color[,1],cex=.7)
plot(clu2@reductions$pca@cell.embeddings[,1],clu2@reductions$pca@cell.embeddings[,2],
     pch=16,col=color[,1],cex=.7)
#### plot P7/P8
color2<-as.matrix(data.frame("color"=rep("grey",1204)))
rownames(color2)<-colnames(clu2@assays$RNA)
color2[rownames(young_p8_score),1]<-rep("red",dim(young_p8_score)[1])
color2[rownames(young_p7_score),1]<-rep("green",dim(young_p7_score)[1])
plot(clu2@reductions$tsne@cell.embeddings[,1],clu2@reductions$tsne@cell.embeddings[,2],
     pch=16,col=color2[,1],cex=.7)

plot(clu2@reductions$pca@cell.embeddings[,1],clu2@reductions$pca@cell.embeddings[,2],
     pch=16,col=color2[,1],cex=.7)

dim(young_p7_score)
dim(young_p8_score)

table(color[,1]=="grey")
table(color2[,1]=="grey")





####plot P2 plate
plate<-read.table(file="/fshare2/Rotation/tianchen/young_hoxb5_cell/hoxb5_96cellinfo",sep=" ",header = F)
rownames(plate)<-gsub("-",".",plate$V1)
p2_info<-as.data.frame(plate[grep("^HOXB5.P2",rownames(plate)),])
p10<-rownames(p2_info)[grep("plate10",p2_info$V2)]
p11<-rownames(p2_info)[grep("plate11",p2_info$V2)]
p7<-rownames(p2_info)[grep("plate7",p2_info$V2)]
p8<-rownames(p2_info)[grep("plate8",p2_info$V2)]
p9<-rownames(p2_info)[grep("plate9",p2_info$V2)]

plot(clu1@reductions$pca@cell.embeddings[,1],clu1@reductions$pca@cell.embeddings[,2],
     pch=16,col=color2[,1],cex=.7)

plot(clu1@reductions$tsne@cell.embeddings[,1],clu1@reductions$tsne@cell.embeddings[,2],
     pch=16,col=color2[,1],cex=.7)

table(sub("(^.*)+_.*$","\\1",rownames(young_p7_score)))
table(sub("(^.*)+_.*$","\\1",rownames(young_p8_score)))
HOXB_P2<-rownames(plate)[grep("^HOXB5.P2",rownames(plate))]
table(plate[grep("^HOXB5.P2",plate$V1),2])


color3<-as.matrix(data.frame("color"=rep("grey",279)))
rownames(color3)<-rownames(tsne_matrix[intersect(P2_tsne,HOXB_P2),])
color3[intersect(p10,rownames(color3)),1]<-rep("red",length(intersect(p10,rownames(color3))))
color3[intersect(p11,rownames(color3)),1]<-rep("green",length(intersect(p11,rownames(color3))))
color3[intersect(p7,rownames(color3)),1]<-rep("green",length(intersect(p7,rownames(color3))))
color3[intersect(p8,rownames(color3)),1]<-rep("orange",length(intersect(p8,rownames(color3))))
color3[intersect(p9,rownames(color3)),1]<-rep("black",length(intersect(p9,rownames(color3))))




tsne_matrix<-as.matrix(clu2@reductions$tsne@cell.embeddings)
pca_matrix<-as.matrix(clu2@reductions$pca@cell.embeddings)
P2_tsne<-rownames(tsne_matrix)[grep("^HOXB5.P2",rownames(tsne_matrix))]
plot(tsne_matrix[intersect(P2_tsne,HOXB_P2),1],
     tsne_matrix[intersect(P2_tsne,HOXB_P2),2],
     pch=16,cex=.7,col=color3[,1])


intersect(P2_tsne,HOXB_P2)


####plot mCherry
table(hoxb_cell_seq$mcherry<0)
hoxb_cell_seq$mcherry[59]
color_mcherry<-as.matrix(data.frame(rep("grey",dim(tsne_matrix)[1])))

color_mcherry[1:500,1]<-c(colorRampPalette(colors = c("blue","lightblue"))(500))
color_mcherry[501:1204,1]<-c(colorRampPalette(colors = c("pink","darkred"))(1204-500))


color_mcherry[,1]<-c(colorRampPalette(colors = c("blue","lightblue","darkred"))(1204))

sort(x = hoxb_cell$mcherry)
hoxb_cell_seq<-hoxb_cell[order(hoxb_cell$mcherry),]
tsne_matrix_seq<-tsne_matrix[intersect(hoxb_cell_seq$X,rownames(tsne_matrix)),]
plot(tsne_matrix_seq[,1],
     tsne_matrix_seq[,2],
     pch=16,cex=.7,col=color_mcherry[,1])

pca_matrix_seq<-pca_matrix[intersect(hoxb_cell_seq$X,rownames(pca_matrix)),]
plot(pca_matrix_seq[,1],
     pca_matrix_seq[,2],
     pch=16,cex=.7,col=color_mcherry[,1])
###HOXB5 gene
hoxb_gene<-rownames(clu2@assays$RNA)[grep("Hoxb5",rownames(clu2@assays$RNA))]
FeaturePlot(object = clu2,features = c("ENSMUST00000049272.4-Hoxb5-201","ENSMUST00000150698.1-Hoxb5os-202",
                                       "ENSMUST00000190470.1-Hoxb5-202"))

table(as.matrix(clu2@assays$RNA[hoxb_gene[1],])>0)
table(as.matrix(clu2@assays$RNA[hoxb_gene[2],])>0)
table(as.matrix(clu2@assays$RNA[hoxb_gene[3],])>0)

expree_hoxb_cell<-colnames(clu2@assays$RNA)[as.matrix(clu2@assays$RNA[hoxb_gene[1],])>0 | as.matrix(clu2@assays$RNA[hoxb_gene[2],])>0 | as.matrix(clu2@assays$RNA[hoxb_gene[2],])>0]
color_hoxb_gene<-as.matrix(rep("grey",1204))
rownames(color_hoxb_gene)<-colnames(clu2@assays$RNA)
color_hoxb_gene[expree_hoxb_cell,1]<-rep("red",length(expree_hoxb_cell))

plot(tsne_matrix_seq[,1],
     tsne_matrix_seq[,2],
     pch=16,cex=.7,col=color_hoxb_gene[,1])

as.matrix(clu2@assays$RNA)[hoxb_gene[1],expree_hoxb_cell]




#####去过批次的去做
load("~/tianchen/sc_rna/hong_data/new/youngcellsubtype.Robj")
young3 = UpdateSeuratObject(object = young2)
qubatch_hoxb <- as.matrix(young2@raw.data)[,grep("^HOXB5",colnames(young2@raw.data))]
#select_feature<- VariableFeatures(object = clu1)
clu_qubatch_hoxb<-CreateSeuratObject(counts = qubatch_hoxb, project = "all_cell_cluster", min.cells = 5, min.features = 500)

clu_qubatch_hoxb <- subset(clu_qubatch_hoxb, subset = nFeature_RNA > 500  & nCount_RNA >200  & nCount_RNA < 5000000  )

clu_qubatch_hoxb <- NormalizeData(clu_qubatch_hoxb, normalization.method = "LogNormalize", scale.factor = 10000)
clu_qubatch_hoxb<- FindVariableFeatures(clu_qubatch_hoxb,selection.method = "vst", nfeatures = 3000)
#### change_pca_gene

###scale
all.genes <- rownames(clu_qubatch_hoxb)
#length(rownames(clu2))
#clu_hgv<- ScaleData(clu2,vars.to.regress = c("nCount_RNA","percent.ercc" ))

clu_qubatch_hoxb<-ScaleData(clu_qubatch_hoxb,features = all.genes)
clu_qubatch_hoxb<- RunPCA(clu_qubatch_hoxb, features = VariableFeatures(object = clu_qubatch_hoxb))
DimPlot(clu_qubatch_hoxb, reduction = "pca",pt.size = 1.5)


clu_qubatch_hoxb <- JackStraw(clu_qubatch_hoxb, num.replicate = 100)
clu_qubatch_hoxb <- ScoreJackStraw(clu_qubatch_hoxb, dims = 1:20)

clu_qubatch_hoxb <- FindNeighbors(clu_qubatch_hoxb, dims = 1:20)
clu_qubatch_hoxb <- FindClusters(clu_qubatch_hoxb, resolution = 1)

clu_qubatch_hoxb <- RunTSNE(clu_qubatch_hoxb, dims.use = 1:10, perplexity = 10)

#clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 30)
TSNEPlot(object = clu_qubatch_hoxb,pt.size = 1.5)

TSNEPlot(object = clu_qubatch_hoxb,pt.size = 1.5,group.by="orig.ident")


clu_qubatch_hoxb_sig<-clu_qubatch_hoxb
clu_qubatch_hoxb_sig <- RunPCA(clu_qubatch_hoxb_sig, features = sig_gene)
DimPlot(clu_qubatch_hoxb_sig, reduction = "pca",group.by="orig.ident",pt.size = 1.5)
DimPlot(clu_qubatch_hoxb_sig, reduction = "tsne",pt.size = 1.5)


clu_qubatch_hoxb_sig <- JackStraw(clu_qubatch_hoxb_sig, num.replicate = 100)
clu_qubatch_hoxb_sig <- ScoreJackStraw(clu_qubatch_hoxb_sig, dims = 1:20)

clu_qubatch_hoxb_sig <- FindNeighbors(clu_qubatch_hoxb_sig, dims = 1:20)
clu_qubatch_hoxb_sig <- FindClusters(clu_qubatch_hoxb_sig, resolution = 1)

clu_qubatch_hoxb_sig <- RunTSNE(clu_qubatch_hoxb_sig, dims.use = 1:10, perplexity = 10)

#### 低分辨率:
clu_qubatch_hoxb_sig_2 <- FindClusters(clu_qubatch_hoxb_sig, resolution = 0.3)
clu_qubatch_hoxb_sig_2 <- RunTSNE(clu_qubatch_hoxb_sig_2, dims.use = 1:10, perplexity = 10)
clu_qubatch_hoxb_sig_2<-RunUMAP(clu_qubatch_hoxb_sig_2,features = intersect(sig_gene,rownames(clu_qubatch_hoxb_sig@assays$RNA@counts)))
DimPlot(clu_qubatch_hoxb_sig_2, reduction = "tsne",pt.size = 1.5)


#clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 30)
TSNEPlot(object = clu_qubatch_hoxb_sig,pt.size = 1.5)

TSNEPlot(object = clu_qubatch_hoxb_sig,pt.size = 1.5,group.by="orig.ident")

### umap
clu_qubatch_hoxb_sig<-RunUMAP(clu_qubatch_hoxb_sig,features = intersect(sig_gene,rownames(clu_qubatch_hoxb_sig@assays$RNA@counts)))
DimPlot(clu_qubatch_hoxb_sig, reduction = "umap",group.by="orig.ident",pt.size = 1.5)

DimPlot(clu_qubatch_hoxb_sig, reduction = "umap",pt.size = 1.5)


color_qubatch<-as.matrix(data.frame("color"=rep("grey",dim(clu_qubatch_hoxb_sig@assays$RNA)[2])))
rownames(color_qubatch)<-rownames(clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings)
color_qubatch[intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(pos_cell)),1]<-rep("red",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(pos_cell))))
color_qubatch[intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(neg_cell)),1]<-rep("green",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(neg_cell))))


###p7,p8
color_qubatch[intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(young_p8_score)),1]<-rep("red",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(young_p8_score))))
color_qubatch[intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(young_p7_score)),1]<-rep("green",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(young_p7_score))))



plot(clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings[,1],clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings[,2],
     pch=16,col=color_qubatch[,1],cex=.8)

plot(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings[,1],clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings[,2],
     pch=16,col=color_qubatch[,1],cex=.8)
abline(h=-1,col="red",lty=2)



####umap:
color_qubatch<-as.matrix(data.frame("color"=rep("grey",dim(clu_qubatch_hoxb_sig@assays$RNA)[2])))
rownames(color_qubatch)<-rownames(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings)
color_qubatch[intersect(rownames(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings),rownames(pos_cell)),1]<-rep("red",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(pos_cell))))
color_qubatch[intersect(rownames(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings),rownames(neg_cell)),1]<-rep("green",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(neg_cell))))


###p7,p8
color_qubatch[intersect(rownames(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings),rownames(young_p8_score)),1]<-rep("red",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(young_p8_score))))
color_qubatch[intersect(rownames(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings),rownames(young_p7_score)),1]<-rep("green",length(intersect(colnames(clu_qubatch_hoxb_sig@assays$RNA),rownames(young_p7_score))))



plot(clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings[,1],clu_qubatch_hoxb_sig@reductions$umap@cell.embeddings[,2],
     pch=16,col=color_qubatch[,1],cex=.8)


###统计比例:
table(clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings[,2]>-1)
shang<-rownames(clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings)[clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings[,2]>-1]
xia<-setdiff(rownames(clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings),shang)
length(intersect(shang,rownames(young_p8_score)))
length(intersect(shang,rownames(young_p7_score)))
#####shang: 396 ,p8:253(63.89%),p7:84(21.2%),other:59

length(intersect(xia,rownames(young_p8_score)))
length(intersect(xia,rownames(young_p7_score)))
##### xia: 398, p8:134(33.67%),p7:229(57.54%)

##### mcherry值

tsne_clu_qubatch_hoxb_sig_seq<-clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings[intersect(hoxb_cell_seq$X,rownames(clu_qubatch_hoxb_sig@reductions$tsne@cell.embeddings)),]
table(hoxb_cell_seq$mcherry<0)
hoxb_cell_seq$mcherry[59]
color_mcherry<-as.matrix(data.frame(rep("grey",dim(tsne_matrix)[1])))

length(intersect(hoxb_cell_seq$X[hoxb_cell_seq$mcherry<0],rownames(tsne_clu_qubatch_hoxb_sig_seq)))

color_mcherry[1:300,1]<-c(colorRampPalette(colors = c("blue","lightblue"))(300))
color_mcherry[301:794,1]<-c(colorRampPalette(colors = c("pink","darkred"))(794-300))


plot(tsne_clu_qubatch_hoxb_sig_seq[,1],tsne_clu_qubatch_hoxb_sig_seq[,2],
     pch=16,col=color_mcherry[,1],cex=.8)



###hoxb5 expression
rownames(clu_qubatch_hoxb_sig@assays$RNA@counts)[grep("Hoxb5",rownames(clu_qubatch_hoxb_sig@assays$RNA@counts))]
FeaturePlot(object = clu_qubatch_hoxb_sig,features = c("ENSMUST00000049272.4-Hoxb5-201","ENSMUST00000150698.1-Hoxb5os-202",
                                       "ENSMUST00000190470.1-Hoxb5-202"))

#### marker
markers1 <- FindAllMarkers(clu_qubatch_hoxb_sig_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers1 %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
write.table(markers1,file = "markers1_young_hoxb5.tsv",row.names = F, sep="\t", quote=F)

VlnPlot(object = clu_qubatch_hoxb_sig_2, features =c("ENSMUST00000023707.10-SOD1-201","ENSMUST00000079812.7-NOTCH2-201"), ncol = 3)

top15 <- markers1 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

DoHeatmap(clu_qubatch_hoxb_sig_2, features = top15$gene, size = 6)+
  theme(
    axis.text.y = element_text(size = 10,
                               family = "Times",face = "bold",colour = "black"))+
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high ="red",
    name = "Expression"
  )

FeaturePlot(object =clu_qubatch_hoxb_sig_2, 
            features =c("ENSMUST00000023707.10-SOD1-201","ENSMUST00000079812.7-NOTCH2-201"),
                        ncol = 3)
length(markers1$gene[markers1$cluster=="1"])
length(markers1$gene[markers1$cluster=="0"])
length(markers1$gene[markers1$cluster=="2"])
length(markers1$gene[markers1$cluster=="3"])



cell_marker<-read.table(file="/fshare2/Rotation/tianchen/new_scRNA/data/GEXC-SurfaceProteinAtlas.csv",header=T,sep="\t")
ss1<-sub("^ENS.*-+(.*)-.*$","\\1",markers1$gene)
ss2<-sub("^ENS.*-+(.*)-.*$","\\1",markers1$gene[markers1$cluster=="1"])

toupper(x=c("apple", "I like your style"))
length(intersect(tolower(cell_marker$GeneSymbol),ss))
length(intersect(tolower(ss2),tolower(cell_marker$GeneSymbol)))

markers1$gene[markers1$cluster=="1"]


FeaturePlot(object =clu_qubatch_hoxb_sig_2, 
            features =c("ENSMUST00000001566.9-TUBB5-201","ENSMUST00000079812.7-NOTCH2-201","ENSMUST00000032754.8-SEMA4B-201",
"ENSMUST00000021864.7-SSR1-201","ENSMUST00000020640.7-RACK1-201","ENSMUST00000139959.7-RACK1-205",
"ENSMUST00000084838.13-CD47-201","ENSMUST00000067880.12-ADAM10-201","ENSMUST00000095078.2-LRRC8A-201",
"ENSMUST00000148993.4-KIT-205","ENSMUST00000050166.13-ADAM22-202","ENSMUST00000070690.7-PTAFR-201",
"ENSMUST00000094361.10-HSP90AA1-202","ENSMUST00000021698.12-HSP90AA1-201"))

VlnPlot(object = clu_qubatch_hoxb_sig_2, features =c("ENSMUST00000001566.9-TUBB5-201","ENSMUST00000079812.7-NOTCH2-201","ENSMUST00000032754.8-SEMA4B-201",
                                                     "ENSMUST00000021864.7-SSR1-201","ENSMUST00000020640.7-RACK1-201","ENSMUST00000139959.7-RACK1-205"
                                                   ),ncol=3)


VlnPlot(object = clu_qubatch_hoxb_sig_2, features =c("ENSMUST00000084838.13-CD47-201","ENSMUST00000067880.12-ADAM10-201",
                                                     "ENSMUST00000095078.2-LRRC8A-201",
                                                     "ENSMUST00000148993.4-KIT-205","ENSMUST00000050166.13-ADAM22-202","ENSMUST00000070690.7-PTAFR-201",
                                                     "ENSMUST00000094361.10-HSP90AA1-202","ENSMUST00000021698.12-HSP90AA1-201"),ncol = 4)

save(clu_qubatch_hoxb_sig_2,file="/fshare2/Rotation/tianchen/young_hoxb5_cell/result/young_aged_hoxb5_cluster.rdata")


####### mcherry_vs_mye
mcherry_vaule<-data.frame("mcherry"=cell_info$mcherry[!is.na(cell_info$mcherry)])
rownames(mcherry_vaule)<-cell_info$X[!is.na(cell_info$mcherry)]
mcherry_vaule <- transform(mcherry_vaule,"mye_vaule"=score_all_cells[rownames(mcherry_vaule),2],"ly_vaule"=score_all_cells[rownames(mcherry_vaule),1])
cor.test(mcherry_vaule$mcherry,mcherry_vaule$mye_vaule)
cor.test(mcherry_vaule$mcherry,mcherry_vaule$ly_vaule)

plot(mcherry_vaule$mcherry,mcherry_vaule$mye_vaule,pch=16,cex=.6,xlab="mcherry_value",ylab="mye_value")
plot(mcherry_vaule$mcherry,mcherry_vaule$ly_vaule,pch=16,cex=.6,xlab="mcherry_value",ylab="ly_value")

plot(mcherry_vaule$mye_vaule,mcherry_vaule$ly_vaule,pch=16,cex=.6,xlab="mcherry_value",ylab="ly_value")

