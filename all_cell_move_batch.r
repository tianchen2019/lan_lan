load("/fshare2/Rotation/tianchen/sc_rna/hong_data/new/origin-RSEM-tRDS.rdata")
y<-as.matrix(y)
load("/fshare2/Rotation/tianchen/sc_rna/hong_data/new/aged_hoxb5.rdata")
rownames(y)<-sub("EF","DQ",rownames(y))
rownames(aged_hoxb5)<-sub("EF","DQ",rownames(aged_hoxb5))
previous_data<-data_quhoxb
GMP<-y[,grep("^GMP",colnames(y))]
aged_hoxb5<-as.matrix(aged_hoxb5)
aged_hoxb5<-aged_hoxb5[rownames(y),]

all_data<-as.matrix(cbind(y[,colnames(data_quhoxb)],GMP,aged_hoxb5))
origin_data<-as.matrix(all_data)
rownames(origin_data)<-sub("_","-",rownames(origin_data))
#####move ERCC
panduan<-grepl("^DQ",rownames(origin_data))
panduan<-!panduan
row_new<-rownames(origin_data)[panduan]
origin_data2<-origin_data[row_new,]
dim(origin_data2)

clu_origin<-CreateSeuratObject(counts = origin_data2, project = "all_cell_cluster", min.cells = 5, min.features = 500)
clu_origin[["percent.mt"]] <- PercentageFeatureSet(clu_origin, pattern = "-mt-")
clu_origin[["percent.ercc"]] <- PercentageFeatureSet(clu_origin, pattern = "^DQ")
VlnPlot(clu_origin, features = c("percent.mt"))
VlnPlot(clu_origin, features = c("nFeature_RNA"))
VlnPlot(clu_origin, features = c("percent.ercc"))
VlnPlot(clu_origin, features = "nCount_RNA")
clu_origin <- subset(clu_origin, subset = nFeature_RNA > 500  & nCount_RNA >200  & nCount_RNA < 5000000 & percent.ercc<30  & percent.mt < 10  )


####### use ERCC to move batch effect
DQ<-rownames(origin_data)[grep("^DQ",rownames(origin_data))]
data_qc<-origin_data[c(DQ,rownames(clu_origin@assays$RNA)),colnames(clu_origin@assays$RNA)]
erccs<-(grepl("^DQ",rownames(data_qc)))         ####92
geterccs<-data_qc[erccs,]
erccsum<-colSums(as.matrix(geterccs)) 
normfactor<-100000/erccsum
data_qc<-as.matrix(data_qc)
data_qc2<-as.matrix(sweep(data_qc,2,normfactor,'*'))
colSums(as.matrix(data_qc2)[(grepl("DQ",rownames(data_qc2))),])
####move ercc
panduan<-grepl("^DQ",rownames(data_qc2))
panduan<-!panduan
row_new<-rownames(data_qc2)[panduan]
data_qc3<-data_qc2[row_new,]
dim(data_qc3)

########## use CCA to move Batch effect(Seurat)
qc_data<-data_qc3
AM<-qc_data[,grep("AM",colnames(qc_data))]
GMP<-qc_data[,grep("GMP",colnames(qc_data))]
HOXB<-qc_data[,grep("HOXB5",colnames(qc_data))]
mHSC<-qc_data[,grep("mHSC",colnames(qc_data))]
mMPPa<-qc_data[,grep("mMPPa",colnames(qc_data))]
AH5<-qc_data[,grep("^AH5",colnames(qc_data))]
rahul_data<-as.matrix(cbind(AM,mHSC,HOXB,mMPPa,GMP,AH5))  

Ryoaged<-qc_data[,grep("^Ryoaged",colnames(qc_data))]
Ryoyoung<-qc_data[,grep("^Ryoyoung",colnames(qc_data))]  
ryo_data<-as.matrix(cbind(Ryoaged,Ryoyoung))         
ryo<-CreateSeuratObject(counts = ryo_data,project = "ryo", min.cells = 0, min.features =0)
rahul<-CreateSeuratObject(counts = rahul_data,project = "rahul", min.cells = 0, min.features =0)

rahul$tech<-as.factor(rep("rahul",dim(rahul_data)[2]))
ryo$tech<-as.factor(rep("ryo",dim(ryo_data)[2]))

two_data.list<-list("rahul"=rahul,"ryo"=ryo)
for (i in 1:length(two_data.list)) {
  two_data.list[[i]] <- NormalizeData(two_data.list[[i]], verbose = FALSE)
  two_data.list[[i]] <- FindVariableFeatures(two_data.list[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

reference.list <- two_data.list
two_data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3000)
two_data.integrated <- IntegrateData(anchorset = two_data.anchors, dims = 1:30,features.to.integrate=rownames(qc_data))

DefaultAssay(two_data.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
two_data.integrated <- ScaleData(two_data.integrated, verbose = FALSE)
two_data.integrated <- RunPCA(two_data.integrated, npcs = 30, verbose = FALSE)
DimPlot(two_data.integrated, reduction = "pca", group.by = "tech")
DimPlot(two_data.integrated, reduction = "pca",group.by = "orig.ident",
        cols = c("#2c5770","#2c5770","#ff6c5f","#ff6c5f","#ff6c5f","#ff0092","#ff0092","#598c14","#0099cc","#0099cc",
                 "#0099cc","#0099cc","#05cc47","#05cc47","orange","orange","#7d3f98","#7d3f98","#7d3f98",
                 "#3be8b0","#3be8b0"))


##################tsne################
######################################
two_data.integrated <- JackStraw(two_data.integrated, num.replicate = 100)
two_data.integrated <- ScoreJackStraw(two_data.integrated, dims = 1:20)

two_data.integrated <- FindNeighbors(two_data.integrated, dims = 1:20)
two_data.integrated<- FindClusters(two_data.integrated, resolution = 0.3)
two_data.integrated <- RunTSNE(two_data.integrated, dims.use = 1:10, perplexity = 10)
TSNEPlot(object = two_data.integrated,pt.size = 1.5, group.by = "orig.ident",
         cols = c("#2c5770","#2c5770","#ff6c5f","#ff6c5f","#ff6c5f","#ff0092","#ff0092","#598c14","#0099cc","#0099cc",
                 "#0099cc","#0099cc","#05cc47","#05cc47","orange","orange","#7d3f98","#7d3f98","#7d3f98",
                 "#3be8b0","#3be8b0"))
TSNEPlot(object = two_data.integrated,pt.size = 1.5)

############## check GMP cell ###############
ss<-data.frame("number"=apply(GMP,2,function(x) sum(x>0)))
ss$group<-rep("GMP",dim(ss)[2])
ggplot(ss,aes(x=group, y=number)) + geom_violin(trim=F,color="black",size=1,width=1)+
  geom_boxplot(width=0.05,position=position_dodge(0.9),fill="white",size=1)+
  scale_fill_manual(values="#56B4E9")+ 
  theme_bw()+geom_hline(aes(yintercept=5500))

cell1<-colnames(GMP)[apply(GMP,2,function(x) sum(x>0))<5500      
cell4<-colnames(GMP)[apply(GMP,2,function(x) sum(x>0))>=5500]
GMP_tsne<-as.data.frame(two_data.integrated@reductions$tsne@cell.embeddings)[colnames(two_data.integrated@assays$RNA)[grep("^GMP",colnames(two_data.integrated@assays$RNA))],1:2]

GMP_tsne$color<-rep("#ff6c5f",dim(GMP_tsne)[1])
GMP_tsne[cell1,3]<-rep("#0099cc",length(cell1))    ###blue
GMP_tsne[cell4,3]<-rep("#05cc47",length(cell4))    ###green
plot(GMP_tsne[,1],GMP_tsne[,2],col=GMP_tsne[,3],pch=16,cex=0.8)


