####### use ERCC to move batch effect
qc_data<-data_quhoxb
DQ<-rownames(origin_data)[grep("^DQ",rownames(origin_data))]
data_qc<-origin_data[c(DQ,rownames(qc_data)),colnames(qc_data)]
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

rahul_data<-as.matrix(cbind(AM,mHSC,HOXB,mMPPa))    
ryo_data<-as.matrix(cbind(Ryoaged,Ryoyoung))         
ryo<-CreateSeuratObject(counts = ryo_data,project = "ryo", min.cells = 0, min.features =0)
rahul<-CreateSeuratObject(counts = rahul_data,project = "rahul", min.cells = 0, min.features =0)

rahul$tech<-as.factor(rep("rahul",dim(rahul)[2]))
ryo$tech<-as.factor(rep("ryo",dim(ryo)[2]))

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
        cols = c("#ff6c5f","#ff6c5f","#ff6c5f","#ff0092","#ff0092","#0099cc","#0099cc",
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
         cols = c("#ff6c5f","#ff6c5f","#ff6c5f","#ff0092","#ff0092","#0099cc","#0099cc",
                 "#0099cc","#0099cc","#05cc47","#05cc47","orange","orange","#7d3f98","#7d3f98","#7d3f98",
                 "#3be8b0","#3be8b0"))
TSNEPlot(object = two_data.integrated,pt.size = 1.5)
