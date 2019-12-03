load("/fshare2/Rotation/tianchen/sc_rna/hong_data/new/origin-RSEM-tRDS.rdata")
noGMP<-y[,setdiff(c(1:ncol(y)),grep("GMP",colnames(y)))]
origin_data<-as.matrix(noGMP)
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

clu1 <- NormalizeData(clu_origin, normalization.method = "LogNormalize", scale.factor = 10000)
clu1<- FindVariableFeatures(clu1,selection.method = "vst", nfeatures = 3000)
###scale
all.genes <- rownames(clu1)
#length(rownames(clu2))
#clu_hgv<- ScaleData(clu2,vars.to.regress = c("nCount_RNA","percent.ercc" ))

clu1<-ScaleData(clu1,features = all.genes)

###PCA
clu1 <- RunPCA(clu1, features = VariableFeatures(object = clu1))
#select_feature<- VariableFeatures(object = clu1)
DimPlot(clu1, reduction = "pca",pt.size = 1)
DimPlot(object = clu1,pt.size = 1,group.by="orig.ident",reduction = "pca",
        cols = c("#ff6c5f","#ff6c5f","#ff6c5f","#ff0092","#ff0092","#0099cc","#0099cc",
                 "#0099cc","#0099cc","#05cc47","#05cc47","orange","orange","#7d3f98","#7d3f98","#7d3f98",
                 "#3be8b0","#3be8b0"))

clu1 <- JackStraw(clu1, num.replicate = 100)
clu1 <- ScoreJackStraw(clu1, dims = 1:20)

clu1 <- FindNeighbors(clu1, dims = 1:20)
clu1 <- FindClusters(clu1, resolution = 1)

clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 20)
clu_origin_weiBatch<-clu_zong
??TSENPlot()
clu1 <- RunTSNE(clu1, dims.use = 1:10, perplexity = 30)
TSNEPlot(object = clu1,pt.size = 1.5)

TSNEPlot(object = clu1,pt.size = 1.5,group.by="orig.ident",
         cols = c("#ff6c5f","#ff6c5f","#ff6c5f","#ff0092","#ff0092","#0099cc","#0099cc",
                  "#0099cc","#0099cc","#05cc47","#05cc47","orange","orange","#7d3f98","#7d3f98","#7d3f98",
                  "#3be8b0","#3be8b0"))
