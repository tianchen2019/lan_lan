##### move batch effect (CCA) 通用code
## 提供 data1_mt.data2_mt
data1<-CreateSeuratObject(counts=data1_mt,
                          project = "data1", min.cells = 0, min.features =0)
data2<-CreateSeuratObject(counts=data2_mt,
                          project = "data2", min.cells = 0, min.features =0)

list<-list("data1"=data1,"data2"=data2
           #"data3"=...
           )
for (i in 1:length(list)) {
  list[[i]] <- NormalizeData(list[[i]], verbose = FALSE)
  list[[i]] <- FindVariableFeatures(list[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}
reference.list <- list
list.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3000)
lsit.integrated <- IntegrateData(anchorset = list.anchors, dims = 1:30,features.to.integrate=rownames(data1@assays$RNA))
DefaultAssay(list.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
list.integrated <- ScaleData(list.integrated, verbose = FALSE)
list.integrated <- RunPCA(list.integrated, npcs = 30, verbose = FALSE)
DimPlot(list.integrated, reduction = "pca", group.by = "tech")
DimPlot(list.integrated, reduction = "pca",group.by = "orig.ident")
list.integrated <- JackStraw(list.integrated, num.replicate = 100)
list.integrated <- ScoreJackStraw(list.integrated, dims = 1:20)
list.integrated <- FindNeighbors(list.integrated, dims = 1:20)
list.integrated<- FindClusters(list.integrated, resolution = 0.3)
list.integrated <- RunTSNE(list.integrated, dims.use = 1:10, perplexity = 10)
TSNEPlot(object = list.integrated,pt.size = 1.5, group.by = "orig.ident")
list.integrated<-RunUMAP(list.integrated,npcs = 30, verbose = FALSE)
DimPlot(list.integrated, reduction = "umap",pt.size = 1.5,group.by="orig.ident")

### save file
