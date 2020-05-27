#### calculate score value
load("/fshare2/Rotation/tianchen/new_move_batch_2020/gene_clp.rdata")
load("/fshare2/Rotation/tianchen/new_move_batch_2020/gene_gmp.rdata")
load("/fshare2/Rotation/tianchen/young_hoxb5_cell/result/hsc_7254_integrate.rdata")
lym.possig<-rownames(gene_clp)
mye.possig<-rownames(gene_gmp)
hsc<-hsc_7254_integrate$raw_data
fin.siggene <- matrix(0,length(lym.possig)+length(mye.possig),2)
colnames(fin.siggene) <- c("lym","mye")
rownames(fin.siggene) <- c(1:dim(fin.siggene)[1])

rownames(fin.siggene)[1:length(lym.possig)] <- lym.possig
fin.siggene[1:length(lym.possig),1] <- rep(1,length(lym.possig))
rownames(fin.siggene)[(length(lym.possig)+1):(length(lym.possig)+length(mye.possig))] <- mye.possig
fin.siggene[(length(lym.possig)+1):(length(lym.possig)+length(mye.possig)),2] <- rep(1,length(mye.possig))

pos.gene.mat <- fin.siggene
select.mhsc <- hsc[intersect(rownames(hsc),rownames(fin.siggene)),] 
select.mhsc[select.mhsc>0] <- 1 
select.posgene.mat <- pos.gene.mat[intersect(rownames(select.mhsc),rownames(pos.gene.mat)),] 
numpro <- 2
pos.score <- matrix(0,dim(select.mhsc)[2],numpro)
colnames(pos.score) <- colnames(select.posgene.mat)[1:numpro]
rownames(pos.score) <- colnames(select.mhsc)
#i:progenitor
for(i in 1:numpro){
  pos.gene.mask <- rownames(select.posgene.mat)[select.posgene.mat[,i]==1]
  pos.score[,i] <- colSums(select.mhsc[pos.gene.mask,])
}
pos.gene.sum <- colSums(select.posgene.mat)
hsc.expressed.gene.amount <- colSums(select.mhsc[rownames(select.posgene.mat),])
new.pos.score <- pos.score/hsc.expressed.gene.amount
hsc.expressed.gene.amount <- colSums(select.mhsc[rownames(select.posgene.mat),])
new.pos.score <- pos.score/hsc.expressed.gene.amount
score_all_cells <- t(t(new.pos.score)/pos.gene.sum[1:numpro])
