load("/fshare2/Rotation/tianchen/sc_rna/hong_data/new/origin-RSEM-tRDS.rdata")
rownames(y)<-sub("EF","DQ",rownames(y))
noGMP<-y[,setdiff(c(1:ncol(y)),grep("GMP",colnames(y)))]
origin_data<-as.matrix(noGMP)

#####move ERCC
panduan<-grepl("^DQ",rownames(origin_data))
panduan<-!panduan
row_new<-rownames(origin_data)[panduan]
origin_data2<-origin_data[row_new,]

###
hsc<-origin_data2

##### ly_mye genes:old genes
load("/fshare2/Rotation/tianchen/new_move_batch_2020/gene_clp.rdata")
load("/fshare2/Rotation/tianchen/new_move_batch_2020/gene_gmp.rdata")
change_rownames<-sub("^ENS.*_+(.*)-.*$","\\1",rownames(hsc))
ly<-rownames(gene_clp)
me<-rownames(gene_gmp)
select<- change_rownames %in% ly
select_ly_gene<-rownames(hsc)[select]

###select_me_genes
select2<- change_rownames %in% me
select_me_gene<-rownames(hsc)[select2]
lym.possig<-select_ly_gene
mye.possig<-select_me_gene
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
save(score_all_cells,file="/fshare2/Rotation/tianchen/score_all_cell/score_all_cells_old_method.rdata")

#### plot score vaule
#### way1:boxplot
clu0_score<-as.data.frame(score_all_cells[intersect(cluster0,rownames(score_all_cells)),])
clu1_score<-as.data.frame(score_all_cells[intersect(cluster1,rownames(score_all_cells)),])

boxplot(clu0_score$lym,clu1_score$lym)
boxplot(clu0_score$mye,clu1_score$mye)

t.test(clu0_score$lym,clu1_score$lym)
t.test(clu0_score$mye,clu1_score$mye)

wilcox.test(clu0_score$mye,clu1_score$mye)
wilcox.test(clu0_score$mye,clu1_score$mye)

##### way2: violin plot
##### 重要的是建立一个data table，group是分组，比如：基因1，基因2，treatment是处理，比如:正常人、患病。
Data_self<-data.frame(Group=factor(c(rep("p7",331),rep("p8",400))),ly_score=c(young_p7_score$lym,young_p8_score$lym),
                      treatment=factor(c(rep("p7",331),rep("p8",400))))
P2<- ggplot(Data_self, aes(x=Group, y=ly_score,fill=treatment)) + 
  geom_violin(trim=F,color="black",size=1,width=1) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.05,position=position_dodge(0.9),fill="white",size=1)+ #绘制箱线图
  scale_fill_manual(values = c("firebrick2", "palegreen3"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 10,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=15),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=15),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("ly_score")+xlab("") #设置x轴和y轴的标题

P2

##### way3: barplot



##### way4: 核密度(density)

