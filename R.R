#############数据合并，批次校正#####################
#加载包
library(limma)
library(sva)
#setwd("D:\\")     #设置工作目录
files=c("TCGA_LUAD_TPMsymbol.txt","GSE30219","GSE3141","GSE50081", "GSE31210.txt","GSE72094.txt")     #输入文件名称

#获取交集基因
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1]!="TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)

#################差异基因####################
#得到差异-lncRNA,差异-mRNA。
##DESeq2

####差异mRNA与外泌体相关mRNA手动取交集得到30个差异的外泌体mRNA。

####################30个差异的外泌体mRNA与437个差异lncRNA做相关性分析#################
library(corrplot)                #引用包
library(ggcorrplot)
inputFile="merge_cancer_30+437_exp.txt"       #输入文件，行为基因名，列为样本名

#读取输入文件，并对数据进行整理
rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)    #读取输入文件
rt=t(rt)
geneNum=ncol(rt)                   #基因数目

#相关性分析
M=cor(rt)
p.mat <- cor_pmat(rt)

M=as.matrix(M)
p.mat=as.matrix(p.mat)

write.table(M, file="merge_cancer_cor.cancer.txt", sep="\t", quote=F, col.names=T)
write.table(p.mat, file="merge_cancer_Pvalue.cancer.txt", sep="\t", quote=F, col.names=T)

#得到的134个差异lncRNA

############################预后分析(单因素cox)############################
#手动合并表达数据和生存数据，第一列为样品名，第二列为生存时间，第三列为生存状态，之后为基因名

library(survival)    

rt=read.table("merge_cancer_exp_time.txt",header=T,sep="\t",check.names=F,row.names=1)   

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="merge_lncRNA_uniCox.txt",sep="\t",row.names=F,quote=F)

##得到19个预后相关的ER-lncRNAs
####热图####
library(pheatmap)        
cellFile="merge_cancer_19_exp.txt"     
scoreFile="scores.txt"   #########免疫微环境评分，因后面单独展示，所以最终用AI去掉了。            
clusterFile="cluster.txt"         
cliFile="clinical.txt"               
#setwd("")     

cell=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
cell=log2(cell+1)

score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)

#
sameSample=intersect(row.names(cell), row.names(score))

cell=cell[sameSample,]
score=score[sameSample,]
immData=cbind(cell, score[,1:2], ICIcluster=cluster[sameSample,])
Project=gsub("(.*?)\\_.*", "\\1", rownames(immData))
rownames(immData)=gsub("(.*?)\\_(.*?)", "\\2", rownames(immData))
immData=cbind(immData, Project)

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(immData), row.names(cli))
immData=immData[sameSample,]
cli=cli[sameSample,]

#
data=cbind(immData, cli)
data=data[order(data$ICIcluster),]

Type=data[,(22:ncol(data))]
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(Type$ICIcluster))
Type$ICIcluster=letter[match(Type$ICIcluster, uniqClu)]
data=t(data[,1:21])

#
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
ICIcol=bioCol[1:length(levels(factor(Type$ICIcluster)))]
names(ICIcol)=levels(factor(Type$ICIcluster))
ann_colors[["ICIcluster"]]=ICIcol

#
pdf("19geneheatmap.pdf", height=5, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

###########聚类分析############

library(ConsensusClusterPlus)        #引用包
cellFile="merge_cancer_19_exp.txt"     #表达输入文件
workDir="D:\\"     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


#输出结果
clusterNum=2        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)

####################聚类后生存分析##############
#引用包
library(survival)
library(survminer)
clusterFile="cluster.txt"     #分类
cliFile="merge_time.txt"               #生存数据文件
#setwd("")      #设置工作目录

#读取输入文件
cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/12

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,], ICIcluster=cluster[sameSample,])
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(rt$ICIcluster))
rt$ICIcluster=letter[match(rt$ICIcluster, uniqClu)]

#生存差异统计
length=length(levels(factor(rt$ICIcluster)))
diff=survdiff(Surv(futime, fustat) ~ ICIcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ ICIcluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="ICI cluster",
                   legend.labs=levels(factor(rt[,"ICIcluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()

#########PCA分析#######
#pca analysis
library(limma)
#setwd("")                                              #设置工作目录
rt=read.table("merge_cancer_19_exp.txt",sep="\t",header=T,check.names=F)                             #读取输入文件
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

type=sapply(strsplit(colnames(data),"\\-"),"[",4)
type=sapply(strsplit(type,""),"[",1)
type=gsub("2","1",type)
data=t(data)


data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)                                  #PCA分析
write.table(predict(data.pca),file="newTab.xls",quote=F,sep="\t")        #输出新表

#可视化
library(ggplot2)
cluster=read.table("cluster.txt",sep="\t",header=F)                      #读取分型文件
group=paste0("",as.vector(cluster[,2]))
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)

pdf(file="PCA.pdf",height=5,width=6.5)             #保存输入出文件
ggplot(data = PCA,aes(PCA1, PCA2))+geom_point(aes(color = group)) + 
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

###############免疫细胞浸润ssGSEA分析################
library(GSVA)
library(limma)
library(GSEABase)
expFile="merge_cancer.txt"    
gmtFile="immune.gmt"      
#setwd("")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
dim(mat)	

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
dim(geneSet)	

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut, file="29ssgseaOut.txt", sep="\t", quote=F, col.names=F)

###################estimate肿瘤微环境#################
library(estimate)        
inputFile="merge_cancer.txt"     
#setwd("")     

filterCommonGenes(input.f=inputFile, 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")


estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")


scores=read.table("estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores=scores[,1:3]
out=rbind(ID=colnames(scores), scores)
write.table(out, file="scores.txt", sep="\t", quote=F, col.names=F)


####################箱式图####################

library(reshape2)
library(ggpubr)
cellFile="29ssgseaOut.txt"     
scoreFile="scores.txt"               
clusterFile="cluster.txt"         
#setwd("")    

cell=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(cell), row.names(score))
cell=cell[sameSample,]
score=score[sameSample,]
immData=cbind(cell, score[,1:2], ICIcluster=cluster[sameSample,])
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(immData$ICIcluster))
immData$ICIcluster=letter[match(immData$ICIcluster, uniqClu)]
immData[,-ncol(immData)]=scale(immData[,-ncol(immData)])

data=melt(immData, id.vars=c("ICIcluster"))
colnames(data)=c("ICIcluster", "Immune", "Fraction")
data$Fraction[data$Fraction>12]=12

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(immData[,"ICIcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="ICIcluster", 
            ylab="Scale of fraction",
            xlab="",
            legend.title="cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="29boxplot.pdf", width=8, height=6.5)                          #????ͼƬ?ļ?
p+stat_compare_means(aes(group=ICIcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

################免疫基因差异##########
library(limma)
library(ggplot2)
library(ggpubr)
expFile="merge_cancer.txt"             
clusterFile="cluster.txt"    
gene="PDCD1LG2"                    
showName="PDCD1LG2"                
#setwd("D")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=rbind(data,gene=data[gene,])
data=t(data[c(gene,"gene"),])

cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,], ICIcluster=cluster[sameSample,])
data=as.data.frame(data)
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(data$ICIcluster))
data$ICIcluster=letter[match(data$ICIcluster, uniqClu)]

group=levels(factor(data$ICIcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"ICIcluster"])))]
pdf(file="PDCD1LG2vioplot.pdf", width=6, height=5)
ggviolin(data, x="ICIcluster", y="gene", fill = "ICIcluster",
         ylab=paste0(showName, " expression"),
         xlab="",
         legend.title="cluster",
         palette=bioCol, 
         add="boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

############################预后模型构建############################
##1.分型间差异基因
library(limma)                
expFile="merge_cancer.txt"           
cluFile="cluster.txt"      
logFCfilter=1                 
adj.P.Val.Filter=0.05        
#setwd("D:\")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

cluster=read.table(cluFile, header=F, sep="\t", check.names=F, row.names=1)

sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample,]

letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster))
Type=letter[match(cluster, uniqClu)]
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()
for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  #print(contrast)
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #
  diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
  allDiffGenes=c(allDiffGenes, row.names(diffSig))
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
}

#
uniqGene=unique(allDiffGenes)
write.table(uniqGene, file="allDiffGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)


##2.分型间预后基因 #######先手动合并生存时间，生存状态，差异基因表达量
library(survival)    
rt=read.table("merge_cancer_time.txt",header=T,sep="\t",check.names=F,row.names=1)   

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="merge_cancer_uniCox.txt",sep="\t",row.names=F,quote=F)


##3.分型间差异-预后基因GO和KEGG分析
#GOKEGG
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

install.packages("R.utils")

R.utils::setOption("clusterProfiler.download.method",'auto') 

rt=read.table("1641venngene.txt",sep="\t",header=T,check.names=F)           #读取id.txt文件

# 读取基因列表
# 如果你的基因列表是gene symbol，可以使用如下函数转换，将SYMBOL转换为基因id
gene.df <- bitr(rt$gene,
                fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)# 物种是人类
head(gene.df)

gene=gene.df$ENTREZID

#GO富集分析
go <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果

#柱状图
pdf(file="GO.pdf",width = 10,height = 12)
barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#kegg富集分析
kegg <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #富集分析
write.table(kegg,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果

showNum=30
if(nrow(kegg)<showNum){
  showNum=nrow(kegg)
}

pdf(file="keggbarplot.pdf",width = 9,height = 7)
barplot(kegg, drop = TRUE, showCategory = showNum, color = "qvalue")
dev.off()

#4.LASSO
library("glmnet")
library("survival")

coxSigFile="merge_sva261.txt"       #输入文件，第一列为样本名，第二列为生存时间，第三列为生存状态，之后为261个基因
geneFile="261.txt"         #基因文件
#setwd("D:")    #设置工作目录
rt=read.table(coxSigFile,header=T,sep="\t",row.names=1)              #读取文件
geneRT=read.table(geneFile,header=F,sep="\t",check.names=F)          #读取基因文件

GSE31210=read.table("GSE31210.expTime.txt",header=T,sep="\t",row.names=1)

sameGene=intersect(as.vector(geneRT[,1]),colnames(GSE31210))

rt=rt[,c("futime","fustat",sameGene)]

#构建模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)

#输出相关基因系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

#输出train组风险值
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="mergeRisk.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值GSE31210
rt=read.table("GSE31210.expTime.txt",header=T,sep="\t",row.names=1)
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="GSE31210Risk.txt",sep="\t",quote=F,row.names=F)

##输出test组风险值-GSE76427
rt=read.table("GSE72094.expTime.txt",header=T,sep="\t",row.names=1)
rt=rt[,lassoGene]
rt=log2(rt+1)

testFinalGeneExp=rt[,lassoGene]

testScore=apply(testFinalGeneExp,1,myFun)

outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="GSE72094Risk.txt",sep="\t",quote=F,row.names=F)

#######5.生存分析
#引用包
library(survival)
library(survminer)
#setwd("")
bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                   #读取输入文件
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=paste0("p=",pValue),
                     pval.size=5,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table=F,
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 5,height =4.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="mergeRisk.txt",outFile="mergeRisk.survival.pdf")
bioSurvival(inputFile="GSE31210Risk.txt",outFile="GSE31210Risk.survival.pdf")
bioSurvival(inputFile="GSE72094Risk.txt",outFile="GSE72094Risk.survival.pdf")

###############ROC#################
#引用包
library(survival)
library(survminer)
library(timeROC)
#定义绘制ROC曲线函数
bioROC=function(inputFile=null,rocFile=null){
  #读取输入文件
  rt=read.table(inputFile,header=T,sep="\t")
  ###ROC曲线
  ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
}

bioROC(inputFile="mergeRisk.txt",rocFile="merge.ROC.pdf")
bioROC(inputFile="GSE31210Risk.txt",rocFile="GSE31210Risk.ROC.pdf")
bioROC(inputFile="GSE72094Risk.txt",rocFile="GSE72094Risk.ROC.pdf")

###########################模型评估#############################
##1.riskplot

library(pheatmap)
#setwd("")         

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)   #??ȡ?????ļ?
  rt=rt[order(rt$riskScore),]    
  
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 8,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("blue",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
  dev.off()
  
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="blue"
  pdf(file=survStatFile,width = 8,height = 6)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
}
bioRiskPlot(inputFile="mergeRisk.txt",riskScoreFile="mergeRisk.suv.pdf",survStatFile="mergeRisk.survStat.pdf")
bioRiskPlot(inputFile="GSE31210Risk.txt",riskScoreFile="GSE31210Risk.suv.pdf",survStatFile="GSE31210Risk.survStat.pdf")
bioRiskPlot(inputFile="GSE72094Risk.txt",riskScoreFile="GSE72094Risk.suv.pdf",survStatFile="GSE72094Risk.survStat.pdf")

##2.PCA

library(Rtsne)
library(ggplot2)
#setwd("")
bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
  
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)
  data=rt[c(3:(ncol(rt)-2))]
  risk=rt[,"risk"]
  
  
  data.pca=prcomp(data, scale. = TRUE)
  pcaPredict=predict(data.pca)
  PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
  #
  pdf(file=pcaFile, height=4.5, width=5.5)   
  p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
    scale_colour_manual(name="Risk",  values =c("red", "blue"))+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
  
  
  
  tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
  tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	
  #
  pdf(file=tsneFile, height=4.5, width=5.5)     
  p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = risk)) +
    scale_colour_manual(name="Risk",  values =c("red", "blue"))+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
}
bioPCA(inputFile="mergeRisk.txt", pcaFile="mergeRisk.PCA.pdf", tsneFile="mergeRisk.t-SNE.pdf")
bioPCA(inputFile="GSE31210Risk.txt", pcaFile="GSE31210Risk.PCA.pdf", tsneFile="GSE31210Risk.t-SNE.pdf")
bioPCA(inputFile="GSE72094Risk.txt", pcaFile="GSE72094Risk.PCA.pdf", tsneFile="GSE72094Risk.t-SNE.pdf")

##3.独立预后因素分析，输入文件分别为mergeRisk，GSE72094Risk

library(survival)    
#setwd("")    
risk=read.table("mergeRisk.txt",header=T,sep="\t",check.names=F,row.names=1)   
cli=read.table("mergeRisk_cli.txt",sep="\t",check.names=F,header=T,row.names=1)    


sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])



uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="mergeRisk_uniCox.txt",sep="\t",row.names=F,quote=F)



uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="mergeRisk_multiCox.txt",sep="\t",row.names=F,quote=F)



bioForest=function(coxFile=null,forestFile=null,height=null){
  
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  
  pdf(file=forestFile, width = 6.3,height = height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, "red", "blue")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}


bioForest(coxFile="mergeRisk_uniCox.txt",forestFile="mergeRisk_uniCoxForest.pdf",height=4.5)
bioForest(coxFile="mergeRisk_multiCox.txt",forestFile="mergeRisk_multiCoxForest.pdf",height=3.5)


##########################列线图#########################

library(survival)
library(regplot)
library(rms)
library(ggplot2)


riskFile="mergeRisk.txt"       #
cliFile="merge_cli.txt"       #每列为临床病理参数
#setwd("")                    #

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)


samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)


res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[1,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)


nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)


pdf(file="calibration.pdf", width=5, height=5)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

#########DCA##########

library(survival)
library(survminer)
library(timeROC)
library(ggDCA)

predictTime=1     #预测时间
riskFile="nomoRisk.txt"        #
cliFile="merge_cli_DCA.txt"         #
#setwd("")     #

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli1=cli[samSample,,drop=F]
data=cbind(risk1, cli1)

#DCA
rt=cbind(risk1[,c("futime","fustat","risk","Nomogram")], cli1)
rt[,"Age"]=ifelse(rt[,"Age"]>65, 1, 0)
rt[,"Nomogram"]=ifelse(rt[,"Nomogram"]>median(rt[,"Nomogram"]), 1, 0)
Nomogram<-coxph(Surv(futime,fustat)~Nomogram,rt)
Risk<-coxph(Surv(futime,fustat)~risk,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
Gender<-coxph(Surv(futime,fustat)~Gender,rt)
Stage<-coxph(Surv(futime,fustat)~Stage,rt)

###DCA图
pdf(file="DCA.pdf", width=6.5, height=5.2)
d_train=dca(Nomogram,Risk,Age,Gender,Stage,times=predictTime)
ggplot(d_train, linetype=1)
print(ggplot)
dev.off()


######ROC
rt=cbind(risk1[,c("futime","fustat","riskScore","Nomogram")], cli1)
aucText=c()
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)
pdf(file="cliROC.pdf", width=6, height=6)


#
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

#
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}


#
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

#################化疗药物敏感性分析#####################
#1Paclitaxel，2Gemcitabine，3Cisplatin，4Docetaxel，5Gefitinib，6Erlotinib，
#7Methotrexate，8Parthenolide，9Rapamycin，10Vinblastine

#####药物列表见2016文件######
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
expFile="merge_cancer.txt"     #表达输入文件
riskFile="cluster.txt"      #分型输入文件
#riskFile="risk_group.txt"      #风险输入文件
drug="Vinblastine"         #药物名称，需要修改
#setwd("")    #设置工作目录

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#预测药物敏感性
senstivity=pRRopheticPredict(data, drug, selection=1,dataset = "cgp2016")
senstivity=senstivity[senstivity!="NaN"]

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#风险文件和药物敏感性合并
sameSample=intersect(row.names(risk), names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample]
rt=cbind(risk, senstivity)
write.table(rt,file="Vinblastine_risk.txt",sep="\t",quote=F,row.names = T)                          #保存结果

#设置比较组
rt$risk=factor(rt$risk, levels=c("A", "B","C"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
                  xlab="cluster",
                  ylab=paste0(drug, " senstivity (IC50)"),
                  legend.title="cluster",
                  palette = "nejm"
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file=paste0(drug, ".pdf"), width=5, height=4.5)
print(boxplot)
dev.off()



