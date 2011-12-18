# analysis of GSE27003
library("doMC")
library("LiblineaR")
library("stringr")
library("Rutilities")
library(edgeR)
library(GenomicFeatures)
# function 
bioc=function(){source("http://www.bioconductor.org/biocLite.R")}
tryRead=function(x,i=5,n=-1){read.table(x,header=TRUE, nrow=n, stringsAsFactors=FALSE)[i]}
tryRead1=function(x,i=1,n=-1){read.table(x,header=TRUE, nrow=n, stringsAsFactors=FALSE)[i]}

tryRead6=function(x,i=6,n=-1){y=read.table(x,header=TRUE, nrow=n, stringsAsFactors=FALSE)[c(1,i)]; values=y[,2];names(values)=y[,1];return(values)}

normalRead=function(x,n=-1,header=TRUE, ...){read.table(x,header=header, nrow=n, stringsAsFactors=FALSE,...)}

rowNoZero=function(m){return(m[rowSums(m>0)==ncol(m),])}
rowFilter=function(m,n){return(m[rowSums(m>n)==ncol(m),])}
topStd=function(x,n){std=apply(x, 1, sd); return(as.matrix(x[which(std>sort(std,decreasing=TRUE)[n]-1),]))}
cluster=function(x, rect=FALSE, k=3, border='red', dist.method="euclidean"){d=dist(scale(t(x)), dist.method);hc=hclust(d);plclust(hc,hang=-1);re=rect.hclust(hc,k,border=border)}
# print the number of genes pass multiple hypothesis testing
mhtGeneNum=function(topTags.object,adjust.method="fdr",alpha=0.05){return(length(which(topTags(de.com, n=length(topTags.object$genes[[1]]), adjust.method=adjust.method)[[1]]$adj.P.Val < alpha)))}
mhtGene=function(topTags.object,adjust.method="fdr",alpha=0.05){
  gn=topTags(de.com, n=length(topTags.object$genes[[1]]), adjust.method=adjust.method)[[1]]$genes
  gn.sig=gn[which(topTags(de.com, n=length(topTags.object$genes[[1]]), adjust.method=adjust.method)[[1]]$adj.P.Val < alpha)]
  return(gn.sig)
  }

# valcano plot
volcanoPlot=function(MPA,m.cut=1,p.cut=0.05,p.transform=log10,ylab="-log 10 Adjust P Value", main=NULL){
  M=MPA[,1]; p=MPA[,2]; A=MPA[,3]
  plot(M, -p.transform(p), pch = 16, col = "darkgrey", main=main, xlab = "log Fold Change", ylab = ylab)   
  points(M[which(M>m.cut  & A<p.cut)], -log10(p[which(M>m.cut  & A<p.cut)]), col = "red",pch = 16)
  points(M[which(M< -m.cut  & A<p.cut)], -log10(p[which(M< -m.cut  & A<p.cut)]), col = "green", pch = 16)

  #abline(v = c(m.cut, -m.cut), col = "black", lty = 2)
  abline(h = -p.transform(max(p[A<p.cut])), col = "black", lty = 2)
  text(x=5, y=-p.transform(0.030), "P value: 0.05", col="blue")
  legend("topleft",c("Up-regulated genes", "Down-regulated genes"), col=c("red","green"), pch=16)
}

getRefseqID=function(x, pattern="_chr", part=1){x.split=unlist(strsplit(x,split=pattern)); id=x.split[seq(part,length(x.split),2)]; return(id)}

# cell lines
cell.lines=unique(matrix(unlist(str_split(dir(),"_")), ncol=4, byrow=TRUE)[,2])
ERminus=c("BT20", "MDAMB468", "MDAMB231")
ERplus=c("MCF7","T47D","BT474","ZR751")
Control="MCF10A"
# read data

rna.count=as.data.frame(lapply(dir()[grep("Raw",dir())],tryRead))
rownames(rna.count)=unlist(tryRead(dir()[grep("Raw",dir())][1],1))
rna=data.frame(control=rna.count[[Control]], ERplus=rna.count[,ERplus], ERminus=rna.count[,ERminus])
rna=rowFilter(rna,1)

millionsMapped <- colSums(rna)/1e+06
rpm <- rna/millionsMapped
# changes from normal
log2FC = log2(rpm[,-1]/rpm[,1])


# rpkm
# 1. get annotation for exon
#hg18ref=loadFeatures("~/work/project/RRBS/preliminary/GSE27003_RAW/refhg18.sqlite")
#exonRanges <- exonsBy(hg18ref, "tx", use.name=TRUE)
load("/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/hg18ref.Rdata")
exonRanges = hg18ref.exon
numBases <- sum(width(exonRanges))
#names(numBases) <- names(exonRanges)
names(numBases) <- getRefseqID(names(exonRanges))
geneLengthsInKB <- numBases/1000

# get genename gene id table
#hg18reftable = normalRead("~/work/project/RRBS/preliminary/GSE27003_RAW/refseq_hg18_113011.txt", header=FALSE)
hg18reftable = normalRead('/scratchLocal01/shl2018/eRRBS/bcdata/GSE/refseq_hg18_113011.txt', header=FALSE)
namesToId = hg18reftable[,1]
names(namesToId) = str_c(hg18reftable[,3], substring(hg18reftable[,2], 4), sep=".")
idToNames = str_c(hg18reftable[,3], substring(hg18reftable[,2], 4), sep=".")
names(idToNames) = hg18reftable[,1]

width = aggregate(numBases, by=list(geneName = idToNames[names(numBases)]), max)
width2 = width[,2]
names(width2) = width[,1]
width3=width2[rownames(rna)]
# filter out rna that donot have a matched gene id 
counts = rna[which(!is.na(width3)),]
millionsMapped <- colSums(rna)/1e+06
rpm <- rna/millionsMapped
rpkm <- rpm/geneLengthsInKB
#rna=rowNoZero(rna)
colnames(rna)=c("Ctrl",str_c("plus",1:4), str_c("minus",1:3))
cond_plus=2:5
cond_minus=6:8
medianNormal=function(count){
  # Compute the geometric mean of the gene counts (rows in lncap_counts) across all samples in the experiment as a pseudo-reference sample.
  geoMeans=as.matrix(exp(apply(log(count),1,mean)))
  # Each library size parameter is computed as the median of the ratio of the sample counts to those of the pseudo-reference sample.
  ratios = rna/geoMeans
  sizeFactors = apply(as.matrix(ratios),2,median)
  # The counts can be transformed to a common scale using size factor adjustment.
  base_count=rna/sizeFactors
  return(base_count)
}

meanNorm=function(count){
  ratios = rna/geoMeans
  sizeFactors = mean(count)
  # The counts can be transformed to a common scale using size factor adjustment.
  base_count=rna/sizeFactors
  return(base_count)
}
# Estimating Negative Binomial Distribution Parameters
estimateBaseParams = function(counts,sizeFactors, cond){return(data.frame(mean=log10(apply(counts[,cond]/sizeFactors[cond],1,mean)),
                                                               var=log10(apply(counts[,cond]/sizeFactors[cond],1,var))))}
plus_base_para = estimateBaseParams(rna,sizeFactors,cond_plus)
minus_base_para = estimateBaseParams(rna,sizeFactors,cond_minus)
par(mfrow=c(1,2))
boxplot(log2(rna))
boxplot(log2(base_count))

i=1
count=as.matrix(base_count[,2:8])
group=as.factor(c(rep('p',4), rep('m',3)))
m1=glm.nb(formula = count[i,] ~ group, init.theta = 1.03271315570559, link = log)
m2=glm(count[i,] ~ group, family="poisson")

topStdRna=topStd(rna[,-1],500)
#heatmap(topStdRna)

#cluster(rna[,-1], k=2)

# edgeR

x = rna[,-1]
adjust.method="fdr"
d = DGEList(x,group=group, genes=rownames(x))
d = calcNormFactors(d)
#plotMDS(d,main="MDS Plot for Breast Cancer Cells", labels=colnames(x))
d <- estimateCommonDisp(d)
print("common dispersion:");print(d$common.dispersion)
de.com <- exactTest(d)
print("Number of genes pass the multiple hypothesis testing:"); print(mhtGeneNum(de.com,adjust.method))
print(topTags(de.com,adjust.method,n=10))

cpm.d=cpm(d,TRUE)
cpm.sig=cpm.d[mhtGene(de.com,"fdr"),]

# filter out rna that do not have a matched gene id 
width3=width2[rownames(cpm.sig)]

cpm.sig.match = cpm.sig[which(!is.na(width3)),]
geneLengthsInKB <- width3[which(!is.na(width3))]/1000
rpkm <- cpm.sig.match/geneLengthsInKB

# load('/scratchLocal01/shl2018/eRRBS/bcdata/GSE/DEGrpkm.Rdata')

# get rpkm for meth feature matrix
gene.names = unlist(normalRead(x="meth.CE.symble.txt",header=FALSE))
rpkmForMeth=rpkm[gene.names,]
write.table(rpkmForMeth, file="rpkmForMeth.txt", quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#eff.libsize <- d$samples$lib.size * d$samples$norm.factor
#norm.count = t(apply(d$counts,1,function(x){return(x/eff.libsize)}))
# volcano plot
getPM=function(topTags.object,adjust.method="fdr"){return(topTags(de.com,adjust.method,n=length(topTags.object$genes[[1]]))[[1]][c("logFC","P.Value","adj.P.Val")])}
vocalnoPlot(getPM(de.com))


# Methylation
# read and formating
meth.raw.count=lapply(dir()[grep("CpGislands",dir())],tryRead6)
# match the island for each sample 
all.islands = unlist(lapply(dir()[grep("CpGislands",dir())], tryRead1))
islands = unique(all.islands)[tabulate(as.factor(all.islands))==length(grep("CpGislands",dir()))]
meth.count=foreach(i=1:length(meth.raw.count)) %do% meth.raw.count[[i]][islands]
# get meth matrix
meth = c(); for (i in 1:length(grep("CpGislands",dir()))){meth = cbind(meth,meth.count[[i]])}
colnames(meth)=unlist(str_split(dir()[grep("CpGislands",dir())], pattern="_"))[seq(2,32,4)]
meth=meth[,c(Control,ERplus,ERminus)]
colnames(meth)=c("Ctrl",str_c("plus",1:4), str_c("minus",1:3))

# function
is.identical=function(x){return(t(x-mean(x)) %*% (x-mean(x)) == 0)}
row.is.identical = function(x){return(unlist(apply(x,1,is.identical)))}
row.shapiro = function(x){return(unlist(apply(x,1,function(x){return(shapiro.test(x)$p.value)})))}
meth.rm.na=meth[apply(meth,1,function(x)all(!is.na(x))),]
meth.filtered=meth.rm.na[-which(row.is.identical(meth.rm.na)),]

# normalization
sizeFactors=apply(meth.rm.na,2,median)/median(meth.rm.na)

# var > mean? yes for 133148, no for 27981 + 242 (==)
x = estimateBaseParams(meth.rm.na,sizeFactors,cond=2:8)
length(which(abs(x$mean)<abs(x$var))) 

# edgeR assume the distribution is the 
meth.com=edgeRfun(meth.rm.na[,2:8])

# linear regression
# feature matrix
featureM = t(rna[mhtGene(de.com),2:8])
meth.list = normalRead("MethList.txt", header=FALSE)
gene.list=normalRead("gene.chr.txt", header=FALSE)

overlap.gene=which(unlist(gene.list) %in% rownames(rna))
overlap.meth=which(unlist(meth.list) %in% rownames(meth.rm.na))
overlap = which(overlap.gene %in% overlap.meth)
  #unique(str_c(EM_list[,1], ".",gsub(EM_list[1:10,3], pattern="chr", replacement="")))
meth_gene = str_c(rownames(rna)[])
featureM = cbind(t(rna[overlap,2:8]), t(meth.rm.na[overlap,2:8]))
status = as.factor(c(rep(1,4),rep(-1,3)))
# glmnet
selectFeature=function(featureM, status, maxit=5000){
  cv.fit <- cv.glmnet(featureM, status, family = "binomial", maxit = maxit, grouped=FALSE)
  fit <- glmnet(featureM, status, family = "binomial", maxit = maxit)
  Coefficients <- coef(fit, s = cv.fit$lambda.min)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  selected_features=colnames(featureM)[Active.Index]
  return(selected_features)
}

glm.bootstrap=function(x){selectFeature(featureM[resamples.all[x,],], status[resamples.all[x,]])}
resamples <- function(data, size, times=20, replace=FALSE){tmp=lapply(1:times, function(i)sample(data, replace = replace,size=size)); return(matrix(unlist(tmp), ncol=size, byrow=TRUE))}
re.plus=resamples(1:4,3,times=100)
re.minus=resamples(5:7,2,times=100)
resamples.all=cbind(re.plus,re.minus)

bootstrapGlmNet=function(resamples.all){
  x=lapply(1:nrow(resamples.all), glm.bootstrap)
  showup=tabulate(as.factor(unlist(x)))
  names(showup)=unique(unlist(x))
  return(sort(showup, decreasing=TRUE))
}

bootstrap=bootstrapGlmNet(resamples.all)
# use the diff expr and meth together table in the 
# 
# AHCTF1
# ABLIM1
# ALOXE3
# ANAPC2
# AGK
# ABHD4

# svm
library("kernlab")
status = as.factor(c(rep("plus",4),rep("minus",3)))
data=cbind(t(rna[mhtGene(de.com),2:8]), status)
model <- ksvm(status ~ ., data = data, type = "C-bsvc", kernel = "rbfdot", kpar = list(sigma = 0.1), C = 10,prob.model = TRUE)
# svm LiblineaR
