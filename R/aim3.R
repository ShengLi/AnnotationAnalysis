# for aim3,
# divide samples into 2 groups
# based on time to R.

# preliminary example
# 10 samples; 1-6: se, 7:10, pe

### y vector generation:

## divide samples into 2 groups by 
## the median of time to R
timeToR=c(3.8,4.3, 3.1,2.6,3.2,1.5,2.1,1.0,3.6,0.8)
names(timeToR)=1:10
cutoff=median(timeToR)

slow.idx=which(timeToR>cutoff)
fast.idx=which(timeToR<cutoff)

status=rep(1,10); status[slow.idx]=0 

### x feature matrix
## cli feature matrix part
gender=c(1,1,-1,1,-1,-1,1,1,-1,1)
race=c(1,1,1,1,-1,0,1,0,1,1)
age=c(16,15.8,14.3,6,17,7.3,1.9,18,13,16)

x.c=cbind(gender,race,age)

## expression feature matrix part

# cuffdiff output 
# source: calculate
require(doMC)
l=10
peall.files=paste(dir("/scratchLocal01/shl2018/peall/data/cufflinks/", pattern="peall",full.names=TRUE),"/gene_exp.diff",sep="")
all.files=paste(dir("/scratchLocal01/shl2018/all/cufflinks", full.names=TRUE),"/gene_exp.diff",sep="")
files=c(all.files, peall.files)
diff.files=foreach(i=1:l)%do% read.table(files[i], stringsAsFactors=FALSE, header=TRUE)
diff.gene_id=foreach(i=1:l)%do% diff.files[[i]]$gene_id[which(diff.files[[i]]$significant=="yes")]
gene.diff=foreach(i=1:l)%do% {x=diff.files[[i]];ok.idx=which(x$status=="OK"); y=(x$value_2-x$value_1)[ok.idx];names(y)=x$gene[ok.idx];return(y)}
common.gene=names(which(table(names(unlist(gene.diff)))==10))
diff.m=c();for(i in 1:l) diff.m=cbind(diff.m, gene.diff[[i]][common.gene])
# apply normal distribution test to see how 
diff.m=diff.m[,-3]
y=apply(diff.m,1,function(x){shapiro.test(x)$p.value})
norm.diff=diff.m[which(y>=0.05),]
sd.diff=apply(norm.diff,1, sd)
se.diff=sd.diff/sqrt(l)
df.diff=l-1
mean.diff=apply(norm.diff,1,mean)
D=0
t=(mean.diff-D)/se.diff
p.value=2*pt(-abs(t),df=l-1)
# end cuffdiff

## apply wilcox.test to get the paired study
# get the D matrix and R matrix
gene.D=c(); gene.R=c();
for(i in 1:l){
  x=diff.files[[i]];
  common.idx=which(x$gene%in%common.gene); 
  D=x$value_1[common.idx];
  names(D)=x$gene[common.idx];
  R=x$value_2[common.idx]
  names(R)=x$gene[common.idx];
  gene.D=cbind(gene.D,D);
  gene.R=cbind(gene.R,R);
  }
colnames(gene.D)=paste("D",1:l,sep="_")
colnames(gene.R)=paste("R",1:l,sep="_")

# because in D_3, the rpkm is all zero,so skip this first
gene.D=gene.D[,-3]
gene.R=gene.R[,-3]
y=apply(cbind(gene.D,gene.R),1,function(x){shapiro.test(x)$p.value})
# for normal run t.test
nidx=which(y>=0.05)
uidx=which(y<0.05)
p.t=sapply(1:length(nidx), function(i){t.test(gene.D[nidx,][i,],gene.R[nidx,][i,],paired=TRUE)}); padj.t=p.adjust(p.t[3,],'fdr')
p.w=sapply(1:length(uidx), function(i){wilcox.test(gene.D[uidx,][i,],gene.R[uidx,][i,],paired=TRUE)}); padj.w=p.adjust(p.w[3,],'fdr')
# wilcox for all data
p.wil=sapply(1:nrow(gene.D), function(i){wilcox.test(gene.D[i,],gene.R[i,],paired=TRUE)}); padj=p.adjust(p.wil[3,],'fdr') 
# cufflinks output
wil.idx=which(p.wil[3,]<0.05)
expr.logRatio=t(log((gene.R[wil.idx,]+1)/(gene.D[wil.idx,]+1)))
expr.diff=t(gene.R[wil.idx,]-gene.D[wil.idx,])

x.expr=expr.logRatio
# end cufflinks output

## feature matrix
#X=t(gene.R-gene.D)
X = cbind(x.c[-3,],expr.logRatio)
#y = status[-3]
y=timeToR[-3]

reg.fit=function(X,y,family='binomial',alpha=0.5,nfolds=3,pmax=10){
  e.cv  <- cv.glmnet( X, y, nfolds=nfolds)
  e.l   <- e.cv$lambda.min
  e.fits <- glmnet( X, y, family=family, alpha=alpha, nlambda=100,pmax=pmax)
  Coefficients <- coef(e.fits, s = e.l)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  names(Active.Coefficients) <- rownames(Coefficients)[Active.Index]
  cor=cor.test( predict(e.fits, X, type="response", s=e.l), y)$estimate
  return(list(coef=Active.Coefficients, cor=cor))
}

 

