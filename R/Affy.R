library(affy)
library("limma")
targets=dir()
abatch=ReadAffy(filenames=targets)
eset=rma(abatch)
e=exprs(eset)

# manu code
files <- dir ('.', pattern='CEL')
cel.data <- ReadAffy (filenames=files) # Loads all data into cel.data
# Background correction and normalization
# expresso function
normalized.data <- expresso(cel.data, bgcorrect.method = "rma", normalize.method = "quantiles", 
        pmcorrect.method = "pmonly", summary.method = "medianpolish")
exprs.data <- exprs (normalized.data) # Matrix : dimensions - no.probesets * samples
# ExpressionSet object: Column names will be filenames and rownames will be probeset names
# Determine gene names for probesetset names
# annaffy package for determining gene names: aafSymbol
library (annaffy)
probesets <- rownames (exprs.data)
chip <- 'hgu133plus2.db'
# Get gene names for all probesets
symbols <- sapply (probesets, function (x) { y<-aafSymbol (x, chip)[[1]]; ifelse (length (y!=0), y , NA) } )
# Meidan summarize to get a gene by gene basis
agg.data <- aggregate (exprs.data, by=list (symbols), median)
# 1st column represents gene names

genes.exprs <- as.matrix (agg.data[,-1])
rownames (genes.exprs) <- agg.data[,1]

save (genes.exprs, file='genes_exprs.Rdata')

# paired data design matrix
pid=factor(rep(1:41,2))
Treat=factor(rep(c("D","R"),each=41))
design=model.matrix(~pid+Treat)
fit=lmFit(genes.exprs, design)
fit=eBayes(fit)
# Create contrasts: The set of comparisons required
contrasts <- c('TreatR')
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
# Linear fit to compare the two sets
fit2 <- contrasts.fit(fit, contrast.matrix)
# Calculates the p-values
fit2 <- eBayes(fit2)

# for feature matrix
top2k=topTable(fit2, number=2000)
sig.gene=top2k$ID[which(top2k$P.Value<0.05)]
sig.exprs=genes.exprs[sig.gene,]
x.exprs=t(log((sig.exprs[,42:82]+1)/(sig.exprs[,1:41]+1)))


library(annotate)
biocLite("hgu133plus2cdf")
biocLite("hgu133plus2.db")
library("hgu133plus2cdf")
library("hgu133plus2.db")

ID=featureNames(eset)
Symbol <- getSYMBOL(ID,"hgu133plus2.db")
fData(eset) <- data.frame(ID=ID,Symbol=Symbol)
fit$genes$NAME=Symbol
topTable(fit,coef=2,number=15,genelist=fit$genes$NAME)

# clinic
clinic=read.table("clinicalStats.txt",header=TRUE,sep="\t")
x.clinic=cbind(Type=clinic$Type,Gender=clinic$Gender,Relapse.type=clinic$Relapse.type)
y=clinic$Months.to.relapse
names(y)=1:length(y)
X=cbind(x.clinic,x.exprs)

# shuffle samples
y.sf=sample(y)
X.sf=X[as.integer(names(y.sf)),]

# binary
status=rep(1,length(y)); status[y<median(y)]=0
y.bi.sf=status[as.integer(names(y.sf))]
# glmnet
reg.fit.test=function(X.sf,y.sf,tn=32, family='gaussian',alpha=0.5,nfolds=8,pmax=30){
  m=length(y.sf)
  X=X.sf[1:tn,]; y=y.sf[1:tn]
  X.test=X.sf[(tn+1):m,]; y.test=y.sf[(tn+1):m]
  
  e.cv  <- cv.glmnet( X, y, nfolds=nfolds)
  e.l   <- e.cv$lambda.min
  e.fits <- glmnet( X, y, family=family, alpha=alpha, nlambda=100,pmax=pmax)
  Coefficients <- coef(e.fits, s = e.l)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  names(Active.Coefficients) <- rownames(Coefficients)[Active.Index]
  cor=cor.test( predict(e.fits, X, type="response", s=e.l), y)$estimate
  test.cor=cor.test( predict(e.fits, X.test, type="response", s=e.l), y.test)$estimate
  summary(e.fits)
  return(list(coef=Active.Coefficients, cor=cor, test.cor=test.cor, e.fits=e.fits))
}

reg.fit.test(X.sf,y.sf, alpha=0.2,pmax=40,tn=35,nfolds=5)

<<<<<<< HEAD

# one sample subtraction a time and see the changes in cor, test.cor and coef

X.rand=matrix(sample(X.sf), ncol=ncol(X.sf))
colnames(X.rand)=colnames(X.sf)
# initialize variable
cor.m=c();test.cor.m=c();
real.coef=list();rand.coef=list()
for (i in -1:-41) {
  real=reg.fit.test(X.sf[i,],y.sf[i], alpha=0.2,pmax=40,tn=35,nfolds=5)
  rand=reg.fit.test(X.rand[i,],y.sf[i], alpha=0.2,pmax=40,tn=35,nfolds=5)
  cor.m=rbind(cor.m,c(real$cor,rand$cor))
  test.cor.m=rbind(test.cor.m,c(real$test.cor,rand$test.cor))
  real.coef[-i]=real$coef
  rand.coef[-i]=rand$coef
}
  
=======
# bootstrapping
boot.huber = function(data=data, indices){
  data = data[indices,]
  mod = reg.fit.test(data[,-1],data[,1], alpha=0.2,pmax=40,tn=35,nfolds=5)
  mod$coef
}
  
data=cbind(X.sf, y.sf)



>>>>>>> d66c6770a75792117b7cb34771a5b963e08a6353
