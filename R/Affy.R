library(affy)
library("limma")
targets=dir()
abatch=ReadAffy(filenames=targets)
eset=rma(abatch)
e=exprs(eset)

# paired data design matrix
pid=factor(rep(1:41,2))
Treat=factor(rep(c("D","R"),each=41))
design=model.matrix(~pid+Treat)
fit=lmFit(eset, design)
fit=eBayes(fit)

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
