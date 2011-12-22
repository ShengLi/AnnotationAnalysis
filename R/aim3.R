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

status=rep(1,10); status[fast.idx]=-1

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
norm.diff=diff.m[which(y>=0.05),]
sd.diff=apply(norm.diff,1, sd)
se.diff=sd.diff/sqrt(l)
df.diff=l-1
mean.diff=apply(norm.diff,1,mean)
D=0
t=(mean.diff-D)/se.diff
p.value=2*pt(-abs(t),df=l-1)
# end cuffdiff

# cufflinks output

# end cufflinks output