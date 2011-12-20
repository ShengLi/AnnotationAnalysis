## Methylation gene-based analysis
library(methylKit)
#
# function
#

# convert a vector to list object
vectorToList=function(x){x.list=list(); for(i in 1:length(x)) x.list[[i]]=x[i]; return(x.list)}

# only let protein coding refseq id pass
filterID=function(GRL,pattern="NM", is.grl=TRUE){if(is.grl)id=names(GRL) else id = elementMetadata(GRL)$name; idx = grep(pattern,id); return(GRL[idx])}

# given refseqid_chr_start_end, return refseqid
getRefseqID=function(x, pattern="_chr", part=1){x.split=unlist(strsplit(x,split=pattern)); id=x.split[seq(part,length(x.split),2)]; return(id)}

# given myobj and GRanges or GRangesList, return a regioal methylation level and unite different samples together. Wrapper for methylKit
methRegion=function(myobj,regions){myobj.region=regionCounts(myobj,regions); meth.region=unite(myobj.region); return(meth.region)}

# given symbol.chrNum, return the matched txUniID with the longest transcript. 
getCEid=function(CE1.symbol, refseq_sc,txUniId,txWidthNM){
  tw1=txWidthNM[which(refseq_sc %in% CE1.symbol)]
  ts1=refseq_sc[which(refseq_sc %in% CE1.symbol)]
  ti1=txUniId[which(refseq_sc %in% CE1.symbol)]
  names(ti1) = str_c(ts1,tw1,sep="_")
  twis1=aggregate(tw1,by=list(ts1),max)
  ceti1=ti1[str_c(twis1[,1], twis1[,2],sep="_")]
  return(ceti1)
}

# get methylation ratio
getRatioByRegion=function(meth.regions, id){
  midx=which(meth.regions$id %in% id)
  ratio=meth.regions[midx,meth.regions@numCs.index]/meth.regions[midx,meth.regions@coverage.index]
  colnames(ratio)=meth.regions@sample.ids
  rownames(ratio)=meth.regions$id[midx]
  return(ratio)
}

# given meth.resions and id return methylation per sample, include a talbe of every regions methylation
getRatioBySample=function(meth.regions.list, id){
  # meth.regions.list=list(exon=meth.exon, intron=meth.intron, promoter=meth.promoter, cpg=meth.cpg, shore=meth.shore)
  # id = CEmeti1
  sample.ids=meth.regions.list[[1]]@sample.ids
  region.ids=names(meth.regions.list)
  # exam region.ids
  if(region.ids[[1]] == "NULL") region.ids=1:length(meth.regions.list) else cat(paste( region.ids,"\n"))
  cat("loop through regions to get regional methylation ratio... \n")
  regions.ratio=list()
  samples.ratio=list()
  for(r in region.ids){
    regions.ratio[[r]] = getRatioByRegion(meth.regions.list[[r]], id)
    for(s in sample.ids){
      samples.ratio[[s]][[r]]= regions.ratio[[r]][[s]]
    }
  }
  samples.ratio.dt=list()
  for(s in sample.ids){
    samples.matrix=matrix(unlist(samples.ratio[[1]]), ncol=length(region.ids))
    rownames(samples.matrix)=id
    colnames(samples.matrix)=region.ids
    samples.ratio.dt[[s]]=as.data.frame(samples.matrix)
    }
  return(samples.ratio.dt)
}
# given refseq ID, a list of myDiff result as well as the regions names according to the sequence in the myDiffList
# return a summary table of myDiff for each refseq ID
# Result format could be byRefseqID of a list or a data.frame contain everything
sumMyDiff=function(refseq.id,myDiff.list, region.names, byRefseqId=TRUE){
#myDiff.list=list(myDiff.exon,myDiff.intron,myDiff.promoter, myDiff.cpg, myDiff.shore)
  # region.names = c("exon","intron","promoter","cpg","shore")
  if(!byRefseqId){
    myDiffSum=data.frame()
    for (i in 1:length(myDiff.list)){
      region.id.myDiff = list()
      for (k in 1:length(refseq.id)){
        myIdDiff=myDiff.list[[i]][grep(refseq.id[k],myDiff.list[[i]]$id),]
        region.id.myDiff = rbind(region.id.myDiff, myIdDiff)
        }
      region.id.myDiff = data.frame(region.id.myDiff, type=rep(region.names[i],nrow(region.id.myDiff)))
      myDiffSum=rbind(myDiffSum, region.id.myDiff)
      }
    }
  else {
    myDiffSum=list()
    for( k in 1:length(refseq.id)){
      region.id.myDiff = data.frame()
      for( i in 1:length(myDiff.list)){
        myIdDiff = myDiff.list[[i]][grep(refseq.id[k],myDiff.list[[i]]$id),]
        myIdDiff = data.frame(myIdDiff, type=rep(region.names[i],nrow(myIdDiff)))
        region.id.myDiff = rbind(region.id.myDiff, myIdDiff)
      }
      myDiffSum[k]=list(region.id.myDiff)
    }
    names(myDiffSum) = refseq.id
  }
  return(myDiffSum)
}

# given regional myobj, return the methylation ratio and match with meid 
regionRatio=function(myobj.region, meid){
  x=myobj.region; y=x$numCs/x$coverage; names(y)=x$id; 
  meth=y[meid]; meth[which(is.na(meth))]=0; return(meth)
}

matchMethRatio=function(txdblist,myobj.sub){
  require(doMC)
  registerDoMC(cores=8)

  t=length(txdblist)
  myobj.txdb=foreach(i = 1:t) %dopar% regionCounts(myobj.sub, txdblist[[i]])
  meid=unique(as.character(unlist(foreach(i=1:t) %dopar% myobj.txdb[[i]]$id)))
  meth=foreach(i=1:t) %dopar% regionRatio(myobj.txdb[[i]], meid)
  meth.m=matrix(unlist(meth), ncol=length(meth),byrow=FALSE)
  colnames(meth.m)=names(txdblist)
  rownames(meth.m)=meid
  return(meth.m)
}

# given a matrix return a binomial feature matrix
polyMap=function(x,k=2){
  require(doMC)
  cbn=combn(1:ncol(x),k)
  numcb=ncol(cbn)
  bn=foreach(i=1:numcb) %do% (x[,cbn[1,i]] * x[,cbn[2,i]])
  bn.m=matrix(unlist(bn),ncol=numcb,byrow=FALSE)
  cn=unlist(foreach(i=1:numcb) %do% {cn=colnames(x)[cbn[1:2,i]]; paste(cn[1],cn[2],sep="*")})
  colnames(bn.m)=cn
  rownames(bn.m)=rownames(x)
  bn.m
}
#
# end of function
#
load("~/annotation/hg18/refseq/hg18regions.Rdata")
# cell line my_CpG file arrangement.
ERplus=c("MCF7","T47D","BT474","ZR751")
ERminus=c("BT20", "MDAMB468", "MDAMB231")
nonmeta=c("MCF7","T47D","BT474","ZR751")
meta=c("MDAMB361", "MDAMB468", "MDAMB231","HS578T")
path='/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/'
pattern='_myCpG.txt'
files = str_c(path,c(ERplus,ERminus),pattern)
#files = str_c(path,c(nonmeta,meta),pattern)
file.list=vectorToList(files)
samples = c(ERplus,ERminus)
#samples = c(nonmeta,meta)
sample.id=vectorToList(samples)

# read meth myCpG file
myobj=read(file.list,sample.id=sample.id, assembly='hg18',treatment=c(1,1,1,1,0,0,0))
regionsMethRatio=foreach(i=1:length(myobj)) %do% matchMethRatio(hg18ref,myobj[[i]])
save(regionsMethRatio, file='/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/regionsMethRatio.Rdata')
save(myobj, file='/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/ERmyobj.Rdata')
reMethRatioMap=foreach(i=1:length(regionsMethRatio))%do% cbind(regionsMethRatio[[i]], polyMap(regionsMethRatio[[i]]))
# metastatic part:
nonmeta=c("MCF7","T47D","BT474","ZR751")
meta=c("MDAMB361", "MDAMB468", "MDAMB231","HS578T")
path='/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/'
pattern='_myCpG.txt'
files = str_c(path,c(nonmeta,meta),pattern)
file.list=vectorToList(files)
samples = c(nonmeta,meta)
sample.id=vectorToList(samples)

meta.myobj=read(file.list,sample.id=sample.id, assembly='hg18',treatment=c(0,0,0,0,1,1,1,1))
save(meta.myobj, file='/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/metamyobj.Rdata')



#
# conver the CpG di into gene based methylation

meth.exon=methRegion(myobj,hg18ref.exon)
meth.intron=methRegion(myobj,hg18ref.intron)
meth.promoter=methRegion(myobj,hg18ref.promoter)
meth.cpg=methRegion(myobj,hg18ref.cpg)
meth.shore=methRegion(myobj,hg18ref.shore)
# get differential profile:
myDiff.exon=calculateDiffMeth(meth.exon);         myDiff25p.exon=get.methylDiff(myDiff.exon)
myDiff.intron=calculateDiffMeth(meth.intron);     myDiff25p.intron=get.methylDiff(myDiff.intron)
myDiff.promoter=calculateDiffMeth(meth.promoter); myDiff25p.promoter=get.methylDiff(myDiff.promoter)
myDiff.cpg=calculateDiffMeth(meth.cpg);           myDiff25p.cpg=get.methylDiff(myDiff.cpg)
myDiff.shore=calculateDiffMeth(meth.shore);       myDiff25p.shore=get.methylDiff(myDiff.shore)

refseq.id=c("NM_004040","NM_006547","NM_032873")
refseq.to.symbol(refseq.id)             
myDiffSum=sumMyDiff(refseq.id,
                    myDiff.list=list(myDiff.exon,myDiff.intron,myDiff.promoter, myDiff.cpg, myDiff.shore),
                    region.names = c("exon","intron","promoter","cpg_island","cpg_island_shore"),
                    byRefseqId=TRUE)
# preparing for matrix
#refseq_id=getRefseqID(names(ref18.exon))
library(sequencingUtils)
library(utilities)
include.eg.db('hg18')
refseq_symbols=refseq.to.symbol(txNameNM)
refseq_chrs=getRefseqID(txChrNM,pattern="chr",2)
refseq_sc = str_c(refseq_symbols, refseq_chrs,sep=".")

# get the txUniId for the co-expressed genes. with this we can get the methylation level.
load("/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/CE.Rdata")
CEti1=getCEid(CE1.symbol,refseq_sc,txUniId,txWidthNM) # 306
CEti2=getCEid(CE2.symbol,refseq_sc,txUniId,txWidthNM) # 510

# get the CEti that have all five regions methylation
filterCEti=function(CEti1){
CEmeti1=CEti1[which(CEti1 %in% meth.exon$id & CEti1 %in% meth.intron$id & 
  CEti1 %in% meth.promoter$id & CEti1 %in% meth.cpg$id & CEti1 %in% meth.shore$id)]
return(CEmeti1)
}
CEmeti1 = filterCEti(CEti1)
CEmeti2 = filterCEti(CEti2)

meth.regions.list=list(exon=meth.exon, intron=meth.intron, promoter=meth.promoter, cpg=meth.cpg, shore=meth.shore)
meth.CE1=getRatioBySample(meth.regions.list, CEmeti1)
meth.CE2=getRatioBySample(meth.regions.list, CEmeti2)

listToBigTable=function(list){x=c(); for(i in 1:length(list)) x=rbind(x,as.matrix(list[[i]])); return(x)}
meth.CE1.table=listToBigTable(meth.CE1) # 385   5
meth.CE2.table=listToBigTable(meth.CE2) # 721   5
write.table(rbind(meth.CE1.table,meth.CE2.table), file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/meth.CE.table.txt",
            quote=FALSE, sep="\t", row.names=FALSE)

getFirstPart=function(x){unlist(strsplit(x,"_"))[seq(1,2*length(x),2)]}
write.table(getFirstPart(names(c(CEmeti1,CEmeti2))), file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/meth.CE.symble.txt", 
            quote=FALSE, sep="\t", col.names=FALSE,row.names=FALSE)

# get regional counts
counts=regionCounts(myobj,hg18ref.exon)
meth.ratio=lapply(1:length(myobj), function(i){
  counts.i=counts[[i]];
  meth.ratio=counts.i$numCs/counts.i$coverage; 
  names(meth.ratio)=counts.i$id;
  rest.id=txUniIdNM[-which(txUniIdNM %in% counts.i$id)]
  return(meth.ratio)})

#
# save
save(myobj, file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/myCpG.Rdata")
save(meth.exon,meth.intron,meth.promoter,meth.cpg,meth.shore,file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/meth.Rdata")
save(myDiff.exon,myDiff.intron,myDiff.promoter,myDiff.cpg,myDiff.shore,file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/myDiff.Rdata")
save(myDiff25p.exon,myDiff25p.intron,myDiff25p.promoter,myDiff25p.cpg,myDiff25p.shore,file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/myDiff25p.Rdata")
save(hg18ref.exon,hg18ref.intron,hg18ref.promoter,hg18ref.cpg,hg18ref.shore,g.cpg,hg18ref.TSS,txPos, file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/hg18ref.Rdata")
save(cpg.dist,file="/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/cpgDist.Rdata")
# cpg.dist from cpgTSSDistTable
