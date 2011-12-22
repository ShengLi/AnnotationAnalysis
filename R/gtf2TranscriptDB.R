# GTF to transcriptDB
# The code is changed based on gist 986084 by DarioS.
# The difference: if exon_number is not assinged to gtf,
# The current script will calculate the exon_number to gtf.db

# awk '$1 !~ "_hap|_random"{print}' genes.gtf|grep 'NM_' > NMgenes.gtf 
filename = "NMgenes.gtf"

# from gist
require(rtracklayer)
require(GenomicRanges)
require(GenomicFeatures)

min.info <- c("gene_id", "transcript_id", "exon_number")

if (verbose) message("Importing ", filename)
gtf <- import.gff(filename, asRangedData=FALSE)

    #parse the per exon attributes
if (verbose) message("Parsing gene/transcript/exon ids.")
e.info <- strsplit(gsub("\"|;", "", values(gtf)$group), split=" ")

ex.info <- lapply(e.info, function(x) {
    data <- x[seq(2, length(x), 2)]
    names(data) <- x[seq(1, length(x), 2)]
    data
})


# end gist

gene.filter='NM'
if(any(names(ex.info[[1]]) %in% "exon_number")){
  exon.info=ex.info
} else {
  
  tx_exonid=c(); 
  groupx=ex.info
  for(i in 1:length(groupx)){
    
    txid=groupx[[i]][3]
    existing.txid=names(tx_exonid)
    
    if (!any(existing.txid %in% txid)){
      # for additional attributes to GRanges
      id = 1; 
      # for counting exon number
      txexon=id; 
      names(txexon)=txid; 
      tx_exonid=c(tx_exonid,txexon); 
    } else {
      id = 1 + tx_exonid[txid];
      
      tx_exonid[txid]=id
    }
  
    names(id)="exon_number";
    groupx[[i]]=c(groupx[[i]],id)

  }
  #groupx.txid=unlist(lapply(1:length(groupx), function(i){groupx[[i]]['transcript_id']}))
  #passed.idx=grep(gene.filter, groupx.txid)
  
  exon.info=groupx
}



attribs <- names(exon.info[[1]])
if (!all(min.info %in% attribs)) stop("Not all required attributes are in this GTF file.")
exon.info <- do.call(data.frame, c(lapply(attribs, function(x)
                                                   sapply(exon.info, "[", x)),
                                          stringsAsFactors = FALSE))
colnames(exon.info) <- attribs

values(gtf) <- exon.info

if (verbose) message("Creating tables.")
#make transcripts table
exons.by.tx <- split(gtf, values(gtf)$transcript_id)
transcripts <- data.frame(
    tx_id = 1:length(exons.by.tx),
    tx_name = names(exons.by.tx),
    tx_chrom = as.character(seqnames(unlist(exons.by.tx))[start(exons.by.tx@partitioning)]),
    tx_strand = as.character(strand(unlist(exons.by.tx))[start(exons.by.tx@partitioning)]),
    tx_start = IRanges::sapply(start(ranges(exons.by.tx)), min),
    tx_end = IRanges::sapply(end(ranges(exons.by.tx)), max),
    stringsAsFactors = FALSE)

#make exons table
exons.ord <- unlist(exons.by.tx)
splicings <- data.frame(
    tx_id = rep(1:length(exons.by.tx), elementLengths(exons.by.tx)),
    exon_rank = as.integer(values(exons.ord)$exon_number),
    exon_chrom = as.character(seqnames(exons.ord)),
    exon_strand = as.character(strand(exons.ord)),
    exon_start = start(exons.ord),
    exon_end = end(exons.ord),
    stringsAsFactors = FALSE)

#make genes table
gene.txs <- tapply(values(gtf)$transcript_id, values(gtf)$gene_id, unique)
genes <- data.frame(
    tx_name = unlist(gene.txs),
    gene_id = rep(names(gene.txs), sapply(gene.txs, length)),
    stringsAsFactors = FALSE)

#create the db
if (verbose) message("Creating TranscriptDb.")
gtf.db <- makeTranscriptDb(transcripts, splicings, genes)

save(gtf.db, file="/scratchLocal01/shl2018/genome/hg18/Homo_sapiens/UCSC/hg18/Annotation/Archives/archive-2011-08-30-21-53-59/Genes")




# read in gtf file directly
NMgenes.gtf=read.table("NMgenes.gtf", colClasses=c("character","character","character", "integer","integer", "character","character","character","character"), stringsAsFactors=FALSE,sep="\t",header=FALSE)
NMgenes.gtf$V9=gsub(";", "",NMgenes.gtf$V9)
stringTotxid=function(x){unlist(strsplit(x," "))[grep('NM_',unlist(strsplit(x," ")))]}
NMgenes.gtf$txid=stringTotxid(NMgenes.gtf$V9)

# function for manu function
.extract.field.from.gtf <- function (gtf, field) {
    vals <- sapply (strsplit (sapply (strsplit (gtf, field),'[[', 2), ';'), '[[', 1)
    return (vals)
}


.extract.field.from.gtf=function(meta, pattern="transcript_id"){
    x=as.character(meta)
    splitedx=unlist(strsplit(x,";"))
    gsub(" ","",gsub(pattern, "",splitedx[grep(pattern,splitedx)]))
}
# end function
# code from Manu
transcriptdb.from.gtf <- function (gtf.file, chrominfo=NULL, file=NULL) {

  show ('Reading gtf file...')
  gtf.data <- read.table (gtf.file, sep='\t')
  if (!is.null (chrominfo))
    gtf.data <- gtf.data[gtf.data[,1] %in% chrominfo$chrom,]

  ## Construct "splicing" data.frame
  show ('Building splicing object...')
  exons <- gtf.data[gtf.data[,3] == 'exon',]
  exons <- GRanges (exons[,1], IRanges (exons[,4], exons[,5]),
                    exons[,7], meta = exons[,9])
  cds <- gtf.data[gtf.data[,3] == 'CDS',]
  cds <- GRanges (cds[,1], IRanges (cds[,4], cds[,5]), cds[,7], meta = cds[,9])
  
  ## Transcript names
  tx.names <- .extract.field.from.gtf (elementMetadata (exons)$meta,'transcript_id ')
  unique.tx.names <- tx.names[!duplicated(tx.names)]
  tx.ids <- 1:length (unique.tx.names); names (tx.ids) <- unique.tx.names

  ## Exon ranks
  strand <- tapply (as.character (strand (exons)), tx.names, '[', 1)[unique.tx.names]
  temp <- sprintf ("%s.%d.%s", tx.names, 1:length (tx.names), strand[tx.names])
  exon.ranks <- tapply (temp, tx.names,
                        function (x) { strand <- strsplit (x[1], '.', fixed=TRUE)[[1]][3]
                                       y <- 1:length(x);
                                       if (strand == '-')
                                         y <- length (x):1
                                       names (y) <- x; y})
  exon.rank.names <- unlist (lapply (exon.ranks, names), use.names=FALSE)
  exon.ranks <- unlist (exon.ranks, use.names=FALSE)
  names (exon.ranks) <- exon.rank.names
  exon.ranks <- exon.ranks[temp]
  names (exon.ranks) <- NULL

  ## Build splicing object
  splice.tx.ids <- tx.ids[tx.names]; names (splice.tx.ids) <- NULL
  splicings <- data.frame (tx_id=splice.tx.ids,
                           exon_rank=exon.ranks,
                           exon_start=start (exons),
                           exon_end=end (exons),
                           cds_id=NA,
                           cds_start=NA,
                           cds_end=NA)

  ## Update cds
  overlaps <- as.matrix (findOverlaps (cds, exons))
  cds.tx.names <- .extract.field.from.gtf (elementMetadata (cds)$meta,
                                           'transcript_id ')
  overlaps <- overlaps[cds.tx.names[overlaps[,1]] == tx.names[overlaps[,2]],]
  cds.overlaps <- pintersect (cds[overlaps[,1]], exons[overlaps[,2]])
  splicings[overlaps[,2], c('cds_id', 'cds_start', 'cds_end')] <-
    cbind (1:nrow (overlaps), start (cds.overlaps), end (cds.overlaps))


  ## Build transcript object
  show ('Building transcript object....')
  chrom <- tapply (as.character (seqnames (exons)), tx.names, '[',
                   1)[unique.tx.names]
  strand <- tapply (as.character (strand (exons)), tx.names, '[',
                    1)[unique.tx.names]
  start <- tapply (start (exons), tx.names, min)[unique.tx.names]
  end <- tapply (end (exons), tx.names, max)[unique.tx.names]
  names (tx.ids) <- names (chrom) <- names (start) <- names (end) <-
    names (strand) <- NULL
  transcripts <- data.frame (tx_id = tx.ids,
                             tx_name = unique.tx.names,
                             tx_chrom = chrom,
                             tx_strand = strand,
                             tx_start = start,
                             tx_end = end)

  genes <- .extract.field.from.gtf (elementMetadata (exons)$meta, 'gene_name ')
  genes <- data.frame (tx_name=unique.tx.names, gene_id=tapply (genes,
                                                  tx.names, '[', 1)[unique.tx.names])

  ## Build txdb
  show ('Building transcript db...')
  txdb <- makeTranscriptDb(transcripts, splicings, genes, chrominfo=chrominfo)

  if (!is.null (file))
    saveFeatures (txdb, file=file)

  return (txdb)
}

# code end 
chrominfo=read.table("ChromInfo.txt", stringsAsFactors=FALSE, header=FALSE)
names(chrominfo)=c("chrom","length")
hg18.txdb=transcriptdb.from.gtf('NMgenes.gtf', chrominfo=chrominfo, file='~/annotation/hg18/refseq/hg18txdb.Rdata')

# get regional ranges
txdb=hg18.txdb
tx.db = transcripts(txdb)
exon.db=exonsBy(txdb,"tx", use.name=TRUE)
intron.db=intronsByTranscript(txdb,use.name=TRUE)
cds.db=cdsBy(txdb,"tx", use.name=TRUE)
threeutr.db = threeUTRsByTranscript(txdb, use.name=TRUE)
fiveutr.db = fiveUTRsByTranscript(txdb, use.name=TRUE)

# first
specifyElement=function(db,pos='first'){
  # assume that the each element in the GRangesList has unique name
  require(data.table)
  if(pos=='first') strands="+"
  else if(pos=='last') strands="-"

  df=as.data.frame(db)
  dt=data.table(df)
  dt1=dt[,list(seqnames=seqnames[1],
         start=ifelse(strand[1]==strands,start[1],start[length(start)]),
         end=ifelse(strand[1]==strands,end[1],end[length(end)]),
         strand=strand[1]) ,by=element]
  first.db=GRanges(seqnames=dt1$seqname, 
                   ranges=IRanges(
                    start=dt1$start,
                    end=dt1$end),
                   strand=dt1$strand,
                   name=dt1$element)
  first.db
}
firstin.db=specifyElement(intron.db,'first')
lastin.db=specifyElement(intron.db,'last')
firstex.db=specifyElement(exon.db,'first')
lastex.db=specifyElement(exon.db,'last')

# basic txdb
txStart=start(tx.db)
txEnd =end(tx.db)
txStrand=as.character(strand(tx.db))
txChr=as.character(seqnames(tx.db))
txName = elementMetadata(tx.db)[,"tx_name"]

# prmoter
extendTss=function(txdb, upFlank, dnFlank){
  txStart=start(tx.db)
  txEnd =end(tx.db)
  txStrand=as.character(strand(tx.db))
  txChr=as.character(seqnames(tx.db))
  txName = elementMetadata(tx.db)[,"tx_name"]
  txTSS = txStart
  txTSS[txStrand=="-"] = txEnd[txStrand=="-"]
  
  extended=GRanges(seqnames=txChr,
                   ranges=IRanges(start=txTSS-upFlank, end=txTSS+dnFlank),
                   strand=txStrand,
                   name = txName)
  extended  
}

tssPM1kb.db=extendTss(txdb,1000,1000)
tssM2kb.db=extendTss(txdb,2000,0)
tssM3kb.db=extendTss(txdb,3000,0)
# tss
tss.db=extendTss(txdb,0,0)

# cpgi by tx-GR
# need to consider strandness?
# need to set a length limitation?
CpG.file="~/software/seqframework/eRRBS_pipeline/annotation_data/hg18/hg18.CpG.UCSCtable.txt"
chrom.info="~/software/seqframework/eRRBS_pipeline/annotation_data/hg18/hg18.chromInfo.txt"
source("/home/shl2018/software/seqframework/eRRBS_pipeline/scripts/CpGmeth.functions.R")
source('/scratchLocal01/shl2018/eRRBS/rpackage/methylkit/R/annotate.R')
g.cpg=read.CpG.UCSC.table(GpG.file,chrom.info)
cpgTSSDistTable=.nearest.2bed(tss.db, g.cpg)

cpgi.db = GRanges(seqnames=cpgTSSDistTable$seqnames.y,
  				ranges = IRanges(
						start= cpgTSSDistTable$start.y, 
						end = cpgTSSDistTable$end.y),
					strand = cpgTSSDistTable$strand.y,
					name= cpgTSSDistTable$name,
					dist = cpgTSSDistTable$dist,
					length= cpgTSSDistTable$length.y,
					cpgNum= cpgTSSDistTable$cpgNum.y,
					gcNum = cpgTSSDistTable$gcNum.y,
					perCpg= cpgTSSDistTable$perCpg.y,
					perGc = cpgTSSDistTable$perGc.y,
					obsExp = cpgTSSDistTable$obsExp.y)

shores=c( flank(cpgi.db,1000),flank(cpgi.db,1000,FALSE) )
shore.db = split(shores, rep(1:length(cpgi.db),2))
names(shore.db)=as.character(values(cpgi.db)[,1])


hg18ref=list(exon=exon.db, intron=intron.db,
             tssPM1kb=tssPM1kb.db, tssM2kb=tssM2kb.db, tssM3kb=tssM3kb.db,
             cpgi=cpgi.db, shore=shore.db, 
             cds=cds.db,threeutr=threeutr.db, fiveutr=fiveutr.db,
             firstin=firstin.db, lastin=lastin.db,
             firstex=firstex.db, lastex=lastex.db)
save(hg18ref,file="~/annotation/hg18/refseq/hg18regions.Rdata")


# get the gene and symbol
exon.width=width(hg18ref$exon)
tx.width=unlist(lapply(1:length(exon.width), function(i){sum(exon.width[[i]])}))
tx.name=names(hg18ref$exon)
names(tx.width)=tx.name
names(txChr)=elementMetadata(tx.db)[,'tx_name']
tx.chr=txChr[tx.name]
dupname=str_c('NM_',sapply(strsplit(tx.name,'_'),'[[',2))
names(dupname)=tx.name

# from refseqid to symbol
library(sequencingUtils)
library(utilities)
include.eg.db('hg18')
symbols=refseq.to.symbol(dupname)

# sybmol + chr + width
CEsymbols=str_c(symbols,sapply(strsplit(tx.chr,'chr'),'[[',2),sep='.')
names(CEsymbols)=tx.name
CEsymbolswidth=str_c(CEsymbols,tx.width,sep="_")
names(CEsymbolswidth)=tx.name
# get the symbol+chr with the max width
max.width=aggregate(tx.width, by=list(CEsymbols), max)
symbols.max.width=str_c(max.width[,1], unlist(max.width[,2] ), sep="_")
matched.tx.name=tx.name[which(CEsymbolswidth %in% symbols.max.width)]
matched.CEsymbols=CEsymbols[which(CEsymbolswidth %in% symbols.max.width)]
names(matched.CEsymbols)=tx.name[which(CEsymbolswidth %in% symbols.max.width)]

# get matched meid and exprid
CEmeti1=matched.CEsymbols[which(matched.CEsymbols %in% CE1.symbol)]
# remove duplicates
CEmeti1=CEmeti1[-which(CEmeti1%in%CEmeti1[which(duplicated(CEmeti1))])]

CEmeti2=matched.CEsymbols[which(matched.CEsymbols %in% CE2.symbol)]
# remove duplicates
CEmeti2=CEmeti2[-which(CEmeti2%in%CEmeti2[which(duplicated(CEmeti2))])]

load('/scratchLocal01/shl2018/eRRBS/bcdata/GSE/DEGrpkm.Rdata')
colnames(rpkm)=sapply(strsplit(colnames(rpkm),'\\.'),'[[',2)

CEmeti1=CEmeti1[which(CEmeti1%in%rownames(rpkm) & names(CEmeti1 )%in%rownames(reMethRatioMap[[1]]))]
CEmeti2=CEmeti2[which(CEmeti2%in%rownames(rpkm) & names(CEmeti2)%in%rownames(reMethRatioMap[[2]]))]

# get rpkm
y1=rpkm[CEmeti1,]
y2=rpkm[CEmeti2,]

# get feature matrix
X1=foreach(i=1:length(reMethRatioMap)) %do% reMethRatioMap[[i]][names(CEmeti1),]
X2=foreach(i=1:length(reMethRatioMap)) %do% reMethRatioMap[[i]][names(CEmeti2),]

X1=foreach(i=1:length(regionsMethRatio)) %do% regionsMethRatio[[i]][names(CEmeti1),]
X2=foreach(i=1:length(regionsMethRatio)) %do% regionsMethRatio[[i]][names(CEmeti2),]

save(X1,X2,y1,y2, file='/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/glmnetXy.Rdata')

coefMatrix=function(X,y, cnames){
  results=c()
  for(i in 1:7) {
    results=cbind(results, selectFeature(X[[i]], y[,i],'gaussian',20000)); 
  }
  colnames(results)=cnames
  results
  #[which(rowSums(results)!=0),]
}

rcoefMatrix=function(X,y, cnames){
  results=c()
  for(i in 1:7) {
    results=cbind(results, selectFeature(X[[i]], sample(y[,i]),'gaussian',20000)); 
  }
  colnames(results)=cnames
  results
}


coefPlot=function(coef.matrix){
  a <- coef.matrix
  #rowv <- as.dendrogram(hclust(as.dist(1-cor(t(a)))))
  colv <- as.dendrogram(hclust(as.dist(1-cor(a))))
  heatmap.2(a, scale="row", Colv=colv)
}

heatmap.wrapper=function(x, rm.zero=TRUE,cexRow=1.2, cexCol=1.5, margins=c(6,6)){
  if(rm.zero) x=x[which(rowSums(x)!=0),]
  heatmap.2(x, col=topo.colors(75), scale="row", ColSideColors=patientcolors,margins=margins,
          key=FALSE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow,cexCol=cexCol, dendrogram="col")
}
colNames=c('plus1', 'plus2', 'plus3', 'plus4', 'minus1', 'minus2', 'minus3')
coef1=coefMatrix(X1,y1,colNames);heatmap.wrapper(coef1[-1,], TRUE)
coef2=coefMatrix(X2,y2,colNames);heatmap.wrapper(coef2[-1,], TRUE)

rcoef1=rcoefMatrix(X1,y1,colNames); heatmap.wrapper(rcoef1[-1,], TRUE)
rcoef2=rcoefMatrix(X2,y2,colNames); heatmap.wrapper(rcoef2[-1,], TRUE)

patientcolors=rep(c("#E890E5","#05F0FC"),c(4,3))

cexRow=1.2; cexCol=1.5; margins=c(6,6)
heatmap.2(coef1[-1,], col=topo.colors(75), scale="row", ColSideColors=patientcolors,margins=margins,
          key=FALSE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow,cexCol=cexCol, dendrogram="col")
heatmap.2(coef2[c(-1,-5),], col=topo.colors(75), scale="row", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow,cexCol=cexRow, dendrogram="col")

heatmap.2(rcoef1[-1,], col=topo.colors(75), scale="row", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow,cexCol=cexRow, dendrogram="col")
heatmap.2(rcoef2[-1,], col=topo.colors(75), scale="row", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow,cexCol=cexRow, dendrogram="col")
# predict
i=1;X=X1[[i]];y=y1[,i];cv.fit = cv.glmnet(X,y,family='gaussian',maxit=20000, grouped=FALSE);ypre=predict(cv.fit,X); cor(t(y),t(ypre))
