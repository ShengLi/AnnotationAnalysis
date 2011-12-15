source("/scratchLocal01/shl2018/eRRBS/rpackage/methylkit/R/annotate.R")
library(GenomicFeatures)
library(GenomicRanges)
# function
getUniID=function(x){
paste(seqnames(x[1,]), min(start(x)), max(end(x)),sep="_")
}


# end of function
hg18ref = loadFeatures("~/database/hg18/refhg18.sqlite")
hg18ref=loadFeatures("~/work/project/RRBS/preliminary/GSE27003_RAW/refhg18.sqlite")


# get the whole list of tx
hg18ref.tx = transcripts(hg18ref)
txChr=as.character(seqnames(hg18ref.tx))
txName = elementMetadata(hg18ref.tx)[,"tx_name"]
txChr=as.character(seqnames(hg18ref.tx))
txStart=start(hg18ref.tx)
txEnd =end(hg18ref.tx)
txStrand=as.character(strand(hg18ref.tx))
txWidth=width(hg18ref.tx) # for later determine use which tx
txSeqlengths=seqlengths(hg18ref.tx)
txUniId=str_c(txName,txChr,txStart,txEnd,sep="_")

## filter non-NM refseq id
txChrNM=txChr[grep("NM",txName)]
txUniIdNM=txUniId[grep("NM",txName)]
txWidthNM=txWidth[grep("NM",txName)]
txNameNM=txName[grep("NM",txName)]

library(utilities)
include.eg.db('hg18')
txSymbolNM=refseq.to.symbol(txNameNM)

# exon by tx-GRL
hg18ref.exon=exonsBy(hg18ref,"tx", use.name=TRUE)
hg18ref.exon
txPos= sapply(hg18ref.exon, getUniID) # may take a long time
exonsByName=names(hg18ref.exon)
exonsByUniId= str_c(exonsByName, txPos, sep="_")
names(hg18ref.exon) =exonsByUniId
genome(hg18ref.exon)=genome(hg18ref)
# intron by tx-GRL
hg18ref.intron=intronsByTranscript(hg18ref,use.name=TRUE)
names(hg18ref.intron) = exonsByUniId
genome(hg18ref.intron)=genome(hg18ref)
# TSS by tx-GR
txTSS = txStart
txTSS[txStrand=="-"] = txEnd[txStrand=="-"]

hg18ref.TSS=GRanges(seqnames=txChr,
            ranges=IRanges(start=txTSS, end=txTSS),
          	strand=txStrand,
			name = txUniId)
genome(hg18ref.TSS)=genome(hg18ref)[names(seqlengths(hg18ref.TSS))]
# promoter by tx-GR
promFlank=1000
hg18ref.promoter=GRanges(seqnames=txChr,
          	ranges=IRanges(start=txTSS-promFlank, end=txTSS+promFlank),
          	strand=txStrand,
			name = txUniId)
genome(hg18ref.promoter)=genome(hg18ref)[names(seqlengths(hg18ref.promoter))]
# cpgi by tx-GR
# need to consider strandness?
# need to set a length limitation?
CpG.file="~/software/seqframework/eRRBS_pipeline/annotation_data/hg18/hg18.CpG.UCSCtable.txt"
chrom.info="~/software/seqframework/eRRBS_pipeline/annotation_data/hg18/hg18.chromInfo.txt"
source("/home/shl2018/software/seqframework/eRRBS_pipeline/scripts/CpGmeth.functions.R")
g.cpg=read.CpG.UCSC.table(GpG.file,chrom.info)

cpgTSSDistTable=.nearest.2bed(hg18ref.TSS, g.cpg)
hg18ref.cpg = GRanges(seqnames=cpgTSSDistTable$seqnames.y,
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

# cpgi shore by tx-GRL
# get CpG island shores
shores=c( flank(hg18ref.cpg,1000),flank(hg18ref.cpg,1000,FALSE) )
hg18ref.shore = split(shores, rep(1:length(hg18ref.cpg),2))
names(hg18ref.shore) = cpgTSSDistTable$name
genome(hg18ref.shore)=genome(hg18ref)[names(seqlengths(hg18ref.cpg))]
# cpgi by tx-GRL
hg18ref.cpg = split(hg18ref.cpg)
names(hg18ref.cpg) = cpgTSSDistTable$name
genome(hg18ref.cpg)=genome(hg18ref)[names(seqlengths(hg18ref.cpg))]

save(hg18ref.exon, hg18ref.intron, hg18ref.promoter, hg18ref.cpg, hg18ref.shore,txPos)