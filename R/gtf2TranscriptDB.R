# GTF to transcriptDB
# The code is changed based on gist 986084 by DarioS.
# The difference: if exon_number is not assinged to gtf,
# The current script will calculate the exon_number to gtf.db
# grep NM_ genes.gtf > NMgenes.gtf

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


# get regional ranges
txdb=gtf.db
tx.db = transcripts(txdb)
exon.db=exonsBy(txdb,"tx", use.name=TRUE)
intron.db=intronsByTranscript(txdb,use.name=TRUE)
cds.db=cdsBy(txdb,"tx", use.name=TRUE)
threeutr.db = threeUTRsByTranscript(txdb, use.name=TRUE)
fiveutr.db = fiveUTRsByTranscript(txdb, use.name=TRUE)

# tss
txStart=start(tx.db)
txEnd =end(tx.db)
txStrand=as.character(strand(tx.db))
txTSS = txStart
txTSS[txStrand=="-"] = txEnd[txStrand=="-"]

tss.db=GRanges(seqnames=txChr,
            ranges=IRanges(start=txTSS, end=txTSS),
            strand=txStrand,
			      name = txUniId)

# cpgi by tx-GR
# need to consider strandness?
# need to set a length limitation?
CpG.file="~/software/seqframework/eRRBS_pipeline/annotation_data/hg18/hg18.CpG.UCSCtable.txt"
chrom.info="~/software/seqframework/eRRBS_pipeline/annotation_data/hg18/hg18.chromInfo.txt"
source("/home/shl2018/software/seqframework/eRRBS_pipeline/scripts/CpGmeth.functions.R")
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

hg18ref=list(assembly='hg18',txdb=txdb, tx=tx.db, exon=exon.db, intron=intron.db,
              tss=tss.db, cpgi=cpgi.db, shore=shore.db)
save(hg18ref,file="~/annotation/hg18/hg18ref.Rdata")

# read in gtf file directly
 NMgenes.gtf=read.table("NMgenes.gtf", colClasses=c("character","character","character", "integer","integer", "character","character","character","character"), stringsAsFactors=FALSE,sep="\t",header=FALSE)
NMgenes.gtf$V9=gsub(";", "",NMgenes.gtf$V9)
stringTotxid=function(x){unlist(strsplit(x," "))[grep('NM_',unlist(strsplit(x," ")))]}
NMgenes.gtf$txid=stringTotxid(NMgenes.gtf$V9)
