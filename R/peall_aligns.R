#========== PEALL aligns from bam ===========

# load library
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(stringr)

# load annotation reference
load("/scratchLocal01/shl2018/eRRBS/bcdata/myCpG/hg18ref.Rdata")

# function

# load alignment from bam file and get read count and annotation for each gene.
bamToCount=function(bamFile,anno=hg18ref.exon){
# from bam to genomicRanges object
	aligns = readBamGappedAlignments(bamFile)
	
# get count
	counts <- countOverlaps(anno, aligns)
	
# get rpkm
	numBases <- sum(width(anno))
	geneLengthsInKB <- numBases/1000
	millionsMapped <- sum(counts)/1e+06
	rpm <- counts/millionsMapped
	rpkm <- rpm/geneLengthsInKB
	
	return(list(counts=counts,rpkm=rpkm, bamFile=bamFile))
	
}

# end of function

# read file
bamDir = '/scratchLocal01/shl2018/peall/data/'
bamFiles = str_c(bamDir,'s_',1:8,'_hg18.bam')


allCounts = lapply(bamFiles, bamToCount)

seBamFiles = dir('/scratchLocal01/shl2018/all/sortedbam',
full.names=TRUE, pattern='bam')[seq(1,24,2)]

seAllCounts = lapply(seBamFiles, bamToCount)
