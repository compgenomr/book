
## ---- createGR ----
library(GenomicRanges)
gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           strand=c("+","-","-")
)
gr
# subset like a data frame
gr[1:2,]


## ---- createGRwMetadata ----
gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           names=c("id1","id3","id2"),
           scores=c(100,90,50)
)
# or add it later (replaces the existing meta data)
mcols(gr)=DataFrame(name2=c("pax6","meis1","zic4"),
                    score2=c(1,2,3))

gr=GRanges(seqnames=c("chr1","chr2","chr2"),
           ranges=IRanges(start=c(50,150,200),end=c(100,200,300)),
           names=c("id1","id3","id2"),
           scores=c(100,90,50)
)

# or appends to existing meta data
mcols(gr)=cbind(mcols(gr),
                          DataFrame(name2=c("pax6","meis1","zic4")) )
gr
# elementMetadata() and values() do the same things
elementMetadata(gr)
values(gr)

# you may also add metadata using the $ operator, as for data frames
gr$name3 = c("A","C", "B")
gr

## ---- convertDataframe2gr ----
# read CpGi data set
cpgi.df = read.table("../data/cpgi.hg19.chr21.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
# remove chr names with "_"
cpgi.df =cpgi.df [grep("_",cpgi.df[,1],invert=TRUE),]

cpgi.gr=GRanges(seqnames=cpgi.df[,1],
                ranges=IRanges(start=cpgi.df[,2],
                              end=cpgi.df[,3]))

## ---- convertDataframe2grTSS ----
# read refseq file
ref.df = read.table("../data/refseq.hg19.chr21.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
ref.gr=GRanges(seqnames=ref.df[,1],
               ranges=IRanges(start=ref.df[,2],
                              end=ref.df[,3]),
               strand=ref.df[,6],name=ref.df[,4])
# get TSS
tss.gr=ref.gr
# end of the + strand genes must be equalized to start pos
end(tss.gr[strand(tss.gr)=="+",])  =start(tss.gr[strand(tss.gr)=="+",])
# startof the - strand genes must be equalized to end pos
start(tss.gr[strand(tss.gr)=="-",])=end(tss.gr[strand(tss.gr)=="-",])
# remove duplicated TSSes ie alternative transcripts
# this keeps the first instance and removes duplicates
tss.gr=tss.gr[!duplicated(tss.gr),]

## ---- importbed_rtracklayer ----
require(rtracklayer)
import.bed("../data/refseq.hg19.chr21.bed")


## ---- importFromUCSC ----
require(rtracklayer)
session <- browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
genome(session) <- "mm9"
## choose CpG island track on chr12
query <- ucscTableQuery(session, track="CpG Islands",table="cpgIslandExt",
        range=GRangesForUCSCGenome("mm9", "chr12"))
## get the GRanges object for the track
track(query)


## ---- genomicFeaturesImport ----
require(rtracklayer)
require(GenomicFeatures)
transcript_ids <- c(
"uc009uzf.1",
"uc009uzg.1",
"uc009uzh.1",
"uc009uzi.1",
"uc009uzj.1"
)

txdb2 <- makeTranscriptDbFromUCSC(genome="mm9", tablename="knownGene",
transcript_ids=transcript_ids)
txdb2


## ---- findPeakwithCpGi ----
library(genomation)
pk1.gr=readBroadPeak("../data/wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz")

pk1.gr=pk1.gr[seqnames(pk1.gr)=="chr21",]
# get the peaks that overlap with CpG
# islands
subsetByOverlaps(pk1.gr,cpgi.gr)


## ---- countOverlaps ----
#count the peaks that
# overlap with CpG islands
counts=countOverlaps(pk1.gr,cpgi.gr)
head(counts)



## ---- findOverlaps ----
findOverlaps(pk1.gr,cpgi.gr)

## ---- findNearest ----
# find nearest CpGi to each TSS
n.ind=nearest(pk1.gr,cpgi.gr)
# get distance to nearest
dists=distanceToNearest(pk1.gr,cpgi.gr,select="arbitrary")
dists

# histogram of the distances to nearest TSS
dist2plot=mcols(dists)[,1]
hist(log10(dist2plot),xlab="log10(dist to nearest TSS)",
     main="Distances")


## ---- findCanonical ----
  


## ---- cannonicalCoverage ----

## ---- qualityCheckFastq ----
library(qrqc)
s.fastq <- readSeqFile(system.file("extdata", "test.fastq", package="qrqc"))
qualPlot(s.fastq)


## ---- countBam ----
# regions of interest
# promoters on chr21
promoter.gr=tss.gr
start(promoter.gr)=start(promoter.gr)-1000
end(promoter.gr)  =end(promoter.gr)+1000
promoter.gr=promoter.gr[seqnames(promoter.gr)=="chr21"]

library(Rsamtools)
bamfile="../data/wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam"

# get reads for regions of interest from the bam file
param <- ScanBamParam(which=promoter.gr)
counts=countBam(bamfile, param=param)

## ---- readGAlignments ----
library(GenomicAlignments)
alns <- readGAlignments(bamfile, param=param)

## ---- getCoverageFromAln ---
covs=coverage(alns) # get coverage vectors
covs

## ---- getCoverageFromBam ---
covs=coverage(bamfile, param=param) # get coverage vectors


## ---- readChrbyChr
library(GenomicAlignments)
chrs=c("chr21")
for( chr in chrs){
  
  
  
  # get the parameters to scan the bam file
  # this this should read everthing in a given chr
  param <- ScanBamParam(which=GRanges(seqnames=chr,
                                      IRanges(start=1L,end=500000000L)
  ))
  alns <- readGAlignments(signal, param=param) # see ?readGAlignments
  aln.df <- as.data.frame(alns) # get data.frame
  aln.gr <- as(alns,"GRanges") # get GRanges
  
}



## ---- getRleFromBigWig ---
library(rtracklayer)
# File from ENCODE ChIP-seq tracks
bwFile="../data/wgEncodeHaibTfbsA549.chr21.bw"
bw.gr=import(bwFile, which=promoter.gr) # get coverage vectors
bw.gr



## ---- BigWigCov ---
cov.bw=coverage(bw.gr,weight = "score")

# or get this directly from
cov.bw=import(bwFile, which=promoter.gr,as = "RleList")


## ---- getViews ---
myViews=Views(cov.bw,as(promoter.gr,"RangesList")) # get subsets of coverage
# there is a views object for each chromosome
myViews
myViews[[1]]
# get the coverage vector from the 5th view and plot
plot(myViews[[1]][[5]],type="l")

## ---- viewMeans ---
# get the mean of the views
head(
  viewMeans(myViews[[1]])
)

# get the max of the views
head(
  viewMaxs(myViews[[1]])
)

