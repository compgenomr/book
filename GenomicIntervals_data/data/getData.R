library(GenomicRanges)
library(rtracklayer)

bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K4me3.fc.signal.bigwig"
k4me3=import(bw.file, 
       which=seqinfo(BigWigFile(bw.file))["chr20"])
k4me3=k4me3[k4me3$score !=0,]
k4me3=keepSeqlevels(k4me3,"chr20")

bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-DNase.fc.signal.bigwig"
dnase=import(bw.file, 
             which=seqinfo(BigWigFile(bw.file))["chr20"])
dnase=dnase[dnase$score !=0,]
dnase=keepSeqlevels(dnase,"chr20")

bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K4me1.fc.signal.bigwig"
k4me1=import(bw.file, 
             which=seqinfo(BigWigFile(bw.file))["chr20"])
k4me1=k4me1[k4me3$score !=0,]
k4me1=keepSeqlevels(k4me1,"chr20")


bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K27ac.fc.signal.bigwig"
H3K27ac=import(bw.file, 
             which=seqinfo(BigWigFile(bw.file))["chr20"])
H3K27ac=H3K27ac[H3K27ac$score !=0,]
H3K27ac=keepSeqlevels(H3K27ac,"chr20")

sml=ScoreMatrix(H3K27ac,prom,strand.aware = TRUE)
plotMeta(sml)


export.bw(k4me3,"H1.ESC.H3K4me3.chr20.bw")
export.bw(k4me1,"H1.ESC.H3K4me1.chr20.bw")
export.bw(dnase,"H1.ESC.dnase.chr20.bw")
export.bw(H3K27ac,"H1.ESC.H3K27ac.chr20.bw")
