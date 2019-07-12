args <- commandArgs(trailingOnly = TRUE)
peaks<- args[1]
filename.bed<- args[2]
print(args)

# load libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(plyr)
library(biovizBase)
library(pbapply)

# genome
print("getting chromosomal length")
genome= BSgenome.Hsapiens.UCSC.hg19
info= as.data.frame(seqinfo(genome))
info= info[1:23,]

# load peaks
print("reading peaks")
peaks<- import.bed(peaks)
print(head(peaks))

# parse info rownames according
if(length(grep("chr",as.character(seqnames(peaks))))==0){
   rownames(info)<- sapply(strsplit(rownames(info), "chr"), "[", 2)
}else{
   info<- info
}
print(info)

# function
windows<- function(i){
   seq.info= Seqinfo(seqnames = rownames(info[i,]),
                     seqlengths = info[i,]$seqlengths,
                     isCircular = info[i,]$isCircular ,
                     genome = "hg19")

   windows<- tileGenome(seq.info,tilewidth = 20,
                        cut.last.tile.in.chrom = TRUE)

   out<- subsetByOverlaps(windows, peaks)
   out<- out[!duplicated(out)]
   out
}


# make clusters
cl<- makeCluster(2)
clusterExport(cl, c("peaks", "windows", "info"))
clusterEvalQ(cl, c(library(BSgenome.Hsapiens.UCSC.hg19), library(biovizBase)))

# run function
window.list<- pblapply(1:nrow(info), windows, cl=cl)
final.windows<- flatGrl(GRangesList(window.list))
print(final.windows)

# stop cluster
stopCluster(cl)
rm(cl)

# export bed
print("writing windows into bed")
export.bed(final.windows, filename.bed)
