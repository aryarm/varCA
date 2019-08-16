args <- commandArgs(trailingOnly = TRUE)
bam<-args[1]
peaks<- args[2]
filename<-args[3]
print(args)

#load dataset
library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)
library(biovizBase)
library(lattice)
library(plyr)
library(dplyr)
library(data.table)
library(reshape)
library(rmutil)
library(lattice)
library(stringr)
library(pbapply)
library(readr)


what = c("seq","qual")
if(is.null(peaks)=="TRUE"){
   param = ScanBamParam(what=what)
}else{
   peaks= import.bed(peaks)
   peaks=reduce(peaks)
   param = ScanBamParam(which = peaks,what=what)
}
msg=paste("reading bam",bam)
print(msg)
reads = readGAlignments(bam, param = param)
reads= as.data.frame(reads)
reads$Ncigar= gsub('[[:digit:]]+', '', reads$cigar)
reads$com= paste(str_sub(reads$Ncigar,1,1),str_sub(reads$Ncigar,-1,-1), sep=";")

# write reads into a tsv file
write_tsv(reads[!duplicated(reads),], filename)
