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

# function to get reads
getreads= function(bam, peaks){
   what = c("seq","qual")
   if(is.null(peaks)=="TRUE"){param = ScanBamParam(what=what)
   }else{
      peaks= import.bed(peaks)
      peaks=reduce(peaks)
      param = ScanBamParam(which = peaks,what=what)}
   msg=paste("reading bam",bam)
   print(msg)
   bam = readGAlignments(bam, param = param)
   mat= as.data.frame(bam)
   mat$Ncigar= gsub('[[:digit:]]+', '', mat$cigar)
   mat$com= paste(str_sub(mat$Ncigar,1,1),str_sub(mat$Ncigar,-1,-1), sep=";")
   return(mat)}

# write reads into a tsv file
reads<- getreads(bam, peaks)
reads<- reads[!duplicated(reads),]
write_tsv(reads, filename)
