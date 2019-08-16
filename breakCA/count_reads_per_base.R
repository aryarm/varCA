args <- commandArgs(trailingOnly = TRUE)
reads<-args[1]
insertion<- args[2]
deletion<- args[3]
pileup<- args[4]
# filename.reads<- args[5]
filename.counts<-args[5]
print(args)

#load dataset
library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
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

print("initialising counts...................")

# colnames for the pileup files
set.cols<- c("seqnames", "start", "pileup")

# read insertion containing reads
print("reading insertion pileup")
if(file.info(insertion)$size >0){
   
   # read insertion pileup file
   print("count insertions per base")
   ins<- as.data.frame(fread(insertion, sep="\t"))
   colnames(ins)<- set.cols
   ins$pileup<- as.character(ins$pileup)
   ins$INS<- str_count(ins$pileup,"\\+")
   
   # add end co-ordinates of the insertion
   ins$end<- ins$start+1
   
   # add insertion sizes
   print("get insertion size")
   ins.size<- function(i){
      p<- as.character(ins$pileup[i])
      return(str_match_all(string = p, "\\+[0-9]+") %>% 
                unlist %>% unique %>% as.numeric %>% 
                abs %>% mean %>% floor)
   }
   ins$size<- unlist(lapply(1:nrow(ins), ins.size))
   
   # add identifier
   ins$bpos<- paste(ins$seqnames, ins$start, sep=":")
}else{
   ins<- data.frame(NULL)
}

# read deletion pileups
print ("reading deletion pileup")
if(file.info(deletion)$size >0){
   
   # read deletions
   print("count deletion per base")
   del<- as.data.frame(fread(deletion, sep="\t"))
   colnames(del)<- set.cols
   del$pileup<- as.character(del$pileup)
   del$DEL<- str_count(del$pileup,"\\-")
   
   # get end position for deletion
   print ("get deletion sizes")
   del.size<- function(i){
      p<- as.character(del$pileup[i])
      return(str_match_all(string = p, "\\-[0-9]+") %>% 
                unlist %>% unique %>% as.numeric %>% 
                abs %>% mean %>% floor)
   }
   del$size<- unlist(lapply(1:nrow(del), del.size))
   
   # add identifier:: DEL sites have 2 identifiers i.e. the start and end position
   del$end= del$start + del$size
   del$start.bpos= paste(del$seqnames, del$start, sep=":")
   del$end.bpos= paste(del$seqnames, del$end, sep=":")
}else{
   del<- data.frame(NULL)
}

# get base positions within peaks
print("getting all pileup positions")
pileup <- as.data.frame(fread(pileup, sep="\t")[,c(1,2)])
colnames(pileup)<- c("seqnames", "start")
pileup$end<- pileup$start
pileup$bpos<- paste(pileup$seqnames, pileup$start, sep=":")

# remove unnecessary chr
pileup<- pileup[!pileup$seqnames %in% c("chrM", "chrY", "chrMT", "M", "MT", "Y"),]
rownames(pileup)<- NULL

# make subject positions
subj<- makeGRangesFromDataFrame(pileup, keep.extra.columns = T)

# reads
# make Granges
print("importing reads and converting to Granges")
f= makeGRangesFromDataFrame(as.data.frame(fread(reads, sep="\t",fill=TRUE)), keep.extra.columns = T)

# # get reads which contain soft-clipped bases
# soft.clipped<- f[grep("S", f$cigar)]

# # process soft-clipped
# process_clipped<- function(query, subject){
   
#    # for each chunk assign SC reads to bases
#    print("overlap soft clipped and base positions")
#    hits<- findOverlaps(query, subject)
#    pos<- query[queryHits(hits)]
#    pos$subj.pos= start(subject[subjectHits(hits)])
   
#    # parse positions and select correct reads
#    df<- as.data.frame(pos)
#    df<- df[!duplicated(df),]
   
#    #calculate distance from start and end of the reads
#    df$dist2start= as.numeric(df$subj.pos) - as.numeric(df$start)
#    df$dist2end=as.numeric(df$subj.pos) - as.numeric(df$end)
#    df$bpos<- paste(df$seqnames, df$subj.pos, sep=":")
   
#    # remove positions where distance to start or end in non-zero
#    sc= subset(df,dist2start==0|dist2end==0)
   
#    # true positions for each type of reads
#    sm= subset(sc,dist2start==0 & com=="S;M")
#    ms= subset(sc,dist2end==0 & com=="M;S")
#    ss= subset(sc,com=="S;S")
   
#    # make new true count dataframe
#    n.sc<- rbind(sm, ms, ss)
   
#    # make sure the important fields are characters
#    n.sc$bpos= as.character(n.sc$bpos)
#    n.sc$com=as.character(n.sc$com)
   
#    # assign S.S type reads correctly
#    print("output soft clipped reads")
#    n.sc$n.com[(n.sc$dist2start==0)]<- "left"
#    n.sc$n.com[(n.sc$dist2end==0)]<- "right"
   
#    return(n.sc)
# }

# # run soft-clipped processing function
# if(length(soft.clipped) > 0){
#    n.sc<- process_clipped(query = soft.clipped, subject = subj)
#    write_tsv(n.sc, filename.reads)
# }else{
#    n.sc<- data.frame(NULL)
# }

# # function count reads per soft-clipped base
# count_clipped<- function(data){
#    counts= data.frame(with(data, table(bpos,n.com)))
#    SC = dcast(counts , bpos~n.com,value.var="Freq")
#    return(SC)
# }

# # running count reads for soft-clipped reads function
# if(nrow(n.sc)> 0){
#    SC<- count_clipped(n.sc)
#    cols<- c("bpos", "left", "right")
#    missing.cols<- cols[!cols %in% colnames(SC)]
#    if(length(missing.cols)> 0){
#       SC[[missing.cols]]<- 0
#    }else{
#       SC<- SC
#    }
# }else{
#    SC<- data.frame(NULL)
# }

# add total counts
print("count total reads per base")
subj$N<- countOverlaps(subj,f)

# make count matrix
mat<- as.data.frame(mcols(subj)[, c("bpos", "N")])

# match and add ins
if(nrow(ins)> 0){
   mat$INS<- ins$INS[match(mat$bpos, ins$bpos)]
}else{
   mat$INS<- 0
}

# match and add del
if(nrow(del)> 0){
   dl<- data.frame(bpos= c(del$start.bpos, del$end.bpos),
                   DEL= rep(del$DEL,2))
   mat$DEL<- dl$DEL[match(mat$bpos, dl$bpos)]
}else{
   mat$DEL<- 0
}

mat[is.na(mat)]<- 0

# # match and add SC
# if (nrow(SC)>0){
#    m= match(mat$bpos, SC$bpos)
#    mat$left<- SC$left[m]
#    mat$right<- SC$right[m]
#    mat[is.na(mat)]<- 0
#    mat$SC<- rowSums(mat[,c("left", "right")])
# }else{
#    mat[is.na(mat)]<- 0
#    mat$left<-0
#    mat$right<- 0
#    mat$SC<- 0
# }

# write count data
print("output counts")
# mat<- mat[, c("bpos", "N", "INS", "DEL", "left", "right", "SC")]
mat<- mat[, c("bpos", "N", "INS", "DEL")]
write_tsv(mat, filename.counts)



























