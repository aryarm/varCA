args <- commandArgs(trailingOnly = TRUE)
soft<-args[1]
filename.reads<- args[2]
filename.info<-args[3]
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

print("initialising info gathering...................")

# read sc read tsv from tmp dir
dat<- as.data.frame(fread(soft))
print("read soft clipped reads")

# cigar manipulations
cigar= dat$cigar
opt=explodeCigarOps(cigar, ops=CIGAR_OPS)
len=explodeCigarOpLengths(cigar, ops=CIGAR_OPS)
print("CIGAR exploded")

# add opt and len info to clipped sequences
add.info = function(x) {
   a = data.frame(opt = unlist(opt[[x]]), length = unlist(len[[x]]))
   a$combine = paste(a$length, a$opt, sep = "")
   b = data.frame(com = paste(a[1, 1], a[nrow(a), 1], sep = ";"),
                  len = paste(a[1, 2], a[nrow(a), 2], sep = ";"))
   return(b)
}
cl<- makeCluster(4)
clusterExport(cl, c("dat", "add.info","opt", "len"))
info.list<- pblapply(1:length(opt), add.info, cl=cl)
stopCluster(cl)
rm(cl)

# add clip length info for cigar manipulations
s<- rbindlist(info.list)
mat= cbind(dat,s)
start= subset(mat, dist2start==0)
nstart= subset(start, com=="S;M"| com=="S;S")
nstart$cliplen= as.numeric(sapply(strsplit(as.character(nstart$len), ";"),"[",1))
end= subset(mat, dist2end==0)
nend= subset(end, com=="M;S" | com=="S;S")
nend$cliplen= as.numeric(sapply(strsplit(as.character(nend$len), ";"),"[",2))
df1= rbind(nstart, nend)
df1$index<- paste("index", rownames(df1),sep=":")
print("clip info added")

# get clip sequences
getsequence <- function(i) {
   dl = df1[i, ]
   if (dl$dist2start == 0) {
      strReverse <- function(x)
         sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
      seq <- substr(as.character(dl$seq), dl$dist2start, dl$cliplen)
      scseq <- strReverse(seq)
   } else{
      scseq <- substr(
         as.character(dl$seq),
         nchar(as.character(dl$seq)) - dl$cliplen + 1,
         nchar(as.character(dl$seq))
      )
   }
   names(scseq)<- dl$index
   scseq
}
cl<- makeCluster(4)
clusterExport(cl, c("df1", "getsequence"))
clusterEvalQ(cl, library(stringr))
clip.seq.list<- unlist(pblapply(1:nrow(df1), getsequence, cl=cl))
print("gathered soft-clip sequences")

# stop cluster
stopCluster(cl)
rm(cl)

# write backup files
m= match(df1$index, names(clip.seq.list))
df1$scseq<- clip.seq.list[m]
dl= df1 # very important step
write_tsv(dl, filename.reads)
print("written soft-clipped sequences in file")

# get base positions
bp= unique(dl$bpos)

# function to calculate info content
getinfo <- function(i){
   # subset matrix by bpos
    mat = dl[dl$bpos %in% bp[i], ]
    mat = mat[order(mat$cliplen, decreasing = T), ]
   
    # get unique clipping types
    clip_type= c("right", "left")
    
    # function to calculate seperate IC for right and left clips
    my.func= function(n){
		clip= clip_type[n]
		cseq= as.character(mat$scseq[mat$n.com %in% clip])
       
		# calculate mean length and sd
		if(length(cseq)>2){
		   mean.length<- mean(nchar(cseq))
		   sd.length<- sd(nchar(cseq))
		}else{
		   mean.length<- mean(nchar(cseq))
		   sd.length<- 0
		}
		
       # loop to create info.content object
       if (length(cseq) > 1) {
          # convert letters to numbers
          func <- function(i) {
             a <- strsplit(cseq[i], "")[[1]]
             n <- as.numeric(chartr("ACGT", "1234", a))
             as.data.frame(t(as.matrix(n)))
          }
          sim.seqs <-do.call(rbind.fill, lapply(1:length(cseq), func))
          sim.seqs<- as.matrix(sim.seqs)
          colnames(sim.seqs)<- NULL
          
          # information needed to calculate infor content
          ref.len <- max(nchar(cseq))
          n.nuc <- 4
          pseudo.count <- 0.1
          n.clipped.reads <- length(cseq)
          n.obs <- apply(!is.na(sim.seqs), 2, sum)
          
          # create entropy and max entropy rows
          entropy <- rep(NA, ref.len)
          max.entropy <- rep(NA, ref.len)
          
          # calculate entropy of each column
          for (i in 1:ref.len) {
             if (n.obs[i] > 1) {
                # calculate entropy of this column, which is sum(-p*log2(p))
                # where p is the proportion of each nucleotide class
                # include a pseudo-count in the calculation of p, because we can't work with 0s
                x <- rep(0.0, n.nuc)
                for (j in 1:n.nuc) {
                   x[j] <- sum(sim.seqs[, i] == j, na.rm = T)
                }
                prp <- (x + pseudo.count) / sum(x + pseudo.count)
                entropy[i] <- sum(-prp * log2(prp))
                
                # calculate maximum possible entropy with this number of observed nucleotides
                # assume that observations are evenly spread across all 4 nucleotide categories
                # but that (n.obs modulu n.categories) have 1 more count.
                obs.per.nuc <- floor(n.obs[i] / n.nuc)
                remainder <- n.obs[i] %% n.nuc
                x <- rep(obs.per.nuc, n.nuc)
                x[1:remainder] <- x[1:remainder] + 1
                prp <- (x + pseudo.count) / sum(x + pseudo.count)
                max.entropy[i] <- sum(-prp * log2(prp))
             }
          }
          # calculate difference from maximum entropy, call this information content
          # of the sequence
          info.content <- max.entropy - entropy
          total.info.content <- sum(info.content, na.rm = T)
          
          # calculate amount of information by nucleotide across sites where
          # info content computed
          n.bases <- sum(as.numeric(!is.na(info.content)) * n.obs)
          info.by.base <- total.info.content / n.bases
       } else{
          info.by.base<- 0
       }
		fin<- c(mean.length, sd.length, info.by.base)
		names<- paste(clip,"_",c("mean_length", "sd_length", "IC"),sep="") 
		names(fin)<- names
		return(fin)
    }
    list<- lapply(1:length(clip_type), my.func)
    all.info<- unlist(list)
    all.info[is.na(all.info)]<- 0
    return(c(all.info, bpos=unique(mat$bpos)))
}

#make cluster, note when running on multiple clusters cannot add base position as bp[i]
cl= makeCluster(4)
clusterExport(cl, c("bp", "getinfo", "dl"))
clusterEvalQ(cl, c(library(data.table), library(plyr)))
list<- pblapply(1:length(bp), getinfo, cl=cl)
dat<- data.frame(do.call(rbind, list))

# stop cluster
stopCluster(cl)
rm(cl)

# outputting clipping information
write.table(dat, filename.info)
print("output clipping information")


