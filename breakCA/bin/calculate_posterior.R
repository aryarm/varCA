args <- commandArgs(trailingOnly = TRUE)
counts<-args[1]
filename.posteriors<- args[2]
clipping.info<- args[3]
filename.all<- args[4]
print(args)

# load libraries
library(data.table)
library(plyr)
library(dplyr)
library(rmutil)
library(readr)
library(rtracklayer)

# read counts file
counts<- as.data.frame(fread(counts, sep="\t"))
print(paste("# of base positions=",nrow(counts), sep=""))

# keep positions with minimum of 10 reads
df= subset(counts, N >= 10)
print(paste("# of base positions=",nrow(df), sep=""))

# functions to calculate posterior mean and standard deviation
#Mean
calcMean <- function(successes, total, a, b) {
   calcBetaMean <-
      function(aa, bb) {
         BetaMean <- (aa) / (aa + bb)
         return(BetaMean)
         
      }
   posterior_a = a + successes
   posterior_b = b + total - successes
   posterior_mean  <- calcBetaMean(posterior_a, posterior_b)
   return(posterior_mean)
}

#Standard deviation
calcSd <- function(successes, total, a, b) {
   calcBetaSd   <-
      function(aa, bb) {
         BetaSd <-
            sqrt((aa * bb) / (((aa + bb) ^ 2) * (aa + bb + 1)))
         return(BetaSd)
         
      }
   posterior_a = a + successes
   posterior_b = b + total - successes
   posterior_sd <- calcBetaSd(posterior_a, posterior_b)
   return(posterior_sd)
}

# estimate beta prior
estBetaParams <- function(mu, var) {
   alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
   beta <- alpha * (1 / mu - 1)
   return(params = list(alpha = alpha, beta = beta))
}

features <- c("right", "left", "INS", "DEL")

# function to estimate posteriors
my.func<- function(i){
   dat<- df[, c("bpos", "N", features[i])]
   dat<- dat[,c(1,3,2)]
   
   # get beta prior
   mu<- mean(dat[,2]/dat[,3])
   var<- var(dat[,2]/dat[,3])
   beta.prior= estBetaParams(mu, var)
   
   # subset dat
   test<- dat[dat[,2]>0,]
   
   # get pm and psd
   test$pm= apply(test[, c(2,3)], 1, function(x) calcMean(x[1], x[2], beta.prior$alpha, beta.prior$beta))
   
   test$psd= apply(test[, c(2,3)], 1, function(x) calcSd(x[1], x[2], beta.prior$alpha, beta.prior$beta))
   
   colnames(test)[c(4,5)]<- paste(features[i], colnames(test)[c(4,5)], sep="_")
   
   final<- test[, c(1,4,5)]
   
   rownames(final)<- NULL
   
   return(final)
}
pos.list<- llply(1:length(features), my.func, .progress = progress_text(char="+"))

# compile posteriors
print("output posteriors")
posteriors= Reduce(function(x, y) merge(x, y, all=TRUE, by="bpos"), pos.list)
posteriors[is.na(posteriors)]<- 0
write_tsv(posteriors,filename.posteriors)

# compile non indel and indel positions
print("add non tested positions")
non<- df$bpos[!df$bpos %in% posteriors$bpos]
non.df<- data.frame(matrix(nrow = length(non), ncol = ncol(posteriors)))
colnames(non.df)<- colnames(posteriors)
non.df$bpos<- non
non.df[is.na(non.df)]<- 0
input<- rbind(posteriors, non.df)
setDT(input, key="bpos")

# load clipping info
print("add clip information")
if(!file.exists(clipping.info)){
   
   # create empty data frame when no clipping is present
   setcol<- c("right_mean_length","right_sd_length","right_IC","left_mean_length",
              "left_sd_length", "left_IC", "bpos")
   clip_info<- as.data.frame(matrix(nrow= nrow(input), ncol= length(setcol),data = 0))
   colnames(clip_info)<- setcol
   clip_info$bpos<- input$bpos
   setDT(clip_info, key = "bpos")
   
   # merge data.table
   new.input<- merge(input, clip_info, by="bpos", all.x= T)
   new.input<- as.data.frame(new.input)
   new.input[is.na(new.input)]<- 0
   write_tsv(new.input, filename.all)
   
}else{
   
   # read clip info file
   clip_info<- read.table(clipping.info)
   setDT(clip_info, key = "bpos")

   # merge data.table
   new.input<- merge(input, clip_info, by="bpos", all.x= T)
   new.input<- as.data.frame(new.input)
   new.input[is.na(new.input)]<- 0
   write_tsv(new.input, filename.all)
   
}




