args <- commandArgs(trailingOnly = TRUE)
all.positions<- args[1]
windows<- args[2]
id<- args[3]
save.file<- args[4]
print(args)

# load library
library(data.table)
library(rtracklayer)
library(biovizBase)
library(dplyr)
library(pbapply)
library(readr)

# read all positions
print("read all positions")
all<- as.data.frame(fread(all.positions, sep="\t"))

# make positions for overlaps
print("convert position to GRanges")
seqnames<- sapply(strsplit(all$bpos, ":"), "[",1)
start<- as.numeric(sapply(strsplit(all$bpos, ":"), "[",2))
end<- start
positions<- makeGRangesFromDataFrame(data.frame(seqnames, start, end,
                                                bpos=all$bpos), keep.extra.columns = T)
print(paste("# of position=", length(positions), sep=""))

# get windows
print("read predefined genomic windows")
win<- import.bed(windows)
win$windows<- paste(seqnames(win), start(win), end(win), sep="_")

# get overlaps
print("find overlap with positions")
hits<- findOverlaps(positions, win, type = "within")
new.pos<- positions[queryHits(hits)]
new.win<- win[subjectHits(hits)]
new.pos$win<- mcols(new.win)$windows

# subset positions
new.all<- all[all$bpos %in% mcols(new.pos)$bpos,]

# match and add windows
m<- match(new.all$bpos, mcols(new.pos)$bpos)
new.all$windows<- mcols(new.pos)$win[m]

# start
new.all$start<- as.numeric(sapply(strsplit(new.all$bpos, ":"), "[", 2))
df<- new.all

# function to select the best value from each window
print("select best value / window")

# divide data.frame to right positions
right<- df[, c("windows", "start", "right_pm", "right_psd",
               "right_IC", "right_mean_length",
               "right_sd_length")]

# select top right positions per 20 windows
right_top= right %>%
   group_by(windows) %>%
   top_n(1, right_pm)
right_top= as.data.frame(right_top)
right_top= right_top[!duplicated(right_top$windows),]
colnames(right_top)[c(2)]<- paste("right",colnames(right_top)[c(2)], sep="_")

# divide data.frame left positions
left <- df[, c("windows", "start", "left_pm", "left_psd",
               "left_IC", "left_mean_length",
               "left_sd_length")]

# select top left clip per 20 bps window
left_top= left %>%
   group_by(windows) %>%
   top_n(1, left_pm)
left_top= as.data.frame(left_top)
left_top= left_top[!duplicated(left_top$windows),]
colnames(left_top)[c(2)]<- paste("left",colnames(left_top)[c(2)], sep="_")

# select top insertions per 20bps windows
INS<- df[, c("windows","INS_pm", "INS_psd")]
INS_top= INS %>%
   group_by(windows) %>%
   top_n(1, INS_pm)
INS_top= as.data.frame(INS_top)
INS_top= INS_top[!duplicated(INS_top$windows),]

# select top deletions per 20bps windows
DEL<- df[, c("windows","DEL_pm", "DEL_psd")]
DEL_top= DEL %>%
   group_by(windows) %>%
   top_n(1, DEL_pm)
DEL_top= as.data.frame(DEL_top)
DEL_top= DEL_top[!duplicated(DEL_top$windows),]

# ranked base-pairs windows
print("merge data frames")
ranked<- Reduce(function(x, y) merge(x, y, all=TRUE, by="windows"),
                list(DEL_top, INS_top, left_top, right_top))

# add if a window in N or 1KG:: only for training data
#print("add is.kg")
#ranked$is.kg<- as.character(mcols(win)$type[match(as.character(ranked$windows), as.character(mcols(win)$windows))])
#ranked$is.kg<- factor(ranked$is.kg, levels = unique(ranked$is.kg))

# absolute distance between top right and top left window
ranked$dist<- abs(ranked$right_start - ranked$left_start)
write.csv(ranked[, c("windows", "left_start", "right_start")], id, row.names = F)
ranked<- ranked[,-c(6,12)]

# write tsv
print("output data")
print(head(ranked))
write_tsv(ranked, save.file)

# conditional split for training data:: only for training data and prec-rec curves
# load dataset
#library(caret)
#ranked_out<- ranked
#rownames(ranked_out)<- ranked_out$windows
#ranked_out= ranked_out[, -1]

# create training and test
#Train <- createDataPartition(ranked_out$is.kg, p=0.5, list=FALSE)
#training <- ranked_out[ Train, ]
#testing <- ranked_out[ -Train, ]

# output training and test
#write.table(training, "training.txt")
#write.table(testing, "testing.txt")
