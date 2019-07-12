args <- commandArgs(trailingOnly = TRUE)
classifier_input<- args[1]
gatk_input<- args[2]
input_windows<- args[3]
filename.output<- args[4]
print(args)

library(data.table)
library(dplyr)
library(plyr)
library(readr)
library(rtracklayer)
library(VariantAnnotation)

# import gatk table
print("import gatk vcf")
gatk<- readVcf(gatk_input)
gr<- granges(gatk)
gr$QD<- info(gatk)$QD
names(gr)<- NULL
gr$name<- paste(seqnames(gr), start(gr), end(gr), sep="_")
var<- gr

# get windows
print("read predefined genomic windows")
win<- import.bed(input_windows)
win$windows<- paste(seqnames(win), start(win), end(win), sep="_")

# overlap
if(length(var)>0){

	hits= findOverlaps(var, win)
	new.var<- as.data.frame(var[queryHits(hits)])
	new.w<- as.data.frame(win[subjectHits(hits)])
	new.w$QD<- new.var$QD

	print("ranking and selecting best QD and FS per window")
	# select best value for QD
	QD.rank<- new.w[, c("windows","QD")]
	QD.top<- QD.rank %>%
		group_by(windows) %>%
		top_n(1, QD) %>% as.data.frame()

	# read classifier input
	input<- fread(classifier_input) %>% as.data.frame()
	ranked<- QD.rank
	print("add QD to input frame")
	input$QD<- ranked$QD[match(input$windows, ranked$windows)]
	input$QD[is.na(input$QD)]<- 0
	print(head(input))}else{
	input<- fread(classifier_input) %>% as.data.frame()
	print("add QD to input frame")
	input$QD<- 0
	print(head(input))
	}	
write_tsv(input, filename.output)

