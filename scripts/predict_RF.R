#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
test.data<- args[1]
model <- args[2]
output<- args[3]

# load libraries
library(data.table)
library(plyr)
library(dplyr)
library(mlr)

# load model
print("loading appropriate model")
load(model)
print(fit)

# load test
print("loading and formatting test data")
test<- read.table(test.data, header=TRUE, sep="\t",, na.strings=c("NA",".","na","N/A"), skipNul=FALSE, row.names=NULL)
print(colnames(test))

# making predictions
print("making predictions and outputting results")
pred= predict(fit, newdata= test, type="prob")
write.table(pred$data, sep='\t', quote=FALSE, row.names=FALSE, na=".", output)
