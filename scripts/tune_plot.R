#!/usr/bin/env RScript

# This R script can be used to visualize the output of the hyperparameter tuning step in train_RF.R

# param1: The path to a TSV containing the results of tuning the hyperparameters
#         For each instance of cross-validation (rows), there should be the following cols: hyperparam vals (2 cols), mean test accuracy
# param2 (optional): The path to a png file in which to store a visualization of the test accuracy for different hyperparam values. If not specified, stdout will be used

args <- commandArgs(trailingOnly = TRUE)
results<- args[1]
out<- args[2]

# load libraries
library(mlr)
library(BBmisc)
library(ggplot2)
library(R.devices)

# load data.frame
print("loading results of hyperparam tuning into R")
results<- read.table(results, header=TRUE, sep="\t", na.strings=c("NA",".","na","N/A"), row.names=NULL)
data = makeS3Obj("HyperParsEffectData", data = results, measures=c("f1.test.mean"), hyperparameters=c("mtry", "min.node.size"), partial=F, nested=F)
plt = plotHyperParsEffect(data, x = "mtry", y = "min.node.size", z = "f1.test.mean", plot.type = "heatmap", interpolate = "regr.earth", show.experiments = TRUE)
# min_plt = min(data$data$acc.test.mean, na.rm = TRUE)
# max_plt = max(data$data$acc.test.mean, na.rm = TRUE)
# med_plt = mean(c(min_plt, max_plt))
# plt = plt + scale_fill_gradient2(breaks = seq(min_plt, max_plt, length.out = 5), low = "blue", mid = "white", high = "red", midpoint = med_plt)
suppressGraphics(ggsave(out, plt))

