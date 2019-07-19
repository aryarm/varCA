#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
training<- args[1]
model<- args[2]

# load libraries
library(plyr)
library(dplyr)
library(mlr)
library(parallelMap)
library(parallel)

# load data.frame
training<- read.table(training, header=TRUE, sep="\t",, na.strings=c("NA",".","na","N/A"), skipNul=FALSE, row.names=NULL)

# optimize hyper parameters
# make training task
traintask <- makeClassifTask(data = training, target = colnames(training)[ncol(training)], positive = 1)

# create learner
rf.lrn <- makeLearner("classif.randomForest", predict.type = "prob")

# create par.vals
rf.lrn$par.vals <- list(importance=TRUE)

# create params
print("creating params to tune")
params <- makeParamSet(makeIntegerParam("mtry",lower = 10,upper = 50),
                       makeIntegerParam("nodesize",lower = 10,upper = 50),
                       makeIntegerParam("ntree", lower = 100, upper=500))

#set validation strategy; 5-fold cross validation
rdesc <- makeResampleDesc("CV",iters=3L)

#set optimization technique
ctrl <- makeTuneControlRandom(maxit = 3L)

#create cluster
parallelStartSocket(20)
print("creating parallel sockets")

#start tuning for accuracy
tune <- tuneParams(learner = rf.lrn, task = traintask, resampling = rdesc,
                   measures = list(acc), par.set = params, control = ctrl,
                   show.info = T)
# tuned params
tune$x
# stop parallel
parallelStop()
print("tuned params for randomForest is obtained")

# make model
rf.lrn$par.vals<- c(list(importance=TRUE), tune$x)
fit= mlr::train(rf.lrn, traintask)

# save.data
save.image( model )
