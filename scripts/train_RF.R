#!/usr/bin/env RScript

args <- commandArgs(trailingOnly = TRUE)
training<- args[1]
model<- args[2]
balance<- args[3] # an integer (0 or 1) indicating whether to balance the data
importance<- args[4] # if specified, importance for each of the variables is saved here

# load libraries
library(plyr)
library(dplyr)
library(mlr)
library(parallelMap)
library(parallel)

# load data.frame
print("loading training data into R")
training<- read.table(training, header=TRUE, sep="\t",, na.strings=c("NA",".","na","N/A"), skipNul=FALSE, row.names=NULL)

print("creating training task and making RF learner")

# optimize hyper parameters
# make training task
traintask <- makeClassifTask(data = training, target = colnames(training)[ncol(training)], positive = 1)

# create learner
rf.lrn <- makeLearner("classif.ranger", predict.type = "prob")

if (as.integer(balance)) {
	print("calculating class weights in order to ensure data is balanced when sampled")
	# first, retrieve the inverse of the counts of each of the labels
	w <- 1/table(training[,ncol(training)])
	# calculate probabilities for each label
	w <- w/sum(w)
	# create a vector containing weights instead of labels
	weights <- rep(0, nrow(training))
	for (val in names(w)) {
		weights[training[,ncol(training)] == val] <- w[val]
	}

	# create par.vals
	rf.lrn$par.vals <- list(importance='impurity', verbose=TRUE, case.weights=weights)
} else {
	# create par.vals
	rf.lrn$par.vals <- list(importance='impurity', verbose=TRUE)
}

# # create params
# print("creating params to tune")
# # mtry default: sqrt(number of features)
# # nodesize default: 5
# # ntree default: 500
# params <- makeParamSet(makeIntegerParam("mtry",lower = 4,upper = 10),
#                        makeIntegerParam("nodesize",lower = 10,upper = 20),
#                        makeIntegerParam("ntree", lower = 50, upper=75))

# #set validation strategy; 5-fold cross validation
# rdesc <- makeResampleDesc("CV",iters=3L)

# #set optimization technique
# ctrl <- makeTuneControlRandom(maxit = 3L)

# #create cluster
# parallelStartSocket(20)
# print("creating parallel sockets")

# #start tuning for accuracy
# tune <- tuneParams(learner = rf.lrn, task = traintask, resampling = rdesc,
#                    measures = list(acc), par.set = params, control = ctrl,
#                    show.info = T)
# # tuned params
# tune$x
# # stop parallel
# parallelStop()
# print("tuned params for randomForest have been obtained")

# # make model
# rf.lrn$par.vals<- c(rf.lrn$par.vals, tune$x)
print("training model")
fit= mlr::train(rf.lrn, traintask)

# print out variable importance
if (importance != "") {
	print("recording variable importance:")
	importance_df = as.data.frame(sort(fit$learner.model$variable.importance, decreasing=TRUE))
	print(importance_df)
	names(importance_df) <- c("variable\timportance")
	write.table(importance_df, sep="\t", file=importance, quote=FALSE)
}

# save.data
save.image( model )
