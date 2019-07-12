args <- commandArgs(trailingOnly = TRUE)
test.data<- args[1]
model <- args[2]
output<- args[3]

# load libraries
library(data.table)
library(dplyr)
library(mlr)
library(plyr)

# load model
print("loading appropriate model")
load(model)
print(fit)

# load test
print("loading and formatting test data")
test<- as.data.frame(fread(test.data))
rownames(test)<- test$windows
test<- test[,-1]
print(head(test))

# subset test data
print("subsetting test data to relevant positions")
selcol<- c("DEL_pm", "INS_pm", "right_pm", "left_pm")
pred.input<- test[rowSums(test[, selcol])>0,]

# making predictions
print("making predictions and outputting results")
pred= predict(fit, newdata= pred.input, type="prob")
write.table(pred$data, output)
