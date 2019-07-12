#!/bin/sh
# help
if [ "$1" == "-h" ] ; then
    echo "Usage= ./predict.bash -a <path to R> -o <output directory> -w < peaks > -m < model >"
    exit 0
fi

# running with options
while getopts ":a:o:w:m:" opt
   do
     case $opt in
		a ) path=$OPTARG;;
		o ) output=$OPTARG;;
    w ) peaks=$OPTARG;;
	m ) model=$OPTARG;;
     esac
done

# predefine regions
$path/Rscript --vanilla ~/BreakCA/bin/predefine_windows.R $peaks $output/windows.bed

# prepare contigs for logistic, elas-net or randomForest
$path/Rscript --vanilla ~/BreakCA/bin/prepare_dataset.R $output/all.positions.tsv $output/windows.bed $output/classifier_id_frame.csv $output/classifier_input.tsv

# make prediction using logistic regression
$path/Rscript --vanilla ~/BreakCA/bin/make_predictions.R $output/classifier_input.tsv $model $output/prediction.txt

echo "predictions made using RF"
