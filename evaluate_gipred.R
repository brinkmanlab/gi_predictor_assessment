# Author: Claire Bertelli (claire.bertelli@gmail.com)
# Last maintained: April 2017

# Functions to evaluate genomic island predictors based on a positive and negative reference dataset,
# as well as a test dataset
###################
# requires three folders containing each the positive, negative and test dataset
# in a 3 column format (name start end) in one file per genome


# load useful libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(ggplot2)

# function to evaluate any prediction in 3 column format 
##############
# All positive, negative, and predicted dataset must have the same structure: 
# a 3 column tab-delimited file with gi_name, start, end - without header for each genome
# each file is named after the id_gis.txt these inputs must be in separate subfolders
# files for the negative dataset must be id_neg.txt

# For several predictors in a vector
evaluate_GRange_pred <- function(predfolders, ids, posfolder, negfolder) {
  allresults <- lapply(predfolders, evaluate_GRange_onePred, ids=ids, posfolder=posfolder, negfolder=negfolder)
  allGRangePred_eval <- do.call(rbind, allresults)
  return(allGRangePred_eval)
}

# For several ids, one predictor
evaluate_GRange_onePred <- function(predfolder, ids, posfolder, negfolder) {
  allresults <- lapply(ids, evaluate_oneGRange_pred, predfolder=predfolder, posfolder=posfolder, negfolder=negfolder)
  allGRange_eval <- do.call(rbind, allresults)
  GRange_eval <- cbind("predictor"=rep(predfolder, nrow(allGRange_eval)), allGRange_eval)
  return(GRange_eval)
}
 
# One id, one predictor 
evaluate_oneGRange_pred <- function(id, predfolder, posfolder, negfolder) {
  print(paste0("Running analysis for : ", id))
  # initialize with no value, so we dont get in trouble if empty
  pred_gis <- data.frame("V1"="GI_0", "V2"=0, "V3"=0)
  posds <- data.frame("V1"="GI_0", "V2"=0, "V3"=0)
  negds <- data.frame("V1"="GI_0", "V2"=0, "V3"=0)
  # now read the real values from file, if they exist
  try(posds <- read.table(paste0(posfolder,"/", id, "_gis.txt"), header=F, sep="\t", as.is=T))
  try(negds <- read.table(paste0(negfolder,"/", id, "_neg.txt"), header=F, sep="\t", as.is=T))
  try(pred_gis <- read.table(paste0(predfolder,"/", id, "_gis.txt"), header=F, sep="\t", as.is=T))
  accstats <- stat_GRanges(pred_gis, posds, negds)
  aid <- rep(id, nrow(accstats))
  result <- cbind(aid, accstats)
  return(result)
}


# function that calculates the recall, precision, accuracy, fscore and Matthew's coefficient 
# for a prediction, given a positive and a negative dataset of genomic ranges, from evaluate_prediction 
#############
# uses IRanges

stat_GRanges <- function(predgis, oneposds, onenegds) {
  #transforming tables into IRanges object
  oneposGR <- IRanges(start=oneposds[,2], end=oneposds[,3])
  predgisGR <- IRanges(start=predgis[,2], end=predgis[,3])
  onenegGR <- IRanges(start=onenegds[,2], end=onenegds[,3])
  
  #finding overlaps between predictions and positive or negative dataset
  gis_overlap <- findOverlaps(predgisGR, oneposGR,minoverlap=1, maxgap = 0L, type="any")
  gis_ranges <- ranges(gis_overlap, predgisGR, oneposGR)
  neg_overlap <- findOverlaps(predgisGR, onenegGR, minoverlap=1, maxgap =0L, type="any")
  neg_ranges <- ranges(neg_overlap, predgisGR, onenegGR)
  
  #calculating bp in confusion matrix as well as evaluation measurements
  truepos <- as.numeric(sum(width(gis_ranges)))
  falsepos <- as.numeric(sum(width(neg_ranges)))
  if (start(oneposGR)[1]!=0) { falseneg <- as.numeric(sum(width(oneposGR))-truepos)
  } else { falseneg <- 0 }
  trueneg <- as.numeric(sum(width(onenegGR))-falsepos)
  accuracy <- (truepos+trueneg)/(truepos+falsepos+falseneg+trueneg)
  recall <- truepos/(truepos+falseneg)
  # we want to give a precision of 1 when there are no truepos but no falsepos either
  # since the software is being conservative. O if the truepos == 0 but there are falsepos.
  if (truepos==0 & falsepos==0) {
    precision <- 1
  } else {
    precision <- truepos/(truepos+falsepos)
  }
  fscore <- 2*truepos/(2*truepos+falsepos+falseneg)
  mcc <- ((truepos*trueneg)-(falsepos*falseneg))/sqrt((truepos+falsepos)*(truepos+falseneg)*(trueneg+falsepos)*(trueneg+falseneg))
  
  # returning results as a dataframe
  metrics <- c("recall", "precision", "accuracy", "fscore", "mcc")
  value <- c(recall, precision, accuracy, fscore, mcc)
  result<- cbind(metrics, value)
  return(result)
}

# Summarizes the accuracy, recall and precision over a list of genomes for the different methods tested vs a reference dataset of gis
# Removes genomes with a different accession number than in the original dataset, indicated by () in the column named IV_jobID_raw
################
# input are:
# metadata, the table with the ivjob
# ref, a string indicating the dataset used as a reference, to adapt the name of the output file

summarize_evaluation_gis <- function(metadata, ref){
  dirs <- dir(".", pattern=paste0("accuracy.+vs_",ref,".tab"))
  evaluation <- metadata[,1:4]
  i <- 0
  recovercol <- 0
  for (dir in dirs)  {
    evaluation <- cbind(evaluation, read.table(dir, header=T, sep="\t", as.is=T))
    recovercol <- c(recovercol,seq(9,13)+9*i)
    i <- i+1
  }
  recovercol <- recovercol[-1]
  write.table(evaluation, paste0("genomes_vs_",ref, ".tab"), row.names=F, col.names=T, sep="\t", quote=F)
  g <- grep("\\(\\)", evaluation$IV_jobID_raw)
  if (length(g)!=0) {
    m <- apply(evaluation[-g,recovercol], 2, mean, na.rm=T)
    s <- apply(evaluation[-g,recovercol], 2, sd, na.rm=T)
  } else {
    m <- apply(evaluation[,recovercol], 2, mean, na.rm=T)
    s <- apply(evaluation[,recovercol], 2, sd, na.rm=T)
  }
  r <- grep("recall", names(m))
  p <- grep("precision", names(m))
  a <- grep("accuracy", names(m))
  f <- grep("score", names(m))
  mcc <- grep("mcc", names(m))
  summarytable <- cbind(m[r], m[p], m[a], m[f], m[mcc],s[r], s[p], s[a],s[f], s[mcc])
  colnames(summarytable)<- c("recall", "precision", "accuracy", "fscore", "mcc", "sd_recall", "sd_precision", "sd_accuracy", "sd_fscore", "sd_mcc")
  summarytable <- summarytable*100
  write.table(summarytable, paste0("summary_evaluationIslandPath_vs_", ref,".tab"), row.names=T, col.names=T, sep="\t", quote=F)
  return(evaluation)
}
