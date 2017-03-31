# Author: Claire Bertelli (claire.bertelli@gmail.com)
# Last maintained 01.04.2017
# Formatting of results from GI prediction software into a 
###################
# Requires the sup tables 2 and 4 of Morgan Langille's BMC bioinformatics paper, 
# formatted to one file per genome in a three column format (name start stop)
#     

# load useful libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)

# Functions to parse the output from the various software - on file per genome -
# and format it to a three column based format - no header
##########
parse_GCProfile <- function() {
  
}

parse_MSGIP <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep("Putative Genomic Island", lines)
  if (length(g>=1)) {
    gis <- paste0("MSGIP_", seq(1, length(g)))
    start <- as.numeric(gsub("Putative Genomic Island: (\\d+\\.\\d+) - \\d+\\.\\d+Mb.", "\\1", lines[g]))*1000000
    end <- as.numeric(gsub("Putative Genomic Island: \\d+\\.\\d+ - (\\d+\\.\\d+)Mb.", "\\1", lines[g]))*1000000
    write.table(cbind(gis, start, end), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  } else {
    writeLines("", outputfile)
  }
}

parse_MTGIpick <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep("genomic_island", lines)
  if (length(g>=1)) {
    gis <- paste0("MTGIpick_", seq(1, length(g)))
    start <- as.numeric(gsub(".+genomic_island\\t(\\d+)\\t\\d+\\t.+", "\\1", lines[g]))
    end <- as.numeric(gsub(".+genomic_island\\t\\d+\\t(\\d+)\\t.+", "\\1", lines[g]))
    write.table(cbind(gis, start, end), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  } else {
    writeLines("", outputfile)
  }
}

parse_ZislandExplorer <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep("GI", lines)
  if (length(g>=1)) {
    gis <- paste0("ZislandExplorer_", seq(1, length(g)))
    start <- as.numeric(gsub("GI\\d+\\s+\\t(\\d+)\\.\\.\\d+\\s+\\t.+", "\\1", lines[g]))
    end <- as.numeric(gsub("GI\\d+\\s+\\t\\d+\\.\\.(\\d+)\\s+\\t.+", "\\1", lines[g]))
    write.table(cbind(gis, start, end), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  } else {
    writeLines("", outputfile)
  }
}

parse_SIGIHMM <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep("Colombo\tPUTAL", lines)
  nogenestart <- grep("PUTAL\t1\t0", lines[g])
  if (length(nogenestart)>=1) g <- g[-nogenestart]
  if (length(g>=1)) {
    start <- as.numeric(gsub(".+PUTAL\\t(\\d+)\\t\\d+\\t.+", "\\1", lines[g]))
    end <- as.numeric(gsub(".+PUTAL\\t\\d+\\t(\\d+)\\t.+", "\\1", lines[g]))
    genetab <- as.data.frame(cbind(start, end))
    difference = diff(g)
    rlediff = rle(difference)
    indexMultigenes <- which(rlediff$lengths>=1 & rlediff$values==1)
    gi <- t(data.frame(lapply(indexMultigenes, getSIGIgi, rlediff, g, genetab)))
    name <- paste0("SIGI_", seq(1, nrow(gi)))
    write.table(cbind(name, gi), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  } else {
    writeLines("", outputfile)
  }
}
getSIGIgi <- function(index, rlediff, w, genetab) {
  if (index==1) {
    startindex <- index
    endindex <- rlediff$length[index]
  } else {
    startindex <- sum(rlediff$lengths[1:(index-1)])+1
    endindex <- sum(rlediff$lengths[1:index])+1
  }
  gi <- c(genetab$start[startindex], genetab$end[endindex])
  return(gi)
}

# NEEDS TO BE IMPLEMENTED
parse_PredictBias <- function(inputfile, outputfile) {

}


#########
# little function to go from a two column of GI pred (start stop) to a three 
# column format (name start stop)
two2threeCol <- function(folder) {
  dir <- dir(folder, pattern=".txt", full.names = T)
  for (file in dir) { 
    tab <- try(read.table(file, header=F, sep="\t", as.is=T))
    if(inherits(tab, "try-error")) {
      tab<-''
    }else{
      write.table(cbind(paste0("GI_", seq(1:nrow(tab))), tab), gsub(".txt","_gis.txt",file), row.names=F, col.names=F, sep="\t", quote=F)
    }
  }
}
two2threeCol("ipath_dimob_v0.3")