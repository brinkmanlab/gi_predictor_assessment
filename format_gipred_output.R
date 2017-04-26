# Author: Claire Bertelli (claire.bertelli@gmail.com)
# Last maintained 01.04.2017
# Functions to parse the output from GI prediction software into a three column 
# format (name start stop) without header. Need to have one file per genome for
# the functions to assess prediction accuracy.
###################


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

parse_mjsd <- function(inputfile,outputfile) {
  lines <- tryCatch(
    {
      read.table(inputfile, sep=" ", colClasses = c("character","numeric"))[1:2]
    },
    error=function(e) {
      cat("",file=outputfile,append=TRUE)
    }
  )
  if(!is.null(lines)) {
    cutoff <- 0.99999
    # removes all segments that are not considered genomic islands
    subset(lines, lines$V2 >= cutoff) -> lines
    if(nrow(lines)){
      lines$V2 <- NULL
      lines$istart <- as.numeric(gsub("-[0-9]*:","",lines$V1))
      lines$iend <- as.numeric(gsub(".*-(.*):", "\\1",lines$V1))
      lines$V1 <- NULL
      # the following while loop merges consecutive islands
      i <- 1
      while(i < (nrow(lines))) {
        while(i < nrow(lines) && lines[i,2] == lines[i+1,1]) {
          lines[i,2] <- lines[i+1,2]
          lines <- lines[-(i+1),]
        }
        i <- i+1
      }
      name <- paste0("MJSD_", seq(1, nrow(lines)))
      options(scipen=999)
      write.table(cbind(name,lines), file=outputfile, sep="\t", row.names=F, col.names=F)
    }
  }
  cat("",file=outputfile,append=TRUE)
}

parse_alien_hunter <- function(inputfile,outputfile) {
  lines <- readLines(inputfile)
  g <- grep("misc_feature", lines, value=TRUE)
  if (length(g>=1)) {
    gis <- paste0("alien_hunter_", seq(1, length(g)))
    indices <- gsub("FT   misc_feature    (\\d+)\\.\\.(\\d+)","\\1\t\\2", g)
    write.table(cbind(gis, indices), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  } else {
    writeLines("", outputfile)
  }
}

parse_centroid <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep(">(\\d+)-(\\d+).*", lines)
  if(length(g)) {
    start <- as.numeric(gsub(">(\\d+)-\\d+.*", "\\1", lines[g]))
    end <- as.numeric(gsub(">\\d+-(\\d+).*", "\\1", lines[g]))
    indices <- cbind(start, end)
    i <- 1
    while(i < (nrow(indices))) {
      j <- i+1
      while(j <= nrow(indices)) {
        if((indices[i,1] <= indices[j,2]) && (indices[i,2] >= indices[j,1])) {
          # The ith and jth entries in indices overlap
          if(indices[j,1] < indices[i,1]) {
            indices[i,1] <- indices[j,1]
          }
          if(indices[j,2] > indices[i,2]) {
            indices[i,2] <- indices[j,2]
          }
          indices <- indices[-j,]
        }
        j <- j+1
      }
      i <- i+1
    }
    gis <- paste0("Centroid_", seq(1, nrow(indices)))
    indices <- as.data.frame(indices)
    options(scipen=999)
    write.table(cbind(gis, indices[order(indices$start), ]), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("", outputfile)
  }
}

parse_pai_ida <- function(inputfile, outputfile) {
  cutoff <- 3.9
  steps <- read.table(inputfile, sep="\t", colClasses = c("V1"="numeric", "V8"="numeric"))[,c(1,8)]
  subset(steps, steps$V8 >= cutoff)[,1] -> steps
  if(length(steps)) {
    steps <- cbind(steps,(steps+5000))
    i <- 1
    while(i < nrow(steps)) {
      if(steps[i,2] == steps[(i+1),1]) {
        steps[i,2] <- steps[(i+1),2]
        steps <- steps[-(i+1),,drop=FALSE] # drop=FALSE is necessary for nrow() to function with 1 row
      }else{
        i <- i+1
      }
    }
    gis <- paste0("PAI-IDA_", seq(1, nrow(steps)))
    options(scipen=999)
    write.table(cbind(gis, steps), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("", outputfile)
  }
}

# NEEDS TO BE IMPLEMENTED
parse_PredictBias <- function(inputfile, outputfile) {

}

parse_GCProfile <- function() {
  
}

#########
# little function to go from a two column of GI pred (start stop) to a three 
# column format (name start stop)
# outputfile is <prefix>_gis.txt
two2threeCol <- function(inputfile, outputfile) {
  tab <- try(read.table(inputfile, header=F, sep="\t", as.is=T))
  if(inherits(tab, "try-error")) {
    tab<-''
  }else{
    write.table(cbind(paste0("GI_", seq(1:nrow(tab))), tab), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }
}