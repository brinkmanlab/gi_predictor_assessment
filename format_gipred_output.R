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
    start <- as.numeric(gsub("GI\\d+\\t(\\d+)\\.\\.\\d+\\t.+", "\\1", lines[g]))
    end <- as.numeric(gsub("GI\\d+\\t\\d+\\.\\.(\\d+)\\t.+", "\\1", lines[g]))
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
    indices <- as.data.frame(cbind(start, end))[order(start),]
    i <- 1
    while(i < (nrow(indices))) {
      if(indices[i,2] >= (indices[(i+1),1]-1)) {
        # The ith entry in indices overlaps with the next entry
        if(indices[(i+1),2] > indices[i,2]) {
          indices[i,2] <- indices[(i+1),2] 
        }
        indices <- indices[-(i+1),,drop=FALSE] # drop=FALSE is necessary for nrow() to function when indices becomes 1 row
      }else{
        i <- i+1
      }
    }
    gis <- paste0("Centroid_", seq(1, nrow(indices)))
    options(scipen=999)
    write.table(cbind(gis, indices), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
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
        steps <- steps[-(i+1),,drop=FALSE] # drop=FALSE is necessary for nrow() to function when steps becomes 1 row
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

parse_SigHunt <- function(inputfile, outputfile, cutoff=5) {
  DIAS=read.table(inputfile,header=TRUE)
  DIAS$pass <- (DIAS$Hurricane >= cutoff)
  DIAS$gi <- 0
  # following for loop adapted from code in SigHuntTransform.r (https://www.iba.muni.cz/index-en.php?pg=research--data-analysis-tools--sighunt)
  # for each window that passes the cutoff, the adjacent two windows on either side are also considered part of the genomic island
  for (i in -2:2){
    move=(which(DIAS$pass==TRUE)+i)
    move=move[move>0 & move<=nrow(DIAS)]
    if(length(move)) {
      DIAS[move,]$gi <- 1
    }
  }
  DIAS <- subset(DIAS, DIAS$gi==1)
  if(nrow(DIAS)>0) {
    start <- as.numeric(gsub(">Region_\\d+_(\\d+)_(\\d+)", "\\1", row.names(DIAS)))
    end <- as.numeric(gsub(">Region_\\d+_(\\d+)_(\\d+)", "\\2", row.names(DIAS)))
    indices <- as.data.frame(cbind(start, end))
    i <- 1
    while(i < (nrow(indices))) {
      if(indices[i,2] >= (indices[(i+1),1])) {
        # The ith entry in indices overlaps with the next entry
        indices[i,2] <- indices[(i+1),2]
        indices <- indices[-(i+1),,drop=FALSE] # drop=FALSE is necessary for nrow() to function when indices becomes 1 row
      }else{
        i <- i+1
      }
    }
    gis <- paste0("SigHunt_", seq(1, nrow(indices)))
    options(scipen=999)
    write.table(cbind(gis, indices), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("",outputfile)
  }
}

parse_GCProfile <- function(inputfile, outputfile) {
  GC <- read.table(inputfile)
  GC <- data.frame(start = GC[seq(1,nrow(GC),2),1], end = GC[seq(2,nrow(GC),2),1], content = GC[seq(2,nrow(GC),2),2])
  GC$length <- GC$end - GC$start + 1
  # Genomic Islands are beyond weighted mean +/- 2 * weighted SD
  wtmean <- weighted.mean(GC$content, GC$length)
  wtsd <- sqrt(sum(GC$length * (GC$content - wtmean)^2)/sum(GC$length))
  cutoff <- wtsd * 2
  GC$GI <- as.factor(ifelse(abs(GC$content - wtmean) >= cutoff, 1, 0))
  GC <- subset(GC, GC$GI == 1)[,c(1,2)]
  if(nrow(GC)>0) {
    i <- 1
    while(i < (nrow(GC))) {
      if(GC[i,]$end == (GC[(i+1),]$start - 1)) {
        # The ith entry in GC is next to the next entry
        GC[i,]$end <- GC[(i+1),]$end
        GC <- GC[-(i+1),,drop=FALSE] # drop=FALSE is necessary for nrow() to function when indices becomes 1 row
      }else{
        i <- i+1
      }
    }
    gis <- paste0("GC-Profile_", seq(1, nrow(GC)))
    write.table(cbind(gis,GC), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("",outputfile)
  }
}

parse_WnSVM <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep("#HGT\\d+#", lines)
  if(length(g)) {
    indices <- gsub(".+\\t(\\d+)\\.\\.(\\d+)\\t.+","\\1\t\\2", lines[g])
    zeros <- grep("\\d+\\t0", indices)
    if(length(zeros)) {
      indices <- indices[-zeros]
    }
    labels <- paste0("Wn-SVM_", seq(1, length(indices)))
    write.table(cbind(labels, indices), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("", outputfile)
  }
}

parse_GI_SVM <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  indices <- gsub("(\\d+\\t\\d+)\\t0\\.\\d+", "\\1", lines)
  if(length(indices)) {
    labels <- paste0("GI_SVM_",seq(1, length(indices)))
    write.table(cbind(labels, indices), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("", outputfile)
  }
}

parse_Sigi_CRF <- function(inputfile, outputfile) {
  lines <- readLines(inputfile)
  g <- grep(".+PUTAL\\t(\\d+\\t\\d+)\\t.+", lines)
  if(length(g)) {
    indices <- gsub(".+PUTAL\\t(\\d+\\t\\d+)\\t.+", "\\1", lines[g])
    labels <- paste0("Sigi_CRF_", seq(1, length(g)))
    write.table(cbind(labels, indices), outputfile, row.names=F, col.names=F, sep="\t", quote=F)
  }else{
    writeLines("", outputfile)
  }
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
