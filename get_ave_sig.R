### get_ave_sig.R
### Given a directory of .rnk files, creates average signature of all .rnk files by taking the average of each gene.
# |dir| is a directory of .rnk files
# |collapse| if set to True, will only return genes that are found in all files. If False, will include all genes. Genes that are 
# not found in a particular .rnk file will be indicated with NA. The average will be the average of all non-NA values. 
# |ranked| if set to True, will average genes by relative rank in their respective lists. Otherwise will average genes by their metric.
# OUTPUT: A data frame containing all metric/rank data for all .rnk files, with genes on rows and .rnk file names on columns. The rightmost column will have the average metric/rank. 
get_ave_sig <- function(dir, collapse = F, ranked = F) {
  files <- list.files(dir, pattern = '.rnk')
  names <- gsub('.rnk','',files)
  allgenes = list()
  for (i in 1:length(files)) {
    temp <- read.delim(paste0(dir,'/',files[i]), row.names = 1, stringsAsFactors = F)
    temp <- temp[order(-temp),,drop = F]
    allgenes[[i]] <- temp
  }
  if (collapse == T) {
    common <- Reduce(intersect,lapply(allgenes, function(x) rownames(x)))
    allgenes.common <- data.frame(row.names = common)
    for (i in 1:length(allgenes)) {
      temp <- allgenes[[i]]
      temp <- temp[rownames(temp) %in% common,,drop = F]
      if (ranked == T) temp[,1] <- 1:nrow(temp)
      temp <- temp[common,,drop = F]
      allgenes.common <- cbind(allgenes.common, temp)
    }
    allgenes.common <- data.frame(allgenes.common, ave = rowMeans(allgenes.common))
  }
  if (collapse == F) {
    common <- unique(unlist(lapply(allgenes, function(x) rownames(x))))
    allgenes.common <- data.frame(row.names = common)
    for (i in 1:length(allgenes)) {
      temp <- allgenes[[i]]
      dif <- common[!common %in% rownames(temp)]
      if (ranked == T) temp[,1] <- 1:nrow(temp)
      temp2 <- data.frame(rep(NA, length(dif)), row.names = dif)
      colnames(temp2) <- colnames(temp)
      temp <- rbind(temp, temp2)
      temp <- temp[common,,drop = F]
      allgenes.common <- cbind(allgenes.common, temp)
    }    
    allgenes.common <- data.frame(allgenes.common, ave = rowMeans(allgenes.common, na.rm = T))
  }
  if (ranked == T) allgenes.common <- allgenes.common[order(allgenes.common$ave),]
  else if (ranked == F) allgenes.common <- allgenes.common[order(-allgenes.common$ave),]
  colnames(allgenes.common) <- c(names,'ave')
  return (allgenes.common)
}