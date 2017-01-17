### clean_data.R 
# Various functions to clean up data

### Fixes names of genes that are automatically converted to dates by Excel
# |names| is a vector of gene names
fix_date_gene <- function(names) {
  dec.nums <- gsub('-Dec|-DEC','',names[grepl('-Dec|-DEC',names)])
  sep.nums <- gsub('-Sep|-SEP','',names[grepl('-Sep|-SEP',names)])
  mar.nums <- gsub('-Mar|-MAR','',names[grepl('-Mar|-MAR',names)])
  names[grepl('-Dec|-DEC',names)] <- paste0('DEC',dec.nums)
  names[grepl('-Sep|-SEP',names)] <- paste0('SEPT',sep.nums)
  names[grepl('-Mar|-MAR',names)] <- paste0('MARCH',mar.nums)
  return(names)
}

### Takes the log of gene counts, removes genes with no variance between counts
# |data| is a matrix of raw or normalized gene counts with genes as the rows and samples as the columns
clean_dat <- function(data){
  novargenes < -which(apply(data,1,var)==0)
  if (length(novargenes) > 0) {
    data <- data[-novargenes,]
  }
  data <- na.omit(data)
  data <- log2(data+1)
  data <- data.matrix(data)
  return (data)
}

### Performs upper quartile normalization on raw gene counts
# |data| is a matrix of raw gene counts with genes as the rows and samples as the columns
quant_norm <- function(data){
  temp <- data
  temp <- as.matrix(temp)
  temp <- temp[rowSums(temp) > 0,]
  quants <- NULL
  for (i in 1:ncol(temp)) {
    column <- temp[,i]
    column <- column[!column == 0]
    f75 <- quantile(column, p = 0.75)
    quants <- c(quants, f75)
  }
  temp <- sweep(temp,2,quants,'/') * 1000
  return(temp)
}

### In instances where genes are ranked by some metric, randomizes the order of genes that are tied (oftentimes the tied genes are sorted alphabetically, which will throw off results)
# |data| is a data frame with two columns. The left column contains gene names and the right column contains the metric. The rows should be sorted by the metric.
randomize_ties <- function(data) {
  dup <- data[duplicated(data[,2]) | duplicated(data[,2], fromLast = T),]
  facts <- as.factor(dup[,2])
  for (level in levels(facts)) {
    indices <- which(data[,2] == level)
    range <- data[indices,]
    range <- range[sample(nrow(range)),]
    data <- rbind(data[1:(min(indices)-1),], range, data[(max(indices) + 1):nrow(data),])
  }
  return(data)
}