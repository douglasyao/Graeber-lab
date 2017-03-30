### Calculate p-value from generalized linear model for each row of matrix
# |data| is a matrix with genes as rows and samples as columns
# |model| is the model matrix
# |contrast| is the contrast
# |output| is the name of the file to be outputted. If NULL, will return pvalues as a data frame.
log_pvalues_from_glm <- function(data, model, contrast, output = NULL) {
  pvalues = NULL
  data <- as.matrix(data)
  for (i in 1:nrow(data)) {
    gene <- data[i,]
    summ <- summary(lm(as.numeric(gene)~model))
    summ <- summ$coefficients[contrast,]
    
    if (is.na(summ[4])) {
      pvalue <- NA
    }
    else {
      pvalue <- -log(summ[4])
      if (summ[3] < 0) {
        pvalue = (pvalue * -1)
      }
    }
    pvalues <- c(pvalues,pvalue)
  }
  pvalues <- data.frame(gene = rownames(data), pvalues = pvalues)
  pvalues <- na.omit(pvalues)
  pvalues <- pvalues[order(-pvalues$pvalue),]
  if (!is.null(output)) {
    write.table(pvalues, file = output, quote = F, sep = '\t', row.names = F)
  }
  else {
    return (pvalues)
  }
}


### Performs a simple t-test for each row of matrix
# |data| is a matrix with genes as rows and samples as columns
# |group1| and |group2| can be a list of the column names in each group or a vector of bools
# |paired| is whether to do paired t-test or not
# |remove.zeroes| filters out any genes where the percent of samples with no reads is greater than a certain threshold in either group. You can set the threshold like 0.2, 0.3 etc. 
# Usually I just set remove.zeroes = T which filters out genes that have 0 reads for all samples in either group. 
# |output| is the name of the file to be outputted. If NULL, will return pvalues as a data frame.
log_pvalues_from_ttest <- function(data, group1, group2, paired = F, remove.zeroes = F, output = NULL) {
  pvalues = NULL
  data <- as.matrix(data)
  for (i in 1:nrow(data)) {
    first <- c(data[i,group1])
    second <- c(data[i,group2])
    if (remove.zeroes == T) {
      if (sum(first) == 0 | sum(second) == 0) {
        pvalues <- c(pvalues,NA)
        next
      }
    }
    
    if (is.numeric(remove.zeroes)) {
      num.zeroes1 <- sum(first == 0)/length(first)
      num.zeroes2 <- sum(second == 0)/length(second)
      if (num.zeroes1 >= remove.zeroes | num.zeroes2 >= remove.zeroes) {
        pvalues <- c(pvalues,NA)
        next
      }
    }
    
    if (paired == T) {
      tTest <- t.test(as.numeric(first),as.numeric(second),paired=TRUE)
    }
    else {
      tTest <- t.test(as.numeric(first),as.numeric(second))
    }
    if (is.na(tTest$statistic)) {
      pvalue <- NA
    }
    else {
      pvalue <- -log(tTest$p.value)
      if (tTest$statistic < 0) {
        pvalue = (pvalue * -1)
      }
    }
    
    pvalues <- c(pvalues,pvalue)
  }
  pvalues <- data.frame(gene = rownames(data), pvalues = pvalues)
  pvalues <- na.omit(pvalues)
  pvalues[pvalues$pvalue == Inf,2] <- 1000
  pvalues[pvalues$pvalue == -Inf,2] <- -1000
  pvalues <- pvalues[order(-pvalues$pvalue),]
  if (!is.null(output)) {
    write.table(pvalues, file = output, quote = F, sep = '\t', row.names = F)
  }
  else {
    return (pvalues)
  }
}

### Same as log_pvalues_from_ttest, but applies Benjimini-Hochberg correction to p-values
log_pvalues_from_ttest_adj <- function(data, group1, group2, paired = F, remove.zeroes = F) {
  pvalues = NULL
  stats = NULL
  for (i in 1:nrow(data)) {
    first <- c(data[i,group1])
    second <- c(data[i,group2])
    if (remove.zeroes == T) {
      if (sum(first) == 0 | sum(second) == 0) {
        pvalues <- c(pvalues,NA)
        next
      }
    }
    
    if (is.numeric(remove.zeroes)) {
      num.zeroes1 <- sum(first == 0)/length(first)
      num.zeroes2 <- sum(second == 0)/length(second)
      if (num.zeroes1 >= remove.zeroes | num.zeroes2 >= remove.zeroes) {
        pvalues <- c(pvalues,NA)
        next
      }
    }
    if (paired == T) {
      tTest <- t.test(as.numeric(first),as.numeric(second),paired=TRUE)
    }
    else {
      tTest <- t.test(as.numeric(first),as.numeric(second))
    }
    if (is.na(tTest$statistic)) {
      pvalue <- NA
      stat <- NA
    }
    else {
      pvalue <- tTest$p.value
      stat <- tTest$statistic
    }
    pvalues <- c(pvalues, pvalue)
    stats <- c(stats, stat)
  }
  pvalues <- p.adjust(pvalues, method = 'BH')
  log.pvalues <- NULL
  for (i in 1:length(pvalues)) {
    if (is.na(pvalues[i])) {log.pvalue <- NA}
    else {
      log.pvalue <- -log(pvalues[i])
      if (stats[i] < 0) {
        log.pvalue <- -log.pvalue
      }
    }
    log.pvalues <- c(log.pvalues, log.pvalue)
  }
  
  log.pvalues <- data.frame(gene = rownames(data), pvalues = log.pvalues)
  log.pvalues <- na.omit(log.pvalues)
  log.pvalues <- log.pvalues[order(-log.pvalues$pvalues),]
  return (log.pvalues)
}


### Takes output from DESeq2 package and generates table of log p-values
# |data| is the output matrix
# |adj| to specify whether to apply multiple hypothesis testing correction
# |output| is the name of the file to be outputted. If NULL, will return pvalues as a data frame.
log_pvalues_from_deseq <- function(data, adj = F, output = NULL) {
  pvalues = NULL
  for (i in 1:nrow(data)) {
    d <- data[i,]
    if (adj == T) {pvalue <- -log(d$padj)}
    else {pvalue <- -log(d$pvalue)}
    if (d$stat < 0) {
      pvalue = (pvalue * -1)
    }
    pvalues <- c(pvalues,pvalue)
  }
  pvalues <- data.frame(gene = rownames(data), pvalue = pvalues)
  pvalues[pvalues$pvalue == Inf,2] <- 1000
  pvalues[pvalues$pvalue == -Inf,2] <- -1000
  pvalues <- pvalues[order(-pvalues$pvalue),]
  if (!is.null(output)) {
    write.table(pvalues, file = output, quote = F, sep = '\t', row.names = F)
  }
  else {
    return (pvalues)
  }
}

### Converts vector of log p-values to multiple hypothesis testing corrected log p-values
# |data| is vector of log p-values
# |base| is what base the log p-values are taken to
log_pvalue_to_adj_log_pvalues <- function(data, base) {
  negs <- which(data < 0)
  data <- abs(data)
  data <- base^(-data)
  adj <- p.adjust(data, method = 'BH')
  adj <- -log(adj)
  adj[negs] <- -adj[negs]
  return (adj)
}
