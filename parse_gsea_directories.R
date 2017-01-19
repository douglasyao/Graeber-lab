### parse_gsea_directories.R
# All functions will parse/analyze directories containing output from Gene Set Enrichment Analysis (GSEA) software from the Broad Institute

library(dplyr)
library(OIdata)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(pheatmap)
library(edgeR)
library(rowr)

### Returns a list of gene sets ordered by average enrichment across multiple GSEA directories
# |dir| is a directory containing multiple subdirectories, each of which is the output from running GSEA on a single .rnk metric file.
# |by| is either 'NES' or 'rank'. 'NES' will sort gene sets by the average normalized enrichment scores across all subdirectories. 'rank' will sort gene sets by rank.
gsea_obtain_top_genesets_across_multiple_gseas <- function(dir, by) {
  gseas <- list.files(dir)
  tot = NULL
  for (i in 1:length(gseas)) {
    files <- list.files(paste0(dir,'/',gseas[i]))
    neg <- files[grepl('gsea_report_for_na_neg_\\d*.xls', files)]
    neg <- read.csv(paste0(dir, '/', gseas[i], '/', neg), sep = '\t', row.names = 1)
    pos <- files[grepl('gsea_report_for_na_pos_\\d*.xls', files)]
    pos <- read.csv(paste0(dir, '/', gseas[i], '/', pos), sep = '\t', row.names = 1)
    comb <- rbind(neg,pos)
    comb <- comb[order(-comb$NES),]
    if (by == 'rank') comb <- data.frame(rank = 1:nrow(comb), row.names = rownames(comb))
    else if (by == 'NES') comb <- data.frame(nes = comb$NES, row.names = rownames(comb))
    if (is.null(tot)) {tot <- comb}
    else {
      comb <- comb[rownames(tot),,drop = F]
      tot <- cbind(tot, comb)
    }
  }
  colnames(tot) <- gseas
  if (by == 'rank') tot <- tot[order(rowMeans(tot)),]
  else if (by == 'NES') tot <- tot[order(-rowMeans(tot)),]
  return(tot)
}

### NOTE: Probably better to use v2 of this function
### Given a directory of metric files, plots a heatmap of the genes according to their metric.
# |genes| is a vector of gene names to plot
# |directory| is a directory containing .rnk metric files
# |outfileprefix| is the name of the output file
# |collapse| if set to True, will collapse duplicate genes and return the gene with the higher absolute value metric while retaining the sign
# OUTPUT: A png file containing a heatmap with genes on the rows and files on the columns. Color of heatmap is determined by scaled metric value of gene. Relative rank of gene is also overlayed. 
gsea_lookup_gene <- function(genes, directory, outfileprefix, collapse = F) {
  require(plyr)
  require(heatmap)
  files <- list.files(directory, pattern = 'rnk')
  combined_rvalues <- data.frame(row.names = genes)
  combined_ranks <- data.frame(row.names = genes)
  names <- gsub('.rnk','',files)

  for (file in files) {
    data <- read.delim(paste0(directory, file), stringsAsFactors = F)
    if (collapse == T) {
      data <- aggregate(data[,2,drop = F], by = data[,1,drop = F], function(x) return (x[which.max(abs(x))]))
    }
    posdata <- data[data[,2] > 0,, drop = F]
    negdata <- data[data[,2] < 0,, drop = F]
    posdata <- posdata[order(-posdata[,2]),, drop = F]
    rownames(posdata) <- 1:nrow(posdata)
    negdata <- negdata[order(negdata[,2]),, drop = F]
    rownames(negdata) <- 1:nrow(negdata)
    temp.genes <- NULL
    temp.ranks <- NULL
    for (gene in genes) {
      temp <- data[data[,1] == gene,]
      rvalue <- temp[,2]
      if (length(rvalue) == 0) {
        rvalue <- 0
        rank <- NA
      }
      
      else if (rvalue == 0) {
        rvalue <- 0
        rank <- 9999
      }
      
      else if (rvalue > 0) {
        rank <- row.names(posdata[posdata[,1] == gene,])
      }
      else if (rvalue < 0) {
        rank <- row.names(negdata[negdata[,1] == gene,])
      }
      temp.genes <- c(temp.genes, rvalue)
      temp.ranks <- c(temp.ranks, rank)
    }
    combined_rvalues <- cbind(combined_rvalues, temp.genes)
    combined_ranks <- cbind(combined_ranks, temp.ranks)
  }
  colnames(combined_rvalues) <- names
  colnames(combined_ranks) <- names
  
  combined_rvalues <- scale(combined_rvalues, center = F)

  maxn <- max(combined_rvalues)
  minn <- min(combined_rvalues)
  cols <- colorRampPalette(c("blue","white","red"))(100)
  if (abs(minn) > abs(maxn)) {
    proportion <- abs(round(maxn/minn * 50))
    cols <- cols[1:(50+proportion)]
  } else if (abs(maxn) > abs(minn)) {
    proportion <- abs(round(minn/maxn * 50))
    cols <- cols[(50-proportion):100]
  }
  
  data.plog = combined_rvalues
  data.plog <- data.plog[order(rowMeans(data.plog)),]
  combined_ranks <- combined_ranks[rownames(data.plog),]
  params <- list(data.plog,annotation=NA,color = cols,border_color="black",cellwidth=7,cellheight=4,fontsize=4, fontsize_number = 2.5, number_color = 'black', show_rownames = T, show_colnames = T, display_numbers = combined_ranks)
  do.call("pheatmap",c(params,cluster_rows = F, cluster_cols = F, filename = paste0(outfileprefix, '.png'),dpi=300))
  
  dev.off()
}


### Better version of gsea_lookup_gene
### Given a directory of metric files, plots a heatmap of the genes according to their metric.
# |genes| is a vector of gene names to plot
# |directory| is a directory containing .rnk metric files
# |outfileprefix| is the name of the output file
# |collapse| if set to True, will collapse duplicate genes and return the gene with the higher absolute value metric while retaining the sign
# |cluster| if set to True, will apply hierarchical clustering to both rows and columns of heatmap
# |transpose| if set to True, the output file with have genes on the columns and files on the rows
# |names| is a vector of names to substitute for the file names in the output heatmap. If set to NULL, the file names will be used.
# OUTPUT: A png file containing a heatmap with genes on the rows and files on the columns (opposite if transpose is set to True). 
# Color of heatmap is determined by scaled metric value taken to the fourth power (janky, but it works). Relative rank of gene and metric value are also overlayed. 
gsea_lookup_gene_v2 <- function(genes, directory, outfileprefix, collapse = F, cluster = F, transpose = F, names = NULL) {
  require(plyr)
  require(pheatmap)
  files <- list.files(directory, pattern = 'rnk')
  combined_pvalues <- data.frame(row.names = genes)
  combined_pvalues_adj <- data.frame(row.names = genes)
  combined_ranks <- data.frame(row.names = genes)
  names.files <- gsub('.rnk','',files)
  
  for (file in files) {
    data <- read.delim(paste0(directory, file), stringsAsFactors = F)
    if (collapse == T) {
      data <- aggregate(data[,2,drop = F], by = data[,1,drop = F], function(x) return (x[which.max(abs(x))]))
    }
    posdata <- data[data[,2] >= 0,, drop = F]
    negdata <- data[data[,2] < 0,, drop = F]
    posdata <- posdata[order(-posdata[,2]),, drop = F]
    rownames(posdata) <- 1:nrow(posdata)
    negdata <- negdata[order(negdata[,2]),, drop = F]
    rownames(negdata) <- 1:nrow(negdata)

    temp.pvalues <- NULL
    temp.pvalues.adj <- NULL
    temp.ranks <- NULL
    for (gene in genes) {
      temp <- data[data[,1] == gene,]
      pvalue <- temp[,2]
      if (length(pvalue) == 0) {
        pvalue <- NA
        pvalue.adj = 0
        rank <- NA
      }
      
      else if (pvalue == 0) {
        pvalue.adj = 0
        rank <- row.names(posdata[posdata[,1] == gene,])
      }
      
      else if (pvalue > 0) {
        rank <- row.names(posdata[posdata[,1] == gene,])
        pvalue.adj = ((-1 * as.numeric(rank) + nrow(posdata)) / nrow(posdata)) ^ 4
      }
      else if (pvalue < 0) {
        rank <- row.names(negdata[negdata[,1] == gene,])
        pvalue.adj = -1*(((as.numeric(rank) - nrow(negdata)) / nrow(negdata)) ^ 4)
      }
      temp.pvalues <- c(temp.pvalues, pvalue)
      temp.pvalues.adj <- c(temp.pvalues.adj, pvalue.adj)
      temp.ranks <- c(temp.ranks, rank)
    }
    combined_pvalues <- cbind(combined_pvalues, temp.pvalues)
    combined_pvalues_adj <- cbind(combined_pvalues_adj, temp.pvalues.adj)
    combined_ranks <- cbind(combined_ranks, temp.ranks)
  }
  colnames(combined_pvalues) <- names.files
  colnames(combined_pvalues_adj) <- names.files
  colnames(combined_ranks) <- names.files
  
  maxn <- max(combined_pvalues_adj)
  minn <- min(combined_pvalues_adj)
  cols <- colorRampPalette(c("blue","white","red"))(100)
  if (abs(minn) > abs(maxn)) {
    proportion <- abs(round(maxn/minn * 50))
    cols <- cols[1:(50+proportion)]
  } else if (abs(maxn) > abs(minn)) {
    proportion <- abs(round(minn/maxn * 50))
    cols <- cols[(50-proportion):100]
  }
  
  data.plog = combined_pvalues_adj
  data.plog <- data.plog[order(rowMeans(data.plog)),]
  combined_ranks <- combined_ranks[rownames(data.plog),]
  combined_pvalues <- combined_pvalues[rownames(data.plog),]
  combined = data.frame(row.names = rownames(combined_ranks))
  for (i in 1:ncol(combined_ranks)) {
    temp <- paste0(combined_ranks[,i],'\n',round(combined_pvalues[,i],digits = 2))
    combined <- cbind(combined,temp)
  }
  colnames(combined) <- colnames(combined_ranks)
  if (!is.null(names)) colnames(data.plog) <- names
  if (transpose == T) {
    combined <- t(combined)
    data.plog <- t(data.plog)
  }

  params <- list(data.plog,annotation=NA,color = cols,border_color="black",cellwidth=8,cellheight=8,fontsize=7, fontsize_number = 2.5, number_color = 'black', show_rownames = T, show_colnames = T, display_numbers = combined)
  
  callback = function(hc, data.plog) {
    colavg <- colSums(abs(data.plog)) / ncol(data.plog)
    rowavg <- rowSums(abs(data.plog)) / nrow(data.plog)
    dend = reorder(as.dendrogram(hc), wts = colavg)
    dend = reorder(as.dendrogram(dend), wts = rowavg)
    as.hclust(dend)
  }
  
  if (cluster == T) {
    do.call("pheatmap",c(params,clustering_method="average", clustering_callback = callback, filename = paste0(outfileprefix, '_euclid.png'),dpi=300))
    do.call("pheatmap",c(params,clustering_method="average",clustering_callback = callback, filename = paste0(outfileprefix, '_pearson.png'),clustering_distance_rows="correlation",clustering_distance_cols="correlation",dpi=300))
  }
  else {
    do.call("pheatmap",c(params,cluster_rows = F, cluster_cols = F, filename = paste0(outfileprefix, '.png'),dpi=300))
  }
  dev.off()
}

### Plots a heatmap of enrichment scores of gene sets from multiple GSEA directories
### calls function get_nes_ranks
# |genesets| is a vector of names of gene sets
# |outfileprefix| is the name of the output file
# |dir| is a directory containing multiple subdirectories, each of which is the output from running GSEA on a single .rnk metric file.
# OUTPUT: A png file containing a heatmap of the normalized enrichment scores (NESs) of the provided gene sets. Gene sets are on the rows, and each GSEA directory is a column. 
# The relative rank and NES of each gene set is also overlayed.
gsea_heatmap <- function(genesets, outfileprefix, dir) {
  require(pheatmap)
  require(pvclust)
  
  results <- get_nes_ranks(dir, genesets)
  data.p <- results[[1]]
  rankdata <- results[[2]]
  nesdata <- results[[3]]
  
  maxn <- max(data.p)
  minn <- min(data.p)
  cols <- colorRampPalette(c("blue","white","red"))(100)
  if (abs(minn) > abs(maxn)) {
    proportion <- abs(round(maxn/minn * 50))
    cols <- cols[1:(50+proportion)]
  } else if (abs(maxn) > abs(minn)) {
    proportion <- abs(round(minn/maxn * 50))
    cols <- cols[(50-proportion):100]
  }
  
  data.plog = data.p
  data.plog <- data.plog[order(rowMeans(data.plog)),,drop = F]
  rankdata <- rankdata[rownames(data.plog),,drop = F]
  nesdata <- nesdata[rownames(data.plog),, drop = F]
  
  combined.ranknes <- data.frame(row.names = rownames(rankdata))
  for (i in 1:ncol(rankdata)) {
    temp <- paste0(rankdata[,i],'\n',round(nesdata[,i],digits = 2))
    combined.ranknes <- cbind(combined.ranknes,temp)
  }
  colnames(combined.ranknes) <- colnames(rankdata)
  
  params <- list(data.plog,annotation=NA,color = cols,border_color="black",cellwidth=7,cellheight=7,fontsize=4, fontsize_number = 2.5, number_color = 'black', show_rownames = T, show_colnames = T, display_numbers = combined.ranknes)
  do.call("pheatmap",c(params,cluster_rows = F, cluster_cols = F, filename = paste0(outfileprefix, '.png'),dpi=300))
  dev.off()
}


### Returns enrichment scores and ranks of gene sets from multiple GSEA directories
# |dir| is a directory containing multiple subdirectories, each of which is the output from running GSEA on a single .rnk metric file.
# |genesets| is a vector of names of gene sets
# OUTPUT: A list containing 3 matrices. The first matrix has modified normalized enrichment scores (NES) of the gene sets. Modified NESs are truncated so that any NES with an absolute value below 1 is set to 0. 
# Otherwise, one is substracted from positive NESs and one is added to negative NESs. Gene sets are on the rows, and each GSEA directory is one of the columns. 
# The second data frame has the relative rank of each gene sets. The third matrix has the unmodified NESs of the gene sets. 
get_nes_ranks <- function(dir, genesets) {
  files <- list.files(dir)
  nesdata.adj <- data.frame(row.names = genesets)
  rankdata <- data.frame(row.names = genesets)
  nesdata <- data.frame(row.names = genesets)
  
  for (j in 1:length(files)) {
    subfile_list <- list.files(paste0(dir, '/', files[j]))
    gsea_neg <- subfile_list[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list)]
    gsea_neg <- read.csv(paste0(dir, '/', files[j], '/', gsea_neg), sep = '\t', row.names = 1)
    gsea_pos <- subfile_list[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list)]
    gsea_pos <- read.csv(paste0(dir, '/', files[j], '/', gsea_pos), sep = '\t', row.names = 1)
    
    gsea.neg.adj <- gsea_neg
    gsea.pos.adj <- gsea_pos
    
    gsea.neg.adj$NES <- gsea.neg.adj$NES + 1
    gsea.neg.adj$NES[gsea.neg.adj$NES > 0] <- 0
    gsea.neg.adj$RANK <- 1:nrow(gsea.neg.adj)
    
    gsea.pos.adj$NES <- gsea.pos.adj$NES - 1
    gsea.pos.adj$NES[gsea.pos.adj$NES < 0] <- 0
    gsea.pos.adj$RANK <- 1:nrow(gsea.pos.adj)
    
    combined.adj <- rbind(gsea.neg.adj, gsea.pos.adj)
    combined <- rbind(gsea_neg, gsea_pos)
    gsets.adj <- combined.adj[genesets,]
    gsets <- combined[genesets,]
    granks <- data.frame(gsets.adj$RANK, row.names = rownames(gsets.adj))
    gnes.adj <- data.frame(gsets.adj$NES, row.names = rownames(gsets.adj))
    gnes <- data.frame(gsets$NES, row.names = rownames(gsets))
    colnames(granks) <- files[j]
    colnames(gnes) <- files[j]
    colnames(gnes.adj) <- files[j]
    
    nesdata.adj <- cbind(nesdata.adj, gnes.adj)
    nesdata <- cbind(nesdata, gnes)
    rankdata <- cbind(rankdata, granks)
  }
  
  nesdata.adj[is.na(nesdata.adj)] <- 0
  nesdata.adj <- as.matrix(nesdata.adj)
  rankdata <- as.matrix(rankdata)
  nesdata <- as.matrix(nesdata)
  return (list(nesdata.adj, rankdata, nesdata))
}


### Renames the default output directories from GSEA to the name of the input file
# |dir| is a directory containing multiple subdirectories, each of which is the output from running GSEA on a single .rnk metric file.
rename_dir_gsea <- function(dir) {
  files <- list.files(dir, pattern = 'my_analysis')
  for (file in files) {
    html <- paste(readLines(paste0(dir,file,'/index.html')), collapse = '\n')
    name <- regmatches(html, regexpr('GSEA Report for Dataset .*</font></h3>\n</div><div><h4>', html))
    name <- gsub('GSEA Report for Dataset ', '', name)
    name <- gsub('</font></h3>\n</div><div><h4>', '', name)
    file.rename(paste0(dir,file), paste0(dir,name))
  }
}


#### OTHER FUNCTIONS
list_genesets_keyword <- function(keyword, type) {
  if (type == 'all') {
    c5 <- list_genesets_keyword(keyword, 'c5')
    c2 <- list_genesets_keyword(keyword, 'c2')
    return (c(c5,c2))
  }
  
  else if (type == 'c5') {
    dir <- '/Volumes/data0/users/dyao/GSEA_c5/GSEA_cell_line_bkpt/AUTONOMIC_GANGLIA_cell_line_bkpt/'
  }
  
  else {
    dir <- '/Volumes/data0/users/dyao/GSEA_c2/cell_line_bkpt/AUTONOMIC_GANGLIA_cell_line_bkpt/'
  }
  subfile_list <- list.files(dir)
  gsea_neg <- subfile_list[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list)]
  gsea_neg <- read.csv(paste0(dir, gsea_neg), sep = '\t')
  gsea_pos <- subfile_list[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list)]
  gsea_pos <- read.csv(paste0(dir, gsea_pos), sep = '\t')
  combined <- rbind(gsea_neg, gsea_pos)
  names <- as.character(combined$NAME)
  
  temp <- names[!grepl('REACTOME',names)]
  temp <- temp[!grepl('KEGG',temp)]
  temp <- temp[!grepl('BIOCARTA',temp)]
  temp <- temp[!grepl('^PID',temp)]
  
  
  if (type == 'cpg') {
    names = temp
  }
  
  if (type == 'c2') {
    names <- names[!(names %in% temp)]
  }
  outputnames = NULL
  for (i in 1:length(keyword)) {
    temp <- names[grepl(keyword[i], names)]
    outputnames <- c(outputnames,temp)
  }
  return (unique(outputnames))
}

list_genes_geneset <- function(genesets, type) {

  no.col <- max(count.fields('/Volumes/data0/users/dyao/GSEA_c2/cell_line_bkpt/AUTONOMIC_GANGLIA_cell_line_bkpt/edb/gene_sets.gmt', sep = '\t'))
  genesets.c2 <- read.delim('/Volumes/data0/users/dyao/GSEA_c2/cell_line_bkpt/AUTONOMIC_GANGLIA_cell_line_bkpt/edb/gene_sets.gmt', header = F, stringsAsFactors = F, fill = T, col.names = 1:no.col)

  no.col <- max(count.fields('/Volumes/data0/users/dyao/GSEA_c5/GSEA_cell_line_bkpt/AUTONOMIC_GANGLIA_cell_line_bkpt/edb/gene_sets.gmt', sep = '\t'))
  genesets.c5 <- read.delim('/Volumes/data0/users/dyao/GSEA_c5/GSEA_cell_line_bkpt/AUTONOMIC_GANGLIA_cell_line_bkpt/edb/gene_sets.gmt', header = F, stringsAsFactors = F, fill = T, col.names = 1:no.col)

  allgenes <- NULL
  for (geneset in genesets) {
    g <- genesets.c2[genesets.c2$X1 == geneset,]
    if (nrow(g) == 0) {
      g <- genesets.c5[genesets.c5$X1 == geneset,]
    }
    g <- g[!g == '']
    genes <- c(g[3:length(g)])
    allgenes <- c(allgenes, genes)
  }
  allgenes <- unique(allgenes)
  return (allgenes)
}

### Comparing GSEAs
compare_gsea <- function(directory1, directory2, label1, label2, type) {
  file_list1 <- list.files(directory1)
  gsea_neg1 <- file_list1[grepl('gsea_report_for_na_neg_\\d*.xls', file_list1)]
  gsea_neg_file1 <- read.csv(paste0(directory1,'/',gsea_neg1), sep = '\t')
  gsea_pos1 <- file_list1[grepl('gsea_report_for_na_pos_\\d*.xls', file_list1)]
  gsea_pos_file1 <- read.csv(paste0(directory1,'/',gsea_pos1), sep = '\t')
  
  file_list2 <- list.files(directory2)
  gsea_neg2 <- file_list2[grepl('gsea_report_for_na_neg_\\d*.xls', file_list2)]
  gsea_neg_file2 <- read.csv(paste0(directory2,'/',gsea_neg2), sep = '\t')
  gsea_pos2 <- file_list2[grepl('gsea_report_for_na_pos_\\d*.xls', file_list2)]
  gsea_pos_file2 <- read.csv(paste0(directory2,'/',gsea_pos2), sep = '\t')
  
  both_neg <- intersect(gsea_neg_file1$NAME, gsea_neg_file2$NAME)
  neg_file1_intersect <- gsea_neg_file1[match(both_neg, gsea_neg_file1$NAME),]
  neg_file2_intersect <- gsea_neg_file2[match(both_neg, gsea_neg_file2$NAME),]
  both_neg_table <- data.frame(both_neg, neg_file1_intersect$NES, neg_file1_intersect$NOM.p.val, neg_file1_intersect$FDR.q.val, neg_file2_intersect$NES, neg_file2_intersect$NOM.p.val, neg_file2_intersect$FDR.q.val)
  both_neg_table <- both_neg_table[(both_neg_table$neg_file1_intersect.FDR.q.val < 0.25 & both_neg_table$neg_file2_intersect.FDR.q.val < 0.25),]
  colnames(both_neg_table) <- c('NAME', paste0(label1, '_NES'), paste0(label1, '_pval'), paste0(label1, '_FDR_qval'), paste0(label2, '_NES'), paste0(label2, '_pval'), paste0(label2, '_FDR_qval'))
  both_neg_table <- both_neg_table[order((both_neg_table[,4] + both_neg_table[,7])/2),]
  
  both_pos <- intersect(gsea_pos_file1$NAME, gsea_pos_file2$NAME)
  pos_file1_intersect <- gsea_pos_file1[match(both_pos, gsea_pos_file1$NAME),]
  pos_file2_intersect <- gsea_pos_file2[match(both_pos, gsea_pos_file2$NAME),]
  both_pos_table <- data.frame(both_pos, pos_file1_intersect$NES, pos_file1_intersect$NOM.p.val, pos_file1_intersect$FDR.q.val, pos_file2_intersect$NES, pos_file2_intersect$NOM.p.val, pos_file2_intersect$FDR.q.val)
  both_pos_table <- both_pos_table[(both_pos_table$pos_file1_intersect.FDR.q.val < 0.25 & both_pos_table$pos_file2_intersect.FDR.q.val < 0.25),]
  colnames(both_pos_table) <- c('NAME', paste0(label1, '_NES'), paste0(label1, '_pval'), paste0(label1, '_FDR_qval'), paste0(label2, '_NES'), paste0(label2, '_pval'), paste0(label2, '_FDR_qval'))
  both_pos_table <- both_pos_table[order((both_pos_table[,4] + both_pos_table[,7])/2),]
  
  pos_neg <- intersect(gsea_pos_file1$NAME, gsea_neg_file2$NAME)
  pos_file1_intersect <- gsea_pos_file1[match(pos_neg, gsea_pos_file1$NAME),]
  neg_file2_intersect <- gsea_neg_file2[match(pos_neg, gsea_neg_file2$NAME),]
  pos_neg_table <- data.frame(pos_neg, pos_file1_intersect$NES, pos_file1_intersect$NOM.p.val, pos_file1_intersect$FDR.q.val, neg_file2_intersect$NES, neg_file2_intersect$NOM.p.val, neg_file2_intersect$FDR.q.val)
  pos_neg_table <- pos_neg_table[(pos_neg_table$pos_file1_intersect.FDR.q.val < 0.25 & pos_neg_table$neg_file2_intersect.FDR.q.val < 0.25),]
  colnames(pos_neg_table) <- c('NAME', paste0(label1, '_NES'), paste0(label1, '_pval'), paste0(label1, '_FDR_qval'), paste0(label2, '_NES'), paste0(label2, '_pval'), paste0(label2, '_FDR_qval'))
  pos_neg_table <- pos_neg_table[order((pos_neg_table[,4] + pos_neg_table[,7])/2),]
  
  neg_pos <- intersect(gsea_neg_file1$NAME, gsea_pos_file2$NAME)
  neg_file1_intersect <- gsea_neg_file1[match(neg_pos, gsea_neg_file1$NAME),]
  pos_file2_intersect <- gsea_pos_file2[match(neg_pos, gsea_pos_file2$NAME),]
  neg_pos_table <- data.frame(neg_pos, neg_file1_intersect$NES, neg_file1_intersect$NOM.p.val, neg_file1_intersect$FDR.q.val, pos_file2_intersect$NES, pos_file2_intersect$NOM.p.val, pos_file2_intersect$FDR.q.val)
  neg_pos_table <- neg_pos_table[(neg_pos_table$neg_file1_intersect.FDR.q.val < 0.25 & neg_pos_table$pos_file2_intersect.FDR.q.val < 0.25),]
  colnames(neg_pos_table) <- c('NAME', paste0(label1, '_NES'), paste0(label1, '_pval'), paste0(label1, '_FDR_qval'), paste0(label2, '_NES'), paste0(label2, '_pval'), paste0(label2, '_FDR_qval'))
  neg_pos_table <- neg_pos_table[order((neg_pos_table[,4] + neg_pos_table[,7])/2),]
  
  if (type == 'all') {
    dir.create(paste0(label1, '_vs_', label2, '_GSEA_comparisons'))
    current_dir <- getwd()
    setwd(paste0(label1, '_vs_', label2, '_GSEA_comparisons'))
    write.table(gsea_pos_file1, file = paste0(label1, '_up.txt'), quote = F, sep = '\t', row.names = F)
    write.table(gsea_neg_file1, file = paste0(label1, '_down.txt'), quote = F, sep = '\t', row.names = F)
    write.table(gsea_pos_file2, file = paste0(label2, '_up.txt'), quote = F, sep = '\t', row.names = F)
    write.table(gsea_neg_file2, file = paste0(label2, '_down.txt'), quote = F, sep = '\t', row.names = F)
    write.table(both_pos_table, file = paste0(label1, '_', label2, '_both_up.txt'), quote = F, sep = '\t', row.names = F)
    write.table(both_neg_table, file = paste0(label1, '_', label2, '_both_down.txt'), quote = F, sep = '\t', row.names = F)
    write.table(pos_neg_table, file = paste0(label1, '_up_', label2, '_down.txt'), quote = F, sep = '\t', row.names = F)
    write.table(neg_pos_table, file = paste0(label1, '_down_', label2, '_up.txt'), quote = F, sep = '\t', row.names = F)
    setwd(current_dir)
  }
  
  if (type == 'simple') {
    total = NULL
    combined <- list(both_pos_table, both_neg_table, pos_neg_table, neg_pos_table)
    for (i in 1:4) {
      number <- nrow(combined[[i]])
      if (number == 0) {
        temp <- c(rep('',10))
        temp <- c(temp, '0 overlapping gene sets')
        total <- cbind(total, temp)
      }
      else if (number < 10) {
        temp <- as.character(combined[[i]]$NAME)[1:number]
        temp <- c(temp, rep('',10-number))
        temp <- c(temp, paste0(number, ' overlapping gene sets'))
        temp <- strtrim(temp, 40)
        total <- cbind(total, temp)
      }
      else if (number >= 10) {
        temp <- as.character(combined[[i]]$NAME[1:10])
        temp <- strtrim(temp, 40)
        temp <- c(temp, paste0(number, ' overlapping gene sets'))
        total <- cbind(total, temp)
      }
    }
    colnames(total) <- c('Both up', 'Both down', paste0('Up ',gsub('_.*','',label1),' down ',gsub('_.*','',label2)), paste0('Down ',gsub('_.*','',label1),' up ',gsub('_.*','',label2)))
    current_dir <- getwd()
    png(paste0(current_dir,'/',label1,'_vs_',label2,'_comparisons.png'), height = 250, width = 1250)
    grid.arrange(tableGrob(total))
    dev.off()
  }
}

compare_gsea_v2 <- function(directory1, directory2, label1, label2, type) {
  file_list1 <- list.files(directory1)
  gsea_neg1 <- file_list1[grepl('gsea_report_for_na_neg_\\d*.xls', file_list1)]
  gsea_neg_file1 <- read.csv(paste0(directory1,'/',gsea_neg1), sep = '\t')
  gsea_neg_file1 <- data.frame(gsea_neg_file1, RANK = rownames(gsea_neg_file1))
  gsea_pos1 <- file_list1[grepl('gsea_report_for_na_pos_\\d*.xls', file_list1)]
  gsea_pos_file1 <- read.csv(paste0(directory1,'/',gsea_pos1), sep = '\t')
  gsea_pos_file1 <- data.frame(gsea_pos_file1, RANK = rownames(gsea_pos_file1))
  
  file_list2 <- list.files(directory2)
  gsea_neg2 <- file_list2[grepl('gsea_report_for_na_neg_\\d*.xls', file_list2)]
  gsea_neg_file2 <- read.csv(paste0(directory2,'/',gsea_neg2), sep = '\t')
  gsea_neg_file2 <- data.frame(gsea_neg_file2, RANK = rownames(gsea_neg_file2))
  gsea_pos2 <- file_list2[grepl('gsea_report_for_na_pos_\\d*.xls', file_list2)]
  gsea_pos_file2 <- read.csv(paste0(directory2,'/',gsea_pos2), sep = '\t')
  gsea_pos_file2 <- data.frame(gsea_pos_file2, RANK = rownames(gsea_pos_file2))
  
  temp_neg_file1 <- gsea_neg_file1[!grepl('REACTOME',gsea_neg_file1$NAME),]
  temp_neg_file1 <- temp_neg_file1[!grepl('KEGG',temp_neg_file1$NAME),]
  temp_neg_file1 <- temp_neg_file1[!grepl('BIOCARTA',temp_neg_file1$NAME),]
  temp_pos_file1 <- gsea_pos_file1[!grepl('REACTOME',gsea_pos_file1$NAME),]
  temp_pos_file1 <- temp_pos_file1[!grepl('KEGG',temp_pos_file1$NAME),]
  temp_pos_file1 <- temp_pos_file1[!grepl('BIOCARTA',temp_pos_file1$NAME),]
  temp_neg_file2 <- gsea_neg_file2[!grepl('REACTOME',gsea_neg_file2$NAME),]
  temp_neg_file2 <- temp_neg_file2[!grepl('KEGG',temp_neg_file2$NAME),]
  temp_neg_file2 <- temp_neg_file2[!grepl('BIOCARTA',temp_neg_file2$NAME),]
  temp_pos_file2 <- gsea_pos_file2[!grepl('REACTOME',gsea_pos_file2$NAME),]
  temp_pos_file2 <- temp_pos_file2[!grepl('KEGG',temp_pos_file2$NAME),]
  temp_pos_file2 <- temp_pos_file2[!grepl('BIOCARTA',temp_pos_file2$NAME),]
  
  if (type == 'cpg') {
    gsea_neg_file1 <- temp_neg_file1
    gsea_neg_file2 <- temp_neg_file2
    gsea_pos_file1 <- temp_pos_file1
    gsea_pos_file2 <- temp_pos_file2
  }
  
  if (type == 'c2') {
    gsea_neg_file1 <- gsea_neg_file1[!(as.character(gsea_neg_file1$NAME) %in% as.character(temp_neg_file1$NAME)),]
    gsea_neg_file1$RANK <- 1:nrow(gsea_neg_file1)
    gsea_neg_file2 <- gsea_neg_file2[!(as.character(gsea_neg_file2$NAME) %in% as.character(temp_neg_file2$NAME)),]
    gsea_neg_file2$RANK <- 1:nrow(gsea_neg_file2)
    gsea_pos_file1 <- gsea_pos_file1[!(as.character(gsea_pos_file1$NAME) %in% as.character(temp_pos_file1$NAME)),]
    gsea_pos_file1$RANK <- 1:nrow(gsea_pos_file1)
    gsea_pos_file2 <- gsea_pos_file2[!(as.character(gsea_pos_file2$NAME) %in% as.character(temp_pos_file2$NAME)),]
    gsea_pos_file2$RANK <- 1:nrow(gsea_pos_file2)
  }
  
  both_neg <- intersect(gsea_neg_file1$NAME, gsea_neg_file2$NAME)
  neg_file1_intersect <- gsea_neg_file1[match(both_neg, gsea_neg_file1$NAME),]
  neg_file2_intersect <- gsea_neg_file2[match(both_neg, gsea_neg_file2$NAME),]
  both_neg_table <- data.frame(both_neg, neg_file1_intersect$NES, neg_file1_intersect$FDR.q.val, neg_file1_intersect$RANK, neg_file2_intersect$NES, neg_file2_intersect$FDR.q.val, neg_file2_intersect$RANK)
  both_neg_table <- both_neg_table[(both_neg_table$neg_file1_intersect.FDR.q.val < 0.25 & both_neg_table$neg_file2_intersect.FDR.q.val < 0.25),]
  colnames(both_neg_table) <- c('NAME', paste0(label1, '_NES'), paste0(label1, '_FDR_qval'), paste0(label1, '_RANK'), paste0(label2, '_NES'), paste0(label2, '_FDR_qval'), paste0(label2, '_RANK'))
  both_neg_table <- both_neg_table[order((both_neg_table[,2] + both_neg_table[,5])/2),]
  
  both_pos <- intersect(gsea_pos_file1$NAME, gsea_pos_file2$NAME)
  pos_file1_intersect <- gsea_pos_file1[match(both_pos, gsea_pos_file1$NAME),]
  pos_file2_intersect <- gsea_pos_file2[match(both_pos, gsea_pos_file2$NAME),]
  both_pos_table <- data.frame(both_pos, pos_file1_intersect$NES, pos_file1_intersect$FDR.q.val, pos_file1_intersect$RANK, pos_file2_intersect$NES, pos_file2_intersect$FDR.q.val, pos_file2_intersect$RANK)
  both_pos_table <- both_pos_table[(both_pos_table$pos_file1_intersect.FDR.q.val < 0.25 & both_pos_table$pos_file2_intersect.FDR.q.val < 0.25),]
  colnames(both_pos_table) <- c('NAME', paste0(label1, '_NES'), paste0(label1, '_FDR_qval'), paste0(label1, '_RANK'), paste0(label2, '_NES'), paste0(label2, '_FDR_qval'), paste0(label2, '_RANK'))
  both_pos_table <- both_pos_table[order(-(both_pos_table[,2] + both_pos_table[,5])/2),]
  
  total = NULL
  combined <- cbind(both_pos_table[1:100,], both_neg_table[1:100,])
  write.table(combined, paste0(label1,'_vs_',label2,'_comparisons_', type, '.txt'), row.names = F, quote = F, sep = '\t')
}

compare_gsea_all_cross_comparisons <- function(directory1, directory2, type) {
  file_list1 <- list.files(directory1)
  file_list2 <- list.files(directory2)
  current_dir <- getwd()
  dir.create(paste0(directory1,'_vs_',directory2))
  setwd(paste0(directory1,'_vs_',directory2))
  combinations <- expand.grid(file_list1, file_list2) 
  for (i in 1:nrow(combinations)) {
    compare_gsea(paste0('../',directory1,'/',combinations[i,1]), paste0('../', directory2,'/',combinations[i,2]), as.character(combinations[i,1]), as.character(combinations[i,2]), type = type)
  }
  setwd(current_dir)
}



### PARSING GSEA DIRECTORIES
gsea_obtain_top <- function(directory) {
  file_list <- list.files(directory)
  combined_neg = NULL
  combined_pos = NULL
  for (i in file_list) {
    subfile_list <- list.files(paste0(directory,'/',i))
    gsea_neg <- subfile_list[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list)]
    gsea_neg_file <- read.csv(paste0(directory,'/',i,'/',gsea_neg), sep = '\t')
    gsea_pos <- subfile_list[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list)]
    gsea_pos_file <- read.csv(paste0(directory,'/',i,'/',gsea_pos), sep = '\t')
    
    combined_neg <- cbind(combined_neg, as.character(gsea_neg_file$NAME[1:20]))
    combined_pos <- cbind(combined_pos, as.character(gsea_pos_file$NAME[1:20]))
  }
  
  colnames(combined_neg) <- file_list
  colnames(combined_pos) <- file_list
  
  write.table(combined_neg, filename = 'up_all.txt', quote = F, row.names = F, sep = '\t')
  write.table(combined_pos, filename = 'down_all.txt', quote = F, row.names = F, sep = '\t')
  
}

gsea_obtain_top_combined <- function(directory1, directory2, label1) {
  
  subfile_list1 <- list.files(directory1)
  gsea_neg1 <- subfile_list1[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list1)]
  gsea_neg1 <- read.csv(paste0(directory1,'/',gsea_neg1), sep = '\t')
  gsea_pos1 <- subfile_list1[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list1)]
  gsea_pos1 <- read.csv(paste0(directory1,'/',gsea_pos1), sep = '\t')
  
  subfile_list2 <- list.files(directory2)
  gsea_neg2 <- subfile_list2[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list2)]
  gsea_neg2 <- read.csv(paste0(directory2,'/',gsea_neg2), sep = '\t')
  gsea_pos2 <- subfile_list2[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list2)]
  gsea_pos2 <- read.csv(paste0(directory2,'/',gsea_pos2), sep = '\t')
  
  nums <- paste0(1:10, '. ')
  nums2 <- paste0(11:20, '. ')
  combined <- data.frame(strtrim(paste0(nums, as.character(gsea_pos1$NAME[1:10])), 36), strtrim(paste0(nums2, as.character(gsea_pos1$NAME[11:20])), 36), strtrim(paste0(nums, as.character(gsea_neg1$NAME[1:10])), 36), strtrim(paste0(nums2, as.character(gsea_neg1$NAME[11:20])), 36), strtrim(paste0(nums, as.character(gsea_pos2$NAME[1:10])), 36), strtrim(paste0(nums2, as.character(gsea_pos2$NAME[11:20])), 36), strtrim(paste0(nums, as.character(gsea_neg2$NAME[1:10])), 36), strtrim(paste0(nums2, as.character(gsea_neg2$NAME[11:20])), 36))
  colnames(combined) <- NULL
  rownames(combined) <- NULL
  png(paste0(label1,'_combined_top_20.png'), height = 250, width = 1875)
  t <- tableGrob(combined, rows = NULL, theme = ttheme_default(base_size = 10, base_colour = 'black'))
  t <- gtable_add_rows(t, unit(1, "line"), 0)
  t <- gtable_add_grob(t, list(textGrob(paste0('Top 20 Up ', label1, ' by breakpoints'), gp = gpar(fontsize = 14, fontface = 'bold')), textGrob(paste0('Top 20 Down ', label1,' by breakpoints'),gp = gpar(fontsize = 14, fontface = 'bold')), textGrob(paste0('Top 20 Up ', label1, ' by ICNA'),gp = gpar(fontsize = 14, fontface = 'bold')), textGrob(paste0('Top 20 Down ', label1, ' by ICNA'),gp = gpar(fontsize = 14, fontface = 'bold'))), t=1,b=1,l=c(1, 3, 5, 7), r=c(2, 4, 6, 8))
  seg1 <- segmentsGrob(x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(0,"npc"), y1 = unit(10,"npc"), gp = gpar(lwd = 3.0))
  t <- gtable_add_grob(t, seg1, t = 2, b = 11, l = 3, r = 3)
  seg1 <- segmentsGrob(x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(0,"npc"), y1 = unit(10,"npc"), gp = gpar(lwd = 3.0))
  t <- gtable_add_grob(t, seg1, t = 2, b = 11, l = 5, r = 5)
  seg1 <- segmentsGrob(x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(0,"npc"), y1 = unit(10,"npc"), gp = gpar(lwd = 3.0))
  t <- gtable_add_grob(t, seg1, t = 2, b = 11, l = 7, r = 7)
  grid.arrange(t)
  dev.off()
}

gsea_obtain_top_single <- function(directory1, label) {
  subfile_list1 <- list.files(directory1)
  gsea_neg1 <- subfile_list1[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list1)]
  gsea_neg1 <- read.csv(paste0(directory1,'/',gsea_neg1), sep = '\t')
  gsea_neg1 <- gsea_neg1[(gsea_neg1$FDR.q.val < 0.25),]
  gsea_neg1 <- gsea_neg1[order(gsea_neg1$FDR.q.val),]
  
  gsea_pos1 <- subfile_list1[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list1)]
  gsea_pos1 <- read.csv(paste0(directory1,'/',gsea_pos1), sep = '\t')
  gsea_pos1 <- gsea_pos1[(gsea_pos1$FDR.q.val < 0.25),]
  gsea_pos1 <- gsea_pos1[order(gsea_pos1$FDR.q.val),]
  
  if (nrow(gsea_pos1) == 0) {
    names1 <- rep('',10)
    scores1 <- rep('',10)
    names2 <- rep('',10)
    scores2 <- rep('',10)
  }
  
  else if (nrow(gsea_pos1) <= 10) {
    nums <- paste0(1:nrow(gsea_pos1), '. ')
    names1 <- strtrim(paste0(nums, as.character(gsea_pos1$NAME[1:nrow(gsea_pos1)])), 50)
    scores1 <- strtrim(as.character(gsea_pos1$FDR.q.val[1:nrow(gsea_pos1)]), 5)
    names1 <- c(names1, rep('',10-nrow(gsea_pos1)))
    scores1 <- c(scores1, rep('',10-nrow(gsea_pos1)))
    names2 <- rep('',10)
    scores2 <- rep('',10)
  }
  
  else if (nrow(gsea_pos1) <= 20 & nrow(gsea_pos1) > 10) {
    nums <- paste0(1:10, '. ')
    nums2 <- paste0(11:nrow(gsea_pos1), '. ')
    names1 <- strtrim(paste0(nums, as.character(gsea_pos1$NAME[1:10])), 50)
    scores1 <- strtrim(as.character(gsea_pos1$FDR.q.val[1:10]), 5)
    names2 <- strtrim(paste0(nums2, as.character(gsea_pos1$NAME[11:nrow(gsea_pos1)])), 50)
    scores2 <- strtrim(as.character(gsea_pos1$FDR.q.val[11:nrow(gsea_pos1)]), 5)
    names2 <- c(names2, rep('',20-nrow(gsea_pos1)))
    scores2 <- c(scores2, rep('',20-nrow(gsea_pos1)))
  }
  
  else if (nrow(gsea_pos1) > 20) {
    nums <- paste0(1:10, '. ')
    nums2 <- paste0(11:20, '. ')
    names1 <- strtrim(paste0(nums, as.character(gsea_pos1$NAME[1:10])), 50)
    scores1 <- strtrim(as.character(gsea_pos1$FDR.q.val[1:10]), 5)
    names2 <- strtrim(paste0(nums2, as.character(gsea_pos1$NAME[11:20])), 50)
    scores2 <- strtrim(as.character(gsea_pos1$FDR.q.val[11:20]), 5)
  }
  
  if (nrow(gsea_neg1) == 0) {
    names3 <- rep('',10)
    scores3 <- rep('',10)
    names4 <- rep('',10)
    scores4 <- rep('',10)
  }
  
  else if (nrow(gsea_neg1) <= 10) {
    nums <- paste0(1:nrow(gsea_neg1), '. ')
    names3 <- strtrim(paste0(nums, as.character(gsea_neg1$NAME[1:nrow(gsea_neg1)])), 50)
    scores3 <- strtrim(as.character(gsea_neg1$FDR.q.val[1:nrow(gsea_neg1)]), 5)
    names3 <- c(names3, rep('',10-nrow(gsea_neg1)))
    scores3 <- c(scores3, rep('',10-nrow(gsea_neg1)))
    names4 <- rep('',10)
    scores4 <- rep('',10)
  }
  
  else if (nrow(gsea_neg1) > 10 & nrow(gsea_neg1) <= 20) {
    nums <- paste0(1:10, '. ')
    nums2 <- paste0(11:nrow(gsea_neg1), '. ')
    names3 <- strtrim(paste0(nums, as.character(gsea_neg1$NAME[1:10])), 50)
    scores3 <- strtrim(as.character(gsea_neg1$FDR.q.val[1:10]), 5)
    names4 <- strtrim(paste0(nums2, as.character(gsea_neg1$NAME[11:nrow(gsea_neg1)])), 50)
    scores4 <- strtrim(as.character(gsea_neg1$FDR.q.val[11:nrow(gsea_neg1)]), 5)
    names4 <- c(names4, rep('',20-nrow(gsea_neg1)))
    scores4 <- c(scores4, rep('',20-nrow(gsea_neg1)))
  }
  
  else if (nrow(gsea_neg1) > 20) {
    nums <- paste0(1:10, '. ')
    nums2 <- paste0(11:20, '. ')
    names3 <- strtrim(paste0(nums, as.character(gsea_neg1$NAME[1:10])), 50)
    scores3 <- strtrim(as.character(gsea_neg1$FDR.q.val[1:10]), 5)
    names4 <- strtrim(paste0(nums2, as.character(gsea_neg1$NAME[11:20])), 50)
    scores4 <- strtrim(as.character(gsea_neg1$FDR.q.val[11:20]), 5)
  }
  
  combined <- data.frame(names1, scores1, names2, scores2, names3, scores3, names4, scores4)
  rownames(combined) <- NULL
  colnames(combined) <- c(rep(c('Name', 'FDR q-val'),4))
  png(paste0(label,'_combined_top_20.png'), height = 250, width = 1875)
  t <- tableGrob(combined, rows = NULL, theme = ttheme_default(base_size = 12, base_colour = 'black'))
  t <- gtable_add_rows(t, unit(1, "line"), 0)
  t <- gtable_add_grob(t, list(textGrob(paste0('Top Upregulated ', label), gp = gpar(fontsize = 14, fontface = 'bold')), textGrob(paste0('Top Downregulated ', label) ,gp = gpar(fontsize = 14, fontface = 'bold'))), t=1,b=1,l=c(1,5), r=c(4,8))
  seg1 <- segmentsGrob(x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(0,"npc"), y1 = unit(11,"npc"), gp = gpar(lwd = 2.0))
  t <- gtable_add_grob(t, seg1, t = 2, b = 12, l = 3, r = 3)
  seg1 <- segmentsGrob(x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(0,"npc"), y1 = unit(11,"npc"), gp = gpar(lwd = 2.0))
  t <- gtable_add_grob(t, seg1, t = 2, b = 12, l = 5, r = 5)
  seg1 <- segmentsGrob(x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(0,"npc"), y1 = unit(11,"npc"), gp = gpar(lwd = 2.0))
  t <- gtable_add_grob(t, seg1, t = 2, b = 12, l = 7, r = 7)
  grid.arrange(t)
  dev.off()
}

### Genes from GSEA
library(plyr)
library(rowr)
get_genes <- function(directory, outfile_prefix, num_genes) {
  file_list <- list.files(directory)
  gsea_neg <- file_list[grepl('gsea_report_for_na_neg_\\d*.xls', file_list)]
  gsea_neg <- read.csv(paste0(directory,'/',gsea_neg), sep = '\t')
  gsea_neg <- gsea_neg[order(gsea_neg$FDR.q.val),]
  gsea_neg_names <- as.character(gsea_neg$NAME)
  gsea_neg_names2 <- paste0(gsea_neg_names, '.xls')
  gsea_neg_names2 <- intersect(gsea_neg_names2, file_list)
  gsea_neg_names3 <- strtrim(gsea_neg_names2, nchar(gsea_neg_names2)-4)
  
  gsea_pos <- file_list[grepl('gsea_report_for_na_pos_\\d*.xls', file_list)]
  gsea_pos <- read.csv(paste0(directory,'/',gsea_pos), sep = '\t')
  gsea_pos <- gsea_pos[order(gsea_pos$FDR.q.val),]
  gsea_pos_names <- as.character(gsea_pos$NAME)
  gsea_pos_names2 <- paste0(gsea_pos_names, '.xls')
  gsea_pos_names2 <- intersect(gsea_pos_names2, file_list)
  gsea_pos_names3 <- strtrim(gsea_pos_names2, nchar(gsea_pos_names2)-4)
  
  neg_genes <- NULL
  pos_genes <- NULL
  if (length(gsea_neg_names2) > 0) {
    for (i in 1:length(gsea_neg_names2)) {
      genes <- read.csv(paste0(directory,'/',gsea_neg_names2[i]), sep = '\t')
      genes <- genes[rev(rownames(genes)),]
      genes <- data.frame(name = genes$PROBE, rank = (num_genes - genes$RANK.IN.GENE.LIST), metric = genes$RANK.METRIC.SCORE)
      colnames(genes)[1] <- gsea_neg_names3[i]
      if (is.null(neg_genes)) {
        neg_genes = genes
      }
      else {
        neg_genes <- cbind.fill(neg_genes, genes, fill = '')
      }
    }
  }
  
  if (length(gsea_pos_names2 > 0)) {
    for (i in 1:length(gsea_pos_names2)) {
      genes <- read.csv(paste0(directory,'/',gsea_pos_names2[i]), sep = '\t')
      genes <- data.frame(name = genes$PROBE, rank = genes$RANK.IN.GENE.LIST, metric = genes$RANK.METRIC.SCORE)
      colnames(genes)[1] <- gsea_pos_names3[i]
      if (is.null(pos_genes)) {
        pos_genes = genes
      }
      else {
        pos_genes <- cbind.fill(pos_genes, genes, fill = '')
      }
    }
  }
  write.table(neg_genes, file = paste0(outfile_prefix, '_down_genes.txt'), row.names = F, quote = F, sep = '\t')
  write.table(pos_genes, file = paste0(outfile_prefix, '_up_genes.txt'), row.names = F, quote = F, sep = '\t')
}

compare_genes <- function(file1, file2, name1, name2) {
  names1 <- read.delim(file1, stringsAsFactors = F)
  names2 <- read.delim(file2, stringsAsFactors = F)
  sets1 <- colnames(names1)[seq(1, length(colnames(names1)), 3)]
  sets2 <- colnames(names2)[seq(1, length(colnames(names2)), 3)]
  
  combined1 <- NULL
  combined2 <- NULL
  i = 1
  while(i < ncol(names1)) {
    temp <- names1[,i:(i+2)]
    temp <- temp[!temp[,1] == '',]
    temp <- cbind(temp, rep(colnames(temp)[1], nrow(temp)), stringsAsFactors = F)
    colnames(temp) <- c('GENE', 'RANK', 'METRIC', 'GENESET')
    combined1 <- rbind(combined1, temp, stringsAsFactors = F)
    i = i + 3
  }
  
  i = 1
  while(i < ncol(names2)) {
    temp <- names2[,i:(i+2)]
    temp <- temp[!temp[,1] == '',]
    temp <- cbind(temp, rep(colnames(temp)[1], nrow(temp)), stringsAsFactors = F)
    colnames(temp) <- c('GENE', 'RANK', 'METRIC', 'GENESET')
    combined2 <- rbind(combined2, temp, stringsAsFactors = F)
    i = i + 3
  }
  
  dup1 <- combined1[duplicated(combined1$GENE),]
  dup_genes1 <- unique(dup1$GENE)
  collapsed1 = NULL
  for (i in dup_genes1) {
    common <- combined1[combined1$GENE == i,]
    set <- paste(common$GENESET, collapse = ', ')
    common <- data.frame(common[1,1:3], GENESET = set)
    collapsed1 <- rbind(collapsed1, common)
  }
  
  combined1 <- rbind(collapsed1, combined1)
  combined1 <- combined1[!duplicated(combined1$GENE),]
  combined1 <- combined1[order(combined1$RANK),]
  
  dup2 <- combined2[duplicated(combined2$GENE),]
  dup_genes2 <- unique(dup1$GENE)
  collapsed2 = NULL
  for (i in dup_genes2) {
    common <- combined2[combined2$GENE == i,]
    set <- paste(common$GENESET, collapse = ', ')
    common <- data.frame(common[1,1:3], GENESET = set)
    collapsed2 <- rbind(collapsed2, common)
  }
  
  combined2 <- rbind(collapsed2, combined2)
  combined2 <- combined2[!duplicated(combined2$GENE),]
  combined2 <- combined2[order(combined2$RANK),]
  
  file_name1 <- gsub('_genes.txt','',file1)
  file_name2 <- gsub('_genes.txt','',file2)
  combined1 <- na.omit(combined1)
  combined2 <- na.omit(combined2)
  write.table(combined1, file = paste0(name1, '_BY_GENE.txt'), quote = F, sep = '\t', row.names = F)
  write.table(combined2, file = paste0(name2, '_BY_GENE.txt'), quote = F, sep = '\t', row.names = F)
  
  common_genes <- intersect(combined1$GENE, combined2$GENE)
  rownames(combined1) <- combined1$GENE
  rownames(combined2) <- combined2$GENE
  colnames(combined1) <- paste0(file_name1, '_', colnames(combined1))
  colnames(combined2) <- paste0(file_name2, '_', colnames(combined2))
  common_genes2 <- data.frame(common_genes, combined1[common_genes, 2:4], combined2[common_genes, 2:4])
  common_genes2 <- common_genes2[order((as.numeric(common_genes2[,2]) + as.numeric(common_genes2[,5]))/2),]
  
  write.table(common_genes2, file = paste(name1, 'vs', name2, 'OVERLAPPING_GENES.txt', sep = '_'), quote = F, sep = '\t', row.names = F)
}

compare_genes_all_cross_comparisons <- function(directory1, directory2) {
  file_list1 <- list.files(directory1)
  up_file_list1 <- file_list1[grepl('_up_', file_list1)]
  down_file_list1 <- file_list1[grepl('_down_', file_list1)]
  
  file_list2 <- list.files(directory2)
  up_file_list2 <- file_list2[grepl('_up_', file_list2)]
  down_file_list2 <- file_list2[grepl('_down_', file_list2)]
  
  current_dir <- getwd()
  dir.create(paste0(directory1,'_vs_',directory2, '_GENE_COMPARISONS'))
  setwd(paste0(directory1,'_vs_',directory2, '_GENE_COMPARISONS'))
  combinations_up <- expand.grid(up_file_list1, up_file_list2) 
  combinations_down <- expand.grid(down_file_list1, down_file_list2)
  for (i in 1:nrow(combinations_up)) {
    compare_genes(paste0('../',directory1,'/',combinations_up[i,1]), paste0('../', directory2,'/',combinations_up[i,2]), gsub('_genes.txt', '', as.character(combinations_up[i,1])), gsub('_genes.txt', '', as.character(combinations_up[i,2])))
  }
  
  for (i in 1:nrow(combinations_down)) {
    compare_genes(paste0('../',directory1,'/',combinations_down[i,1]), paste0('../', directory2,'/',combinations_down[i,2]), gsub('_genes.txt', '', as.character(combinations_down[i,1])), gsub('_genes.txt', '', as.character(combinations_down[i,2])))
  }
  setwd(current_dir)
}


gsea_lookup_keyword <- function(word) {
  dir <- '/Volumes/data0/users/dyao/GSEA/'
  genelist <- list(cellbkpt = rep(18900, 23), cellicna = rep(18900, 23), stem = c(15518, 24693, 18802, 19104, 18653), tumorbkpt = rep(18990, 35), tumoricna = rep(18990, 35))
  files <- list.files(dir)
  setnames <- NULL
  fullgenes <- list()
  for (i in 1:length(files)) {
    file_list <- list.files(paste0(dir, files[i]))
    for (j in 1:length(file_list)) {
      subfile_list <- list.files(paste0(dir, files[i], '/', file_list[j]))
      gsea_neg <- subfile_list[grepl('gsea_report_for_na_neg_\\d*.xls', subfile_list)]
      gsea_neg <- read.csv(paste0(dir, files[i], '/', file_list[j], '/', gsea_neg), sep = '\t')
      gsea_neg <- gsea_neg[gsea_neg$FDR.q.val < 0.25,]
      gsea_neg <- gsea_neg[grepl(word, gsea_neg$NAME),]
      gsea_neg_names <- as.character(gsea_neg$NAME)
      gsea_neg_names2 <- paste0(gsea_neg_names, '.xls')
      gsea_neg_names2 <- intersect(gsea_neg_names2, subfile_list)
      gsea_neg_names3 <- substr(gsea_neg_names2, 0, nchar(gsea_neg_names2)-4)
      
      gsea_pos <- subfile_list[grepl('gsea_report_for_na_pos_\\d*.xls', subfile_list)]
      gsea_pos <- read.csv(paste0(dir, files[i],'/', file_list[j], '/', gsea_pos), sep = '\t')
      gsea_pos <- gsea_pos[gsea_pos$FDR.q.val < 0.25,]
      gsea_pos <- gsea_pos[grepl(word, gsea_pos$NAME),]
      gsea_pos_names <- as.character(gsea_pos$NAME)
      gsea_pos_names2 <- paste0(gsea_pos_names, '.xls')
      gsea_pos_names2 <- intersect(gsea_pos_names2, subfile_list)
      gsea_pos_names3 <- substr(gsea_pos_names2, 0, nchar(gsea_pos_names2)-4)
      
      
      if (length(gsea_neg_names2) > 0) {
        temp <- data.frame(NAME = gsea_neg$NAME, NES = gsea_neg$NES, FDRqval = gsea_neg$FDR.q.val, TYPE = file_list[j])
        setnames <- rbind(setnames, temp)
        for (k in 1:length(gsea_neg_names2)) {
          genes <- read.csv(paste0(dir, files[i],'/', file_list[j], '/', gsea_neg_names2[k]), sep = '\t')
          genes <- genes[rev(rownames(genes)),]
          if (!(gsea_neg_names3[k] %in% names(fullgenes))) {
            temp <- data.frame((genelist[[i]][j] - genes$RANK.IN.GENE.LIST), genes$RANK.METRIC.SCORE, row.names = genes$PROBE)
            colnames(temp) <- c(paste0(file_list[j], '_RANK'), 'METRIC')
            fullgenes[[gsea_neg_names3[k]]] <- temp
          }
          
          else if (gsea_neg_names3[k] %in% names(fullgenes)) {
            temp <- data.frame((genelist[[i]][j] - genes$RANK.IN.GENE.LIST), genes$RANK.METRIC.SCORE, row.names = genes$PROBE)
            colnames(temp) <- c(paste0(file_list[j], '_RANK'), 'METRIC')
            geneorder <- rownames(fullgenes[[gsea_neg_names3[k]]])
            temp <- temp[geneorder,]
            fullgenes[[gsea_neg_names3[k]]] <- cbind(fullgenes[[gsea_neg_names3[k]]], temp)
          }
        }
      }
      
      if (length(gsea_pos_names2) > 0) {
        temp <- data.frame(NAME = gsea_pos$NAME, NES = gsea_pos$NES, FDRqval = gsea_pos$FDR.q.val, TYPE = file_list[j])
        setnames <- rbind(setnames, temp)
        for (k in 1:length(gsea_pos_names2)) {
          genes <- read.csv(paste0(dir, files[i],'/', file_list[j], '/', gsea_pos_names2[k]), sep = '\t')
          if (!(gsea_pos_names3[k] %in% names(fullgenes))) {
            temp <- data.frame(genes$RANK.IN.GENE.LIST, genes$RANK.METRIC.SCORE, row.names = genes$PROBE)
            colnames(temp) <- c(paste0(file_list[j], '_RANK'), 'METRIC')
            fullgenes[[gsea_pos_names3[k]]] <- temp
          }
          
          else if (gsea_pos_names3[k] %in% names(fullgenes)) {
            temp <- data.frame(genes$RANK.IN.GENE.LIST, genes$RANK.METRIC.SCORE, row.names = genes$PROBE)
            colnames(temp) <- c(paste0(file_list[j], '_RANK'), 'METRIC')
            geneorder <- rownames(fullgenes[[gsea_pos_names3[k]]])
            temp <- temp[geneorder,]
            fullgenes[[gsea_pos_names3[k]]] <- cbind(fullgenes[[gsea_pos_names3[k]]], temp)
          }
        }
      }
    }
  }
  write.table(setnames, file = paste0('GSEA_genesets_', word, '.txt'), row.names = F, sep = '\t', quote = F)
  for (x in 1:length(fullgenes)) {
    temp <- fullgenes[[x]]
    if (ncol(temp) > 2) {
      avgranks <- rowMeans(temp[,seq(1,ncol(temp),2)])
      avgranks <- avgranks[order(avgranks)]
      fullgenes[[x]] <- temp[order(rowMeans(temp[,seq(1,ncol(temp),2)])),]
      fullgenes[[x]] <- data.frame(fullgenes[[x]], AVGRANK = avgranks)
    }
    write.table(fullgenes[[x]], file = paste0(names(fullgenes)[x], '_genes.txt'), sep = '\t', quote = F, col.names = NA)
  }
}
