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


### Returns all MSigDB gene set names containing a given keyword
# |keyword| is the keyword contained in the gene set name (ex. 'SECRETION')
# |type| is 'c2', 'c5', 'cpg', or 'all'. 'c2', 'c5', and 'cpg' will return gene set names from their respective groups from MSigDB. 'all' will return gene set names from all three groups. 
# OUTPUT: A vector of gene set names
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


### Returns all genes in given gene sets from MSigDB
# |genesets| is a vector of gene set names
# OUTPUT: A single vector containing all the genes in all the gene sets 
list_genes_geneset <- function(genesets) {
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


### Outputs a summary of comparisons between two GSEA directories.
# |directory1| and |directory2| are the two GSEA directories
# |label1| is the label to use for GSEA directory 1 in table titles and file name
# |label2| is the label to use for GSEA directory 2 in table titles and file name
# |type| is either 'simple' or 'all'.
# OUTPUT: If |type| is set to 'simple', the output file will be a png containing a table with four columns.
# The leftmost column will have the top gene sets upregulated in both GSEA directories. Gene sets are sorted by average FDR qval. 
# The bottommost row has the number of 'overlapping' gene sets, or gene sets that are upregulated in both GSEA directories and both have a qval < 0.25.
# The remaining three columns are top downregulated gene sets, top gene sets upregulated in 1 and downregulated in 2, and vice versa for the last column.
# If |type| is set to 'all', the output is a directory containing 8 files. Four of the files are essentially the columns of the table from the 'simple' output,
# but with NESs, p-values, and FDR q-values listed as well. The other four files are simply all the upregulated and downregulated gene sets from 
# both GSEA directories with all associated statistics.
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


### Outputs a summary of comparisons between two GSEA directories.
# |directory1| and |directory2| are the two GSEA directories
# |label1| is the label to use for GSEA directory 1 in table titles and file name
# |label2| is the label to use for GSEA directory 2 in table titles and file name
# |type| is 'all', 'cpg', or 'c2'. If type is 'cpg' or 'c2', then summaries will only contain gene sets from those groups.
# OUTPUT: A txt file containing a table. The left half of the table will contain the top 100 upregulated gene sets between the two directories.
# Gene sets are sorted by average NES between both directories. The NES, FDR, and rank of each geneset for each directory are listed as well.
# The right half of the table is the same as the left for for the top 100 downregulated gene sets.
compare_gsea_v2 <- function(directory1, directory2, label1, label2, type = 'all') {
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
  temp_neg_file1 <- temp_neg_file1[!grepl('^PID',temp_neg_file1$NAME),]
  
  temp_pos_file1 <- gsea_pos_file1[!grepl('REACTOME',gsea_pos_file1$NAME),]
  temp_pos_file1 <- temp_pos_file1[!grepl('KEGG',temp_pos_file1$NAME),]
  temp_pos_file1 <- temp_pos_file1[!grepl('BIOCARTA',temp_pos_file1$NAME),]
  temp_pos_file1 <- temp_pos_file1[!grepl('^PID',temp_pos_file1$NAME),]
  
  temp_neg_file2 <- gsea_neg_file2[!grepl('REACTOME',gsea_neg_file2$NAME),]
  temp_neg_file2 <- temp_neg_file2[!grepl('KEGG',temp_neg_file2$NAME),]
  temp_neg_file2 <- temp_neg_file2[!grepl('BIOCARTA',temp_neg_file2$NAME),]
  temp_neg_file2 <- temp_neg_file2[!grepl('^PID',temp_neg_file2$NAME),]
  
  temp_pos_file2 <- gsea_pos_file2[!grepl('REACTOME',gsea_pos_file2$NAME),]
  temp_pos_file2 <- temp_pos_file2[!grepl('KEGG',temp_pos_file2$NAME),]
  temp_pos_file2 <- temp_pos_file2[!grepl('BIOCARTA',temp_pos_file2$NAME),]
  temp_pos_file2 <- temp_pos_file2[!grepl('^PID',temp_pos_file2$NAME),]
  
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


### Calls the function compare_gsea for all possible cross comparisons of GSEAs within two directories of GSEAs.
# |directory1| and |directory2| are both directories containing multiple subdirectories, each of which is the output from running GSEA on a single .rnk metric file.
# |type| is either 'simple' or 'all'. See compare_gsea.
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


### Returns the top gene sets for each individual GSEA directory in a directory
# |directory| is a directory containing multiple subdirectories, each of which is the output from running GSEA on a single .rnk metric file.
# |outfileprefix| is the prefix for the output file name
# OUTPUT: Two txt files containing tables. Each table has each of the GSEA directories as columns and the top 20 upregulated/downregulated gene sets by NES for each GSEA listed underneath.
# NOTE: no cross-GSEA comparisons happen, this is simply to list the names
gsea_obtain_top <- function(directory, outfileprefix) {
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
  
  write.table(combined_neg, file = paste0(outfileprefix, '_up_all.txt'), quote = F, row.names = F, sep = '\t')
  write.table(combined_pos, file = paste0(outfileprefix, '_down_all.txt'), quote = F, row.names = F, sep = '\t')
}

### Outputs a summary of statistics from a GSEA directory
# |directory1| is a GSEA directory
# |label| is a label for the directory to use in table titles and file name
# OUTPUT: A png file containing a table. The left half contains the top 20 upregulated gene sets ordered by FDR. The FDR is also listed.
# The right half contains the top 20 downregulated gene sets.
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