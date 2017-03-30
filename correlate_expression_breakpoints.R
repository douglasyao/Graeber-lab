library('tools')
library('waterfall')
library('infotheo')

### Computes the Pearson correlation of the expression of each gene with breakpoints
# |prefix| is name of prefix for output file
# |log2_flag| is whether to log transform expresion matrix
# |temp| is the gene expression matrix with columns as samples and rows as genes
# |breakpoints| is the output file from breakpoints_by_chrom.py which contains breakpoint information
compute_icna_bkpt <- function(prefix, log2_flag, temp, breakpoints) {

  if(log2_flag==T) {
    temp<-round(log2(temp+1),5)
    prefix<- paste0(prefix,"_log2")
  } 
  
  countTable<-data.matrix(temp)
  countTable<-countTable[!rowSums(!is.finite(countTable)),]
  
  breakpoints_transpose <-t(breakpoints)
  breakpoints_transpose_names<- colnames(breakpoints_transpose)
  cnt_tble_t_names<- colnames(countTable)
  common_names_corr <- intersect(breakpoints_transpose_names,cnt_tble_t_names)
  breakpoints.idx <- match(common_names_corr,breakpoints_transpose_names)
  countTable.idx <-match(common_names_corr,cnt_tble_t_names)
  breakpoints_common <-breakpoints_transpose[,breakpoints.idx]
  countTable_common <-countTable[,countTable.idx]
  
  removeNAs<-function(dataset){
    genes_with_NA <- apply(dataset,1, function(x) any(is.na(x))) #remove features with NAs
    dataset <- dataset[!genes_with_NA,]
    return(dataset)
  }

  correlation_table <-data.matrix(matrix(0,nrow(countTable_common),6))
  rownames(correlation_table) <- rownames(countTable_common)
  colnames(correlation_table) <- c("bkpt_vs_RNA_rsquared","bkpt_vs_RNA_p_val","bkpt_vs_RNA_rvalue","ICNA_vs_RNA_rsquared","ICNA_vs_RNA_p_val","ICNA_vs_RNA_rvalue")
  for (i in 1:nrow(countTable_common)) {
    if(log2_flag==T) {    
      correlation_table[i,1] <- summary(lm(log2(breakpoints_common[c("bkpt_samp"),]+1)  ~ countTable_common[i,]))$adj.r.squared
      correlation_table[i,2] <- cor.test(log2(breakpoints_common[c("bkpt_samp"),]+1),countTable_common[i,])$p.value
      correlation_table[i,3] <- cor.test(log2(breakpoints_common[c("bkpt_samp"),]+1), countTable_common[i,])$estimate[[1]]
      correlation_table[i,4] <- summary(lm(log2(breakpoints_common[c("ICNA_score"),]+1)~ countTable_common[i,]))$adj.r.squared
      correlation_table[i,5] <- cor.test(log2(breakpoints_common[c("ICNA_score"),]+1),countTable_common[i,])$p.value
      correlation_table[i,6] <- cor.test(log2(breakpoints_common[c("ICNA_score"),]+1), countTable_common[i,])$estimate[[1]]
      }else{
        correlation_table[i,1] <- summary(lm(breakpoints_common[c("bkpt_samp"),]  ~ countTable_common[i,]))$adj.r.squared
        correlation_table[i,2] <- cor.test(breakpoints_common[c("bkpt_samp"),],countTable_common[i,])$p.value
        correlation_table[i,3] <- cor.test(breakpoints_common[c("bkpt_samp"),],countTable_common[i,])$estimate[[1]]
        correlation_table[i,4] <- summary(lm(breakpoints_common[c("ICNA_score"),]~ countTable_common[i,]))$adj.r.squared
        correlation_table[i,5] <- cor.test(breakpoints_common[c("ICNA_score"),],countTable_common[i,])$p.value
        correlation_table[i,6] <- cor.test(breakpoints_common[c("ICNA_score"),],countTable_common[i,])$estimate[[1]]   
      }   
  }
  correlation_table<-data.frame(correlation_table)
  correlation_table_bkpt<-correlation_table[order(correlation_table$bkpt_vs_RNA_rvalue,decreasing=TRUE), ]
  write.table(removeNAs(correlation_table_bkpt),paste0(prefix,"_RNA_pearson_correlation.txt"), col.names=NA,row.names=T,sep="\t", quote = F)
}
