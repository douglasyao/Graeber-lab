### small_cell.R
# Perform data processing for specific datasets related to small cell cancer



### merging stuff
tcga.icna <- read.delim('~/ICNA_vs_expression_correlation_surface_genes_from_bryan_only.txt', stringsAsFactors = F, row.names = 1)
tcga.icna <- tcga.icna[,'PRAD',drop = F]
tcga.bkpt <- read.delim('~/bkpt_vs_expression_correlation_surface_genes_from_bryan_only.txt', stringsAsFactors = F, row.names = 1)
tcga.bkpt <- tcga.bkpt[,'PRAD',drop = F]

beltran.icna <- read.delim('~/beltran2016_instability/CRPC_NEPC_ICNA_surface_genes_from_bryan.rnk', stringsAsFactors = F, row.names = 1)
beltran.bkpt <- read.delim('~/beltran2016_instability/CRPC_NEPC_bkpt_surface_genes_from_bryan.rnk', stringsAsFactors = F, row.names = 1)
allgenes <- union(rownames(tcga.bkpt), rownames(beltran.bkpt))

all <- data.frame(row.names = allgenes)
all <- cbind(all, tcga.bkpt[allgenes,])
all <- cbind(all, beltran.bkpt[allgenes,])
all <- cbind(all, rowMeans(all))
all <- cbind(rownames(all),all)
colnames(all) <- c('Gene', 'TCGA_PRAD_bkpt','beltran_NEPC_Adeno_bkpt','Ave')
all <- all[order(-all$Ave),]
write.table(all, 'tcga_beltran_merged_surface_genes_bkpt.txt', quote = F, sep = '\t', row.names = F)

###stuff
setwd('~/beltran2016_instability/')
files <- list.files(pattern = 'rnk')

surface <- read.csv('~/Downloads/Gene Lists for Cell Surface Markers_Kinases_Nuclear Hormone Receptors.csv', stringsAsFactors = F, header = F)
surface <- surface$V1
uni <- read.delim('~/Downloads/uniprot-organism%3A-homo+sapiens-+AND+locations%3A%28location%3A-cell+memb--.tab', stringsAsFactors = F)
uni <- gsub('_HUMAN','',uni$Entry.name)

for (file in files) {
  f <- read.delim(file, stringsAsFactors = F)
  f.surface <- f[f[,1] %in% surface,]
  name <- gsub('.rnk', '_surface_genes_from_bryan.rnk', file)
  write.table(f.surface, name, quote = F, sep = '\t', row.names = F)
  f.uni <- f[f[,1] %in% uni,]
  name <- gsub('.rnk', '_surface_genes_from_uniprot.rnk', file)
  write.table(f.uni, name, quote = F, sep = '\t', row.names = F)
}




bkpt <- bkpt[rownames(bkpt) %in% union,]
icna <- icna[rownames(icna) %in% union,]

write.table(bkpt, file = 'bkpt_vs_expression_correlation_surface_genes_union_only.txt', quote = F, col.names = NA, sep = '\t')
write.table(icna, file = 'ICNA_vs_expression_correlation_surface_genes_union_only.txt', quote = F, col.names = NA, sep = '\t')

###for grant
library(pheatmap)
beltran <- read.delim('~/Dropbox/Doug/1.10.17/Beltran_2016_rsem_genes_upper_norm_counts.txt', stringsAsFactors = F, row.names = 1)
beltran <- na.omit(data.matrix(beltran))
beltran <- beltran[rowSums(beltran) > 0,]
beltran <- clean_dat(beltran)

john <- read.delim('~/Dropbox/Doug/1.10.17/john_lee_rsem_genes_upper_norm_counts.txt', stringsAsFactors = F, row.names = 1)
john <- na.omit(data.matrix(john))
john <- john[rowSums(john) > 0,]
john <- clean_dat(john)

jung <- read.delim('~/Dropbox/Doug/1.10.17/jung_wook_rsem_genes_upper_norm_counts.txt', stringsAsFactors = F, row.names = 1)
jung <- na.omit(data.matrix(jung))
jung <- jung[rowSums(jung) > 0,]
jung <- clean_dat(jung)

sig <- read.delim('/Volumes/data0/users/dyao/lung_prostate_small_cell_rank_files/ttest_Beltran2016_NEPC_vs_CRPC_genes.rnk', stringsAsFactors = F)
sig <- read.delim('/Volumes/data0/users/dyao/lung_prostate_small_cell_rank_files/ttest_beltran2011_sc_vs_organ.rnk', stringsAsFactors = F)
sig <- get_ave_sig_by_rank('~/Documents/small_cell/')
genes <- intersect_all(sig[1:50,1], rownames(beltran), rownames(john), rownames(jung))

a <- beltran[genes,]
b <- john[genes,]
c <- jung[genes,]

all <- cbind(a[,grepl('_C', colnames(a))], b[,1:4], a[,grepl('_N', colnames(a))], b[,5:9], c)
allv2 <- cbind(a[,grepl('_C', colnames(a))], b[,1:4], a[,grepl('_N', colnames(a))], c)
all <- t(scale(t(all)))
allv2 <- t(scale(t(allv2)))
pheatmap(all, cluster_rows = F, cluster_cols = F, gaps_col = c(34,38,53,58), main = 'Log upper quartile normalized counts of top 50 genes in Beltran 2016')
pheatmap(allv2, cluster_rows = F, cluster_cols = F, gaps_col = c(34,38,53), main = 'Log upper quartile normalized counts of top 50 genes in Beltran 2016')

###genomic instability
setwd('~/CCLE/')
files <- list.files(pattern = 'correlation')
names <- gsub('_RNA_pearson_correlation.txt','',files)
names.icna <- paste0(names,'_ICNA.rnk')
names.bkpt <- paste0(names,'_bkpt.rnk')
for (i in 1:length(files)) {
  temp <- read.delim(files[i])
  temp.icna <- data.frame(gene = temp$X, rvalue = temp$ICNA_vs_RNA_rvalue)
  temp.icna <- temp.icna[order(-temp.icna$rvalue),]
  write.table(temp.icna, names.icna[i], quote = F, sep = '\t', row.names = F)
  temp.bkpt <- data.frame(gene = temp$X, rvalue = temp$bkpt_vs_RNA_rvalue)
  temp.bkpt <- temp.bkpt[order(-temp.bkpt$rvalue),]
  write.table(temp.bkpt, names.bkpt[i], quote = F, sep = '\t', row.names = F)
}

##NBL
annot <- read.delim('~/Downloads/beltran2016goodsamples.txt', header = F, stringsAsFactors = F, sep = '\n')
prefix<-"NEPC"
breakpoints <- read.delim('~/Downloads/NBL_bkpts_sorted_by_bkpts_per_chrom.txt', row.names=1,check.names=F,stringsAsFactors = F)
rownames(breakpoints) <- strtrim(rownames(breakpoints),16)
annot <- read.csv('~/Downloads/TARGET_NBL_ClinicalData_20151124.csv', stringsAsFactors = F)
annot <- data.frame(name = annot$TARGET.USI, type = annot$MYCN.status)
amp <- as.character(annot[annot$type == 'Amplified',1])
notamp <- as.character(annot[annot$type == 'Not Amplified',1])
temp <- read.table("~/Dropbox/Doug/NBL_normalized_results_processed.txt", header=TRUE, row.names=NULL, sep='\t', check.names = FALSE,stringsAsFactors = FALSE)
source('~/Dropbox/Doug/Sanaz/collapse_data.R')
temp <- collapse_data(temp,group = 'sample')
colnames(temp) <- strtrim(colnames(temp),16)
temp <- temp[,colnames(temp) %in% notamp | colnames(temp) %in% amp]

##CCLE
temp <- read.delim('~/Documents/CCLE_Expression_Entrez_2012-09-29.gct', header=TRUE, row.names=NULL, check.names = FALSE,stringsAsFactors = FALSE)
temp <- temp[,-1]
temp <- data.frame(collapse_data(temp, type = 'maxavg', group = 'Description'))
copy <- temp
temp <- copy
temp <- temp[,colnames(temp) %in% sc]

annot <- read.delim('~/Documents/CCLE_sample_info_file_2012-10-18.txt', stringsAsFactors = F)
annot <- annot[annot$Site.Primary == 'lung',]
adeno <- annot[annot$Hist.Subtype1 == 'adenocarcinoma',1]
squam <- annot[annot$Hist.Subtype1 == 'squamous_cell_carcinoma',1]
sc <- annot[annot$Hist.Subtype1 == 'small_cell_carcinoma',1]

breakpoints <- read.delim(file = '~/Dropbox/Doug/ICNA_BKPT/cell_line_bkpts_sorted_by_bkpts_per_chrom.txt',header=T,row.names=1,check.names=F,stringsAsFactors = F)



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
  correlation_table[i,1] <- summary(lm(breakpoints_common[c("bkpt_samp"),]  ~ countTable_common[i,]))$adj.r.squared
  correlation_table[i,2] <- cor.test(breakpoints_common[c("bkpt_samp"),],countTable_common[i,])$p.value
  correlation_table[i,3] <- cor.test(breakpoints_common[c("bkpt_samp"),],countTable_common[i,])$estimate[[1]]
  correlation_table[i,4] <- summary(lm(breakpoints_common[c("ICNA_score"),]~ countTable_common[i,]))$adj.r.squared
  correlation_table[i,5] <- cor.test(breakpoints_common[c("ICNA_score"),],countTable_common[i,])$p.value
  correlation_table[i,6] <- cor.test(breakpoints_common[c("ICNA_score"),],countTable_common[i,])$estimate[[1]]   
}
correlation_table<-data.frame(correlation_table)
correlation_table_bkpt<-correlation_table[order(correlation_table$bkpt_vs_RNA_rvalue,decreasing=TRUE), ]

write.table(removeNAs(correlation_table_bkpt),"CCLE_SCLC_RNA_pearson_correlation.txt", col.names=NA,row.names=T,sep="\t", quote = F)

  
  
  
  
  
} else if (method=="spearman"){   
  correlation_table <-data.matrix(matrix(0,nrow(countTable_common),4))
  rownames(correlation_table) <- rownames(countTable_common)
  colnames(correlation_table) <- c("bkpt_vs_RNA_p_val","bkpt_vs_RNA_rho_val","ICNA_vs_RNA_pval","ICNA_vs_RNA_rho_val")  
  for (i in 1:nrow(countTable_common)) {
    if(log2_flag==T) {    
      correlation_table[i,1] <- cor.test(log2(breakpoints_common[c("bkpt_samp"),]+1),countTable_common[i,],method="spearman",exact=F)$p.value
      correlation_table[i,2] <- cor.test(log2(breakpoints_common[c("bkpt_samp"),]+1), countTable_common[i,],method="spearman",exact=F)$estimate[[1]]
      correlation_table[i,3] <- cor.test(log2(breakpoints_common[c("ICNA_score"),]+1),countTable_common[i,],method="spearman",exact=F)$p.value
      correlation_table[i,4] <- cor.test(log2(breakpoints_common[c("ICNA_score"),]+1), countTable_common[i,],method="spearman",exact=F)$estimate[[1]]
    }else{    
      correlation_table[i,1] <- cor.test(breakpoints_common[c("bkpt_samp"),],countTable_common[i,],method="spearman",exact=F)$p.value
      correlation_table[i,2] <- cor.test(breakpoints_common[c("bkpt_samp"),],countTable_common[i,],method="spearman",exact=F)$estimate[[1]]     
      correlation_table[i,3] <- cor.test(breakpoints_common[c("ICNA_score"),],countTable_common[i,],method="spearman",exact=F)$p.value
      correlation_table[i,4] <- cor.test(breakpoints_common[c("ICNA_score"),],countTable_common[i,],method="spearman",exact=F)$estimate[[1]]    
    }   
  }
  correlation_table<-data.frame(correlation_table)
  correlation_table_bkpt<-correlation_table[order(correlation_table$bkpt_vs_RNA_rho_val,decreasing=TRUE), ]
  filter_set_names=rownames(filter_set)
  corr_names=rownames(correlation_table)
  filter_corr_common_names=intersect(corr_names,filter_set_names)
  filter_corr.idx <- match(filter_corr_common_names,corr_names)
  filter_correlation_table <- correlation_table[filter_corr.idx,]
  filter_correlation_table_bkpt<-filter_correlation_table[order(filter_correlation_table$bkpt_vs_RNA_rho_val,decreasing=TRUE), ] 
  write.table(removeNAs(correlation_table_bkpt),paste0(prefix,"_RNA_spearman_correlation.txt"), col.names=NA,row.names=T,sep="\t")
  write.table(removeNAs(filter_correlation_table_bkpt),paste0(prefix,"_RNA_genefilter_spearman_correlation.txt"),col.names=NA,row.names=T,sep='\t')
}else if (method=="mutual"){
  correlation_table <-data.matrix(matrix(0,nrow(countTable_common),2))
  rownames(correlation_table) <- rownames(countTable_common)
  colnames(correlation_table) <- c("bkpt_mutinf","ICNA_mutinf") 
  for (i in 1:nrow(countTable_common)) {
    if(log2_flag==T) {    
      discretized_bkpt<-infotheo::discretize(cbind(log2(breakpoints_common[c("bkpt_samp"),]+1), countTable_common[i,]))
      discretized_ICNA<-infotheo::discretize(cbind(log2(breakpoints_common[c("ICNA_score"),]+1), countTable_common[i,]))
      correlation_table[i,1] <- mutinformation(discretized_bkpt[,1],discretized_bkpt[,2])   
      correlation_table[i,2] <- mutinformation(discretized_ICNA[,1],discretized_ICNA[,2])  
    }else{
      discretized_bkpt<-infotheo::discretize(cbind(breakpoints_common[c("bkpt_samp"),], countTable_common[i,]))
      discretized_ICNA<-infotheo::discretize(cbind(breakpoints_common[c("ICNA_score"),], countTable_common[i,]))
      correlation_table[i,1] <- mutinformation(discretized_bkpt[,1],discretized_bkpt[,2])   
      correlation_table[i,2] <- mutinformation(discretized_ICNA[,1],discretized_ICNA[,2])  
    } 
  }  
  correlation_table<-data.frame(correlation_table)
  correlation_table_bkpt<-correlation_table[order(correlation_table$bkpt_mutinf,decreasing=TRUE), ] 
  filter_set_names=rownames(filter_set)
  corr_names=rownames(correlation_table)
  filter_corr_common_names=intersect(corr_names,filter_set_names)
  filter_corr.idx <- match(filter_corr_common_names,corr_names)
  filter_correlation_table <- correlation_table[filter_corr.idx,]
  filter_correlation_table_bkpt<-filter_correlation_table[order(filter_correlation_table$bkpt_mutinf,decreasing=TRUE), ]
  write.table(removeNAs(correlation_table_bkpt),paste0(prefix,"_RNA_mut_info_correlation.txt"), col.names=NA,row.names=T,sep="\t") 
  write.table(removeNAs(filter_correlation_table_bkpt),paste0(prefix,"_RNA_genefilter_mutinfo_correlation.txt"),col.names=NA,row.names=T,sep='\t')  
}




# ############
# #writes out rank files of correlations
# #corr_bkpt.rnk <- data.frame(rownames(correlation_table_bkpt),correlation_table_bkpt$bkpt_vs_RNA_rquared,row.names=1)
# corr_bkpt.rnk <- data.frame(rownames(correlation_table_bkpt),correlation_table_bkpt$bkpt_vs_RNA_rsquared,row.names=1)
# corr_ICNA.rnk <- data.frame(rownames(correlation_table_ICNA),correlation_table_ICNA$ICNA_vs_RNA_rsquared,row.names=1)
# if(log2_flag==T) {
#   write.table(corr_bkpt.rnk,paste0(prefix,"_RNA_log2_bkpt_correlation.rnk"), col.names=F,row.names=T,sep="\t",quote=F)  
#   write.table(corr_ICNA.rnk,paste0(prefix,"_RNA_log2_ICNA_correlation.rnk"), col.names=F,row.names=T,sep="\t",quote=F)
# } else{
#   write.table(corr_bkpt.rnk,paste0(prefix,"_RNA_bkpt_correlation.rnk"), col.names=F,row.names=T,sep="\t",quote=F)  
#   write.table(corr_ICNA.rnk,paste0(prefix,"_RNA_ICNA_correlation.rnk"), col.names=F,row.names=T,sep="\t",quote=F)
# }




#remove_samples("bkpts_sorted_by_bkpts_per_chrom_noBRCA.txt","BRCA_TP53_MUTATED_FROM_MUTSIG.txt",substr_flag = T,substring_sample = 1,start = 1,end = 15,output_file_name = "bkpts_sorted_by_bkpts_per_chrom_noBRCA_noP53.txt",data_header = T)
#remove_samples("bkpts_sorted_by_ICNA_per_sample_noBRCA.txt","BRCA_TP53_MUTATED_FROM_MUTSIG.txt",substr_flag = T,substring_sample = 1,start = 1,end = 15,output_file_name = "bkpts_sorted_by_ICNA_per_sample_noBRCA_noP53.txt",data_header = T)
# 
# 
# breakpoints_common["bkpt_samp",]
# countTable_common["ANKZF1",]
# 
# plot(breakpoints_common["bkpt_samp",],countTable_common["ANKZF1",])






####VIPER

library(viper)
library(mixtools)

aracne.dir <- '/Volumes/data0/users/dyao/aracne_small_cell/'
data.dir <- '/Volumes/data0/users/dyao/small_cell_datasets/'
aracne.files <- list.files(aracne.dir)
data.files <- list.files(data.dir)

data <- as.matrix(read.delim(paste0(data.dir, data.files[1]), row.names = 1, stringsAsFactors = F))
colnames(data) <- gsub('Bel\\.','',colnames(data))
reg <- aracne2regulon(paste0(aracne.dir, aracne.files[1]), data)
annot <- read.delim('~/Dropbox/Doug/Prostate_Lung_Comparisons/Beltran.2011_Annotations.txt', stringsAsFactors = F)
org <- annot[annot$Diagnosis == 'PCa',]
sc <- annot[annot$Diagnosis == 'small cell PCa',]

sig <- rowTtest(data[,sc$SampleID], data[,org$SampleID])
sig <- (qnorm(sig$p.value/2, lower.tail = F) * sign(sig$statistic)[,1])
null <- ttestNull(data[,sc$SampleID], data[,org$SampleID])

ms <- msviper(sig,reg,null)
shadow <- shadow(ms, regulators = 25)

avesig <- get_ave_sig_by_rank('~/Documents/small_cell/')
avesig$gene <- as.character(avesig$gene)


### PCA bullshit
source('~/Dropbox/Doug/intersect_log_doPCA_from_file_and_project_second_dataset.R')
source('~/Dropbox/Doug/plot_pca.R')
library(Rtsne)
library(ggplot2)
library(dbscan)
library(DESeq2)
library(pheatmap)

up <- c('ASCL1','DLL3','FOXA2','SYP','NCAM1','SEZ6','SEZ6L','CHGA','CHGB','KCNH2','NKX2-1','PEG10','CEACAM5','DDC','NEUROD1','UCHL1','INSM1','BEX1','SOX2','SLITRK6')
down <- c('AR','TACSTD2','KLK3','PSCA','NKX3-1','REST','FOLH1','STEAP2')

cell.annot <- read.csv('~/Downloads/aad0501_Table_S3.csv', stringsAsFactors = F)
all.annot <- list()
for (i in 1:ncol(cell.annot)) {
  temp <- cell.annot[,i]
  temp <- gsub('\'','',temp)
  temp <- temp[!temp == '']
  all.annot[[colnames(cell.annot)[i]]] <- temp
}


#### calculating TPM from raw counts
library(edgeR)
library(GenomicRanges)
library(plyr)

#get gene lengths to calculate RPKM
a <- read.delim('~/Dropbox/mart_export1.txt')
b <- read.delim('~/Dropbox/mart_export2.txt')
genelengthData <- rbind(a,b)
genelengthData <- unique(genelengthData)
genelengthData <- na.omit(genelengthData)
seqNames <- genelengthData$Associated.Gene.Name
chrStart <- genelengthData$Exon.Chr.Start..bp.
chrEnd <- genelengthData$Exon.Chr.End..bp.
exonLengths <- GRanges(seqnames = Rle(seqNames), ranges = IRanges(start = chrStart, end = chrEnd))
exonLengths <- reduce(exonLengths)
individualseqnames <- as.factor((seqnames(exonLengths)))
individualexonLengths <- width(ranges(exonLengths))
totalgeneData <- data.frame(individualseqnames,individualexonLengths)
genelengthData <- ddply(totalgeneData, 'individualseqnames', numcolwise(sum))
colnames(genelengthData) <- c('gene','sum_exon_lengths')
write.table(genelengthData, file = 'Biomart_gene_lengths.txt', sep = '\t', quote = F, row.names = F)

lengths <- read.delim('~/Dropbox/Doug/Biomart_gene_lengths.txt', stringsAsFactors = F, row.names = 1)
a <- read.delim('~/Dropbox/Doug/DropSeq LungCancerBiopsy20160822.collapsed.human.dge.tsv', row.names = 1)
b <- read.delim('~/Downloads/Lung_Biopsy_v2.collapsed.human.dge.tsv', row.names = 1)
#ss <- read.delim('~/Downloads/Liver_9_4.collapsed.human.dge.tsv', row.names = 1)
ss <- read.delim('~/Downloads/Liver_9_22.collapsed.human.dge.tsv', row.names = 1)

common <- intersect(rownames(a), rownames(b))
a <- a[common,]
b <- b[common,]
ss <- cbind(a,b)


rownames(ss) <- fix_date_gene(rownames(ss))
ss <- quant_norm(ss)
ss <- clean_dat(ss)
tpm <- ss


ss <- ss[,order(-colSums(ss))]
#ss <- ss[,1:156]
common <- intersect(rownames(ss), rownames(lengths))
ss <- ss[common,]
lengths <- lengths[common,]

rpkm <- rpkm(ss,lengths)
tpm <- sweep(rpkm,2,colSums(rpkm),'/') * 1000000
tpm <- clean_dat(tpm)
#tpm <- t(scale(t(tpm),center = F))


callback = function(hc, data.plog) {
  colavg <- colSums(abs(data.plog)) / ncol(data.plog)
  rowavg <- rowSums(abs(data.plog)) / nrow(data.plog)
  dend = reorder(as.dendrogram(hc), wts = colavg)
  dend = reorder(as.dendrogram(dend), wts = rowavg)
  as.hclust(dend)
}

allgenes <- unlist(all.annot)
allgenes2 <- unlist(all.annot[1:5])
temp <- all.annot[1:3]
temp[[4]] <- rev(avesig$gene[(nrow(avesig)-49):(nrow(avesig))])
temp[[5]] <- avesig$gene[1:50]
common <- intersect(unlist(temp), rownames(tpm))
temp <- lapply(temp, function(x) return(x[x %in% common]))
indices <- cumsum(c(sapply(temp, function(x) return (length(x)))))
indices <- indices[-length(indices)]


sc.genes <- tpm[rownames(tpm) %in% allgenes,]
#sc.genes <- sc.genes[,colSums(sc.genes) > 0]
#pheatmap(sc.genes[1:10,1:10], clustering_callback = callback, cluster_rows = F, gaps_row = c(5,6))


sc.genes <- tpm[common,]
#sc.genes <- sc.genes[,colSums(sc.genes) > 0]
temp <- pheatmap(sc.genes, cluster_rows = F, show_rownames = F, show_colnames = F, clustering_distance_cols = 'correlation', gaps_row = indices, main = '9.22 Liver Biopsy Upper Quartile Normalized Counts')
macrop <- colnames(sc.genes)[temp$tree_col$order[1:32]]
tcell <- colnames(sc.genes)[temp$tree_col$order[33:92]]
rest <- colnames(sc.genes)[temp$tree_col$order[93:156]]



top50 <- c(rownames(pc1)[1:100],rownames(pc1)[(length(rownames(pc1))-100):length(rownames(pc1))])
sc.genes <- tpm[rownames(tpm) %in% top50,]
sc.genes <- sc.genes[,colSums(sc.genes) > 0]
pheatmap(sc.genes, fontsize = 5)


###single cell rna seq
ss <- read.delim('~/Dropbox/Doug/DropSeq LungCancerBiopsy20160822.collapsed.human.dge.tsv', row.names = 1)
rownames(ss) <- fix_date_gene(rownames(ss))
ss <- ss[,order(-colSums(ss))]
ss <- ss[,1:156]
sizes <- colSums(ss)
ss.raw <- ss
ss <- quant_norm(ss)
ss <- clean_dat(ss)
ss <- unique(ss)
pca.ss <- prcomp(t(ss), scale = T)
pc1 <- data.frame(gene = rownames(pca.ss$rotation), PC1 = pca.ss$rotation[,1])
pc1 <- pc1[order(pc1$PC1),]
pc2 <- data.frame(gene = rownames(pca.ss$rotation), PC2 = pca.ss$rotation[,2])
pc2 <- pc2[order(pc2$PC2),]

write.table(pc1, file = '/Volumes/data0/users/dyao/single_cell_rnaseq_pc1.rnk', quote = F, sep = '\t', row.names = F)
write.table(pc2, file = '/Volumes/data0/users/dyao/single_cell_rnaseq_pc2.rnk', quote = F, sep = '\t', row.names = F)

tsne.ss <- Rtsne(t(ss), initial_dims = 5, theta = 0, max_iter = 2500, perplexity = 30, verbose = T)
qplot(x = tsne.ss$Y[,1], y = tsne.ss$Y[,2], colour = as.factor(dbscan.ss$cluster), size = sizes)
dbscan.ss <- dbscan(tsne.ss$Y, eps = 2)

gene.cluster <- sapply(colnames(sc.genes), function(x) if (x %in% macrop) 'macrophage' else if (x %in% tcell) 'tcell' else 'rest')
qplot(x = PC1, y = PC2, data = data.frame(pca.ss$x[,1:2]), colour = as.factor(dbscan.ss$cluster), size = sizes)
qplot(x = PC2, y = PC3, data = data.frame(pca.ss$x[,2:3]), colour = as.factor(dbscan.ss$cluster), size = sizes)
qplot(x = PC2, y = PC4, data = data.frame(pca.ss$x[,c(2,4)]), colour = as.factor(dbscan.ss$cluster), size = sizes)
qplot(x = PC3, y = PC4, data = data.frame(pca.ss$x[,3:4]), colour = as.factor(dbscan.ss$cluster), size = sizes)


####sc gene based
cluster <- sapply(colnames(ss), function(x) if (x %in% scsc) 'sc' else 'nonsc')
qplot(x = PC1, y = PC2, data = data.frame(pca.ss$x[,1:2]), colour = as.factor(cluster))

column.data <- data.frame(group = dbscan.ss$cluster, row.names = colnames(ss.raw))
column.data <- column.data[column.data$group == 1 | column.data$group == 2 | column.data$group == 3,,drop = F]
column.data$group <- as.factor(column.data$group)
ss.raw <- ss.raw[,colnames(ss.raw) %in% rownames(column.data)]

dds <- DESeqDataSetFromMatrix(ss.raw, column.data, ~group)
dds <- DESeq(dds)
res1v2 <- results(dds, contrast = c('group','1','2'))
res2v3 <- results(dds, contrast = c('group','2','3'))
res1v3 <- results(dds, contrast = c('group','1','3'))
res1v2 <- na.omit(res1v2)
res2v3 <- na.omit(res2v3)
res1v3 <- na.omit(res1v3)
log_pvalues_from_deseq(res1v2, output = '/Volumes/data0/users/dyao/single_cell_rnaseq_group1v2.rnk')
log_pvalues_from_deseq(res1v3, output = '/Volumes/data0/users/dyao/single_cell_rnaseq_group1v3.rnk')
log_pvalues_from_deseq(res2v3, output = '/Volumes/data0/users/dyao/single_cell_rnaseq_group2v3.rnk')


##CCLE 
ccle <- read.delim('~/Downloads/CCLE_RNAseq_TPM_gene_count_matrix_lung.txt', stringsAsFactors = F, row.names = 2)
ccle <- ccle[,!grepl('PROSTATE', colnames(ccle))]
ccle <- ccle[2:ncol(ccle)]
ccle <- ccle[rowSums(ccle) > 0,] 
rownames(ccle) <- fix_date_gene(rownames(ccle))
ccle <- quant_norm(ccle)
ccle <- clean_dat(ccle)

annot <- read.delim('~/Dropbox/Doug/Prostate_Lung_Comparisons/CCLE_sample_info_file_2012-10-18.txt', stringsAsFactors = F, row.names = 1)
annot <- annot[grepl('LUNG',rownames(annot)),]
annot <- annot[rownames(annot) %in% colnames(ccle),]
annot <- annot[colnames(ccle),]
annot <- annot[grepl('^adenocarcinoma$|^small_cell_carcinoma$|^squamous_cell_carcinoma$', annot$Hist.Subtype1),]

ccle <- ccle[,colnames(ccle) %in% rownames(annot)]


ccle.pca <- prcomp(t(ccle))
qplot(x = PC1, y = PC2, data = data.frame(ccle.pca$x[,1:2]), colour = as.factor(annot$Hist.Subtype1))

common <- intersect(rownames(ccle), rownames(ss))
ccle.common <- ccle[common,]
ss.common <- ss[common,]
ccle.pca.common <- prcomp(t(ccle.common))
ss.pca.common <- prcomp(t(ss.common))
ss.project.ccle <- scale(t(ss.common), ccle.pca.common$center, ccle.pca.common$scale) %*% ccle.pca.common$rotation 
ccle.project.ss <- scale(t(ccle.common), ss.pca.common$center, ss.pca.common$scale) %*% ss.pca.common$rotation 

qplot(x = PC1, y = PC2, data = data.frame(ss.project.ccle[,1:2]), colour = as.factor(dbscan.ss$cluster), size = sizes)
qplot(x = PC1, y = PC2, data = data.frame(rbind(ss.project.ccle[,1:2], ccle.pca.common$x[,1:2])), colour = as.factor(c(dbscan.ss$cluster,annot$Hist.Subtype1)))
qplot(x = PC1, y = PC2, data = data.frame(ccle.project.ss[,1:2]), colour = as.factor(annot$Hist.Subtype1))
qplot(x = PC1, y = PC2, data = data.frame(rbind(ccle.project.ss[,1:2], ss.pca.common$x[,1:2])), colour = as.factor(c(annot$Hist.Subtype1,dbscan.ss$cluster)))


### single cell group 1 only
ss1 <- ss[,dbscan.ss$cluster == 1]
ss1.sizes <- sizes[dbscan.ss$cluster == 1]
ss1.pca <- prcomp(t(ss1))
qplot(x = PC1, y = PC2, data = data.frame(ss1.pca$x[,1:2]), main = 'PCA of single cell non-macrophage/t-cell', size = ss1.sizes)

ss1.common <- ss1[common,]
ss1.pca.common <- prcomp(t(ss1.common))
ss1.project.ccle <- scale(t(ss1.common), ccle.pca.common$center, ccle.pca.common$scale) %*% ccle.pca.common$rotation 
ccle.project.ss1 <- scale(t(ccle.common), ss1.pca.common$center, ss1.pca.common$scale) %*% ss1.pca.common$rotation 
qplot(x = PC1, y = PC2, data = data.frame(ss1.project.ccle[,1:2]), size = ss1.sizes)
qplot(x = PC1, y = PC2, data = data.frame(rbind(ss1.project.ccle[,1:2], ccle.pca.common$x[,1:2])), colour = as.factor(c(dbscan.ss1$cluster,annot$Hist.Subtype1)))
qplot(x = PC1, y = PC2, data = data.frame(ccle.project.ss1[,1:2]), colour = as.factor(annot$Hist.Subtype1))
qplot(x = PC1, y = PC2, data = data.frame(rbind(ccle.project.ss1[,1:2], ss1.pca.common$x[,1:2])), colour = as.factor(c(annot$Hist.Subtype1,dbscan.ss1$cluster)))

pc1 <- data.frame(gene = rownames(ss1.pca$rotation), PC1 = ss1.pca$rotation[,1])
pc1 <- pc1[order(pc1$PC1),]
pc2 <- data.frame(gene = rownames(ss1.pca$rotation), PC2 = ss1.pca$rotation[,2])
pc2 <- pc2[order(pc2$PC2),]

write.table(pc1, file = '/Volumes/data0/users/dyao/single_cell_rnaseq_group1_pc1.rnk', quote = F, sep = '\t', row.names = F)
write.table(pc2, file = '/Volumes/data0/users/dyao/single_cell_rnaseq_group1_pc2.rnk', quote = F, sep = '\t', row.names = F)



### fixing stuff
temp <- read.delim('ttest_phillips_ICAon_vs_INAon_adj.rnk')
temp$pvalues <- -temp$pvalues
temp <- temp[order(-temp$pvalues),]
temp <- randomize_ties(temp)
write.table(temp, 'ttest_phillips_INAon_vs_ICAon.rnk', quote = F, row.names = F, sep = '\t')

t <- read.delim('etabm389_small_intestine_NE_vs_normal.rnk')
t$pvalues <- log_pvalue_to_adj_log_pvalues(t$pvalues, base = exp(1))
write.table(t, 'etabm389_small_intestine_NE_vs_normal_adj.rnk', quote = F, row.names = F, sep = '\t')


total <- gsea_obtain_top_genesets_across_multiple_gseas('/Volumes/data0/users/dyao/lung_prostate_small_cell_rank_files_GSEA/', by = 'NES')
top.pos <- rownames(total)[1:100]
top.neg <- rev(rownames(total)[(nrow(total)-100):nrow(total)])
gsea_heatmap_v3(top.pos, 'Top_up_genesets', '/Volumes/data0/users/dyao/lung_prostate_small_cell_rank_files_GSEA/')
gsea_heatmap_v3(top.neg, 'Top_down_genesets', '/Volumes/data0/users/dyao/lung_prostate_small_cell_rank_files_GSEA/')



### ZHANG 2015
library(GEOquery)
zhang <- getGEO(GEO = 'GSE66187', destdir = getwd())
zhang <- zhang[[1]]
annot <- pData(zhang)
annot <- annot[,c('title','characteristics_ch2.3')]
annot <- data.frame(paste0(annot$title, annot$characteristics_ch2.3), row.names = rownames(annot))

labels <- fData(zhang)
labels <- labels[,'GENE_SYMBOL',drop = F]
zhang.data <- exprs(zhang)
zhang.data <- data.frame(gene = labels$GENE_SYMBOL, zhang.data)
zhang.data <- zhang.data[!(zhang.data$gene) == '',]
zhang.data <- na.omit(zhang.data)
zhang.data <- data.frame(collapse_data(zhang.data, s.start = 2, 'maxavg', group = 'gene'))

types <- factor(as.character(sapply(annot$paste0.annot.title..annot.characteristics_ch2.3., function(x) if(grepl('xenograft.*ne group: CHGA negative group$',x)) "xenograft.CRPC" else if(grepl('xenograft.*ne group: CHGA positive/ SYP positive/ ',x)) "xenograft.NEPC" else if(grepl('metastasis.*ne group: CHGA negative group$',x)) "meta.CRPC" else if(grepl('metastasis.*ne group: CHGA positive/ SYP positive/ ',x)) "meta.NEPC" else 'na')),levels=c('xenograft.CRPC','xenograft.NEPC','meta.CRPC','meta.NEPC','na'))
log_pvalues_from_ttest(zhang.data, types == 'xenograft.NEPC', types == 'xenograft.CRPC', remove.zeroes = T, output = 'zhang2015_xenograft_NEPC_vs_CRPC.rnk')
log_pvalues_from_ttest(zhang.data, types == 'meta.NEPC', types == 'meta.CRPC', remove.zeroes = T, output = 'zhang2015_metastasis_NEPC_vs_CRPC.rnk')



### ewing's sarcoma
rabbit <- read.csv('~/Downloads/GSE73610_ews-msc-sailfish.csv')
pvalues = NULL
for (i in 1:nrow(rabbit)) {
  pvalue = -log(rabbit$adjusted.p.value[i])
  stat = rabbit$log2FoldChange[i]
  if (stat < 0) {
    pvalue <- -pvalue
  }
  pvalues <- c(pvalues, pvalue)
}
pvalues <- data.frame(gene = rabbit$symbol, pvalues = pvalues)
pvalues <- pvalues[order(-pvalues$pvalues),]
pvalues <- na.omit(pvalues)
write.table(pvalues, 'Rabbitt_Ewings_sarcoma_vs_mesenchymal_stem_cell.rnk', quote = F, sep = '\t', row.names = F)


###brca basal vs luminal
brca <- as.data.frame(fread('/Volumes/data0/users/dyao/TCGA_genes_norm/BRCA_expected_count_gene_norm.txt', stringsAsFactors = F))
rownames(brca) <- brca$sample
brca <- brca[,-1]
brca <- brca[,!grepl('-11$',colnames(brca))]

annot <- read.delim('~/Dropbox/Doug/BRCA.datafreeze.20120227.txt', stringsAsFactors = F)
basal <- annot[annot$PAM50 == 'Basal',]
basal <- strtrim(basal$Sample, 15)
basal <- intersect(basal, colnames(brca))
luminal <- annot[annot$PAM50 == 'LumA' | annot$PAM50 == 'LumB',]
luminal <- strtrim(luminal$Sample, 15)
luminal <- intersect(luminal, colnames(brca))

pvalues <- log_pvalues_from_ttest(brca, basal, luminal, remove.zeroes = T)
write.table(pvalues, 'TCGA_BRCA_basal_vs_luminal.rnk', quote = F, sep = '\t', row.names = F)


#### blca papillary vs basal
blca <- as.data.frame(fread('/Volumes/data0/users/dyao/TCGA_genes_norm/BLCA_expected_count_gene_norm.txt', stringsAsFactors = F))
rownames(blca) <- blca$sample
blca <- blca[,-1]
blca <- blca[,!grepl('-11$',colnames(blca))]

annot <- read.delim('~/Downloads/BLCA_cluster-assign-k4.tsv', stringsAsFactors = F, header = F, sep = ' ')
basal <- annot[annot$V2 == 3,]
basal <- strtrim(basal$V1, 15)
basal <- intersect(basal, colnames(blca))
papillary <- annot[annot$V2 == 1,]
papillary <- strtrim(papillary$V1, 15)
papillary <- intersect(papillary, colnames(blca))

pvalues <- log_pvalues_from_ttest(blca, basal, papillary, remove.zeroes = T)
write.table(pvalues, 'TCGA_BLCA_basal_vs_papillary.rnk', quote = F, sep = '\t', row.names = F)




### small intestine NE
library(dplyr)
library(org.Hs.eg.db)
library(limma)
source('~/Dropbox/Doug/Sanaz/collapse_data.R')

nsc <- read.delim('~/Downloads/E-GEOD-6272-processed-data-1628897961.txt', stringsAsFactors = F)
nsc <- nsc[-1,]
nsc[2:ncol(nsc)] <- sapply(nsc[2:ncol(nsc)], as.numeric)
nsc.annot <- read.delim('~/Downloads/A-AFFY-33.adf.txt', stringsAsFactors = F, skip = 86, header = F)
s <- AnnotationDbi::select(org.Hs.eg.db, keys= keys(org.Hs.eg.db), columns=c("SYMBOL","ENTREZID"))
genes <- s$SYMBOL[match(nsc.annot$V6, s$ENTREZID)]
nsc.annot <- data.frame(probe = nsc.annot$V1, gene = genes)
nsc.annot <- na.omit(nsc.annot)

nsc <- nsc[nsc[,1] %in% nsc.annot$probe,]
nsc.annot <- nsc.annot[match(nsc[,1], nsc.annot$probe),]
nsc[,1] <- nsc.annot$gene
nsc <- as.data.frame(collapse_data(nsc, type = 'maxavg', group = 'Scan.REF'))
nsc <- nsc[,grepl('_SI',colnames(nsc))]
nsc <- clean_dat(nsc)


sc <- read.delim('~/Downloads/E-TABM-389-processed-data-1604615423.txt', stringsAsFactors = F, check.names = F)
sc <- sc[,c(1,5,11,17,23,29,35,41,47,53,59,65,71)]
sc <- sc[-1,]
sc[2:ncol(sc)] <- sapply(sc[2:ncol(sc)], as.numeric)

sc.annot <- read.delim('~/Downloads/A-AFFY-44.adf.txt', stringsAsFactors = F, skip = 83, header = F)
genes <- s$SYMBOL[match(sc.annot$V6, s$ENTREZID)]
sc.annot <- data.frame(probe = sc.annot$V1, gene = genes)
sc.annot <- na.omit(sc.annot)

sc <- sc[sc[,1] %in% sc.annot$probe,]
sc.annot <- sc.annot[match(sc[,1], sc.annot$probe),]
sc[,1] <- sc.annot$gene
colnames(sc) <- make.names(colnames(sc))
sc <- as.data.frame(collapse_data(sc, type = 'maxavg', group = 'Scan.REF'))
sc <- clean_dat(sc)


pvalues <- log_pvalues_from_ttest(sc, grepl('L|P',colnames(sc)), grepl('N|S',colnames(sc)), remove.zeroes = T)
write.table(pvalues, 'etabm389_small_intestine_NE_vs_normal.rnk', quote = F, row.names = F, sep = '\t')

### subclasses
files = list.files(pattern = 'rnk')
for (file in files) {
  temp <- read.delim(file)
  temp[,2] <- log_pvalue_to_adj_log_pvalues(temp[,2], exp(1))
  write.table(temp, paste0(gsub('\\.rnk','',file), '_adj.rnk'), quote = F, sep = '\t', row.names = F)
}

surface <- read.delim('~/Downloads/Combined Surface Gene List HUGO.txt', header = F, stringsAsFactors = F)
surface <- surface$V1
kinase <- read.delim('~/Dropbox/Doug/Prostate_Lung_Comparisons/uniprot_kinases.txt', header = F, stringsAsFactors = F)
kinase <- kinase$V1
ncrna <- read.delim('~/Dropbox/Doug/Prostate_Lung_Comparisons/gencode_ncrnas_transcript_id_names_types.txt', stringsAsFactors = F)
ncrna.t <- ncrna$transcript.name
lincrna <- ncrna[ncrna$biotype == 'lincRNA',]
lincrna <- lincrna$transcript.name


for (file in files) {
  temp <- read.delim(file)
  temp <- temp[temp$gene %in% lincrna,]
  write.table(temp, paste0(gsub('\\.rnk','',file), '_lincrnas.rnk'), quote = F, sep = '\t', row.names = F)
}


### Philips Nymc cmyc
data.raw <- read.delim('~/Downloads/092216 for Nikov2.txt', row.names = 1)
data <- quant_norm(data.raw)
data <- clean_dat(data)
write.table(data, 'Phillips_data_upper_quartile_normalized.txt', sep = '\t', col.names = NA, quote = F)
pvalue <- log_pvalues_from_ttest(data, grepl('INAon', colnames(data)), grepl('INAoff', colnames(data)), remove.zeroes = T)
pvalue <- randomize_ties(pvalue)
write.table(pvalue, 'ttest_phillips_ICAon_vs_ICAoff.rnk', quote = F, row.names = F, sep = '\t')

# c <- read.delim('~/Dropbox/Doug/do_RRHO_heatmap_for_Niko/ttest_phillips_ICAon_vs_ICAoff.rnk', stringsAsFactors = F, row.names = 1)
# dup.c <- c[duplicated(c$pvalues),]
# write.table(dup.c, 'test.txt', quote = F, sep = '\t')
# n <- read.delim('~/Dropbox/Doug/do_RRHO_heatmap_for_Niko/ttest_phillips_INAon_vs_INAoff.rnk', stringsAsFactors = F, row.names = 1)
# dup.n <- n[duplicated(n$pvalues),]
# write.table(dup.n, 'testn.txt', quote = F, sep = '\t')
# 
# common <- intersect(rownames(c), rownames(n))
# c <- c[rownames(c) %in% common,,drop = F]
# n <- n[rownames(n) %in% common,,drop = F]
# 
# t <- cbind(data[,grepl('ICAon', colnames(data))], data[,grepl('ICAoff', colnames(data))], data[,grepl('INAon', colnames(data))], data[,grepl('INAoff', colnames(data))])
# t <- t[rownames(tog),]
# write.table(t, 'test4.txt', quote = F, col.names = NA, sep = '\t')
# 
# s <- cbind(data.raw[,grepl('ICAon', colnames(data.raw))], data.raw[,grepl('ICAoff', colnames(data.raw))], data.raw[,grepl('INAon', colnames(data.raw))], data.raw[,grepl('INAoff', colnames(data.raw))])
# s <- s[rownames(tog),]
# 
# c <- cbind(c, rank = 1:nrow(c))
# n <- cbind(n, rank.n = 1:nrow(n))
# n <- n[rownames(c),]
# 
# tog <- cbind(c,n)
# write.table(tog, 'test3.txt', quote = F, sep = '\t', col.names = NA)
# plot(tog$rank, tog$rank.n)

### BELTRAN 2016

source('~/Dropbox/Doug/generate_RRHO_plots_all_cancers.R')

### transcript
data <- as.data.frame(fread('~/Dropbox/Datasets_Small_Cell/Beltran.2016/Beltran_2016_rsem_transcripts_upper_norm_counts.txt', stringsAsFactors = F))
txs <- data$gene
data <- data[,-1]
data <- sapply(data, as.numeric)
rownames(data) <- txs
data <- na.omit(data)

data <- data[rowSums(data) > 0,]
data <- clean_dat(data)

pvalue <- log_pvalues_from_ttest(data, grepl('_N',colnames(data)), grepl('_C',colnames(data)), remove.zeroes = T)
pvalue <- randomize_ties(pvalue)
write.table(pvalue, file = 'ttest_Beltran2016_NEPC_vs_CRPC_transcripts.rnk',quote = F, sep = '\t', row.names = F)


### gene
data <- read.delim('~/Dropbox/Datasets_Small_Cell/Beltran.2016/Beltran_2016_rsem_genes_upper_norm_counts.txt', stringsAsFactors = F, row.names = 1)
txs <- rownames(data)
data <- sapply(data, as.numeric)
rownames(data) <- txs
data <- na.omit(data)

data <- data[rowSums(data) > 0,] 
data <- clean_dat(data)

pvalue <- log_pvalues_from_ttest(data, grepl('_N',colnames(data)), grepl('_C',colnames(data)), remove.zeroes = T)
pvalue <- randomize_ties(pvalue)
write.table(pvalue, file = 'ttest_Beltran2016_NEPC_vs_CRPC_genes.rnk',quote = F, sep = '\t', row.names = F)

##### CCLE
data <- read.delim('~/Dropbox/Doug/Prostate_Lung_Comparisons/CCLE_RNAseq_raw_gene_count_matrix_lung_prostate.txt', stringsAsFactors = F, row.names = 2)
data <- data[,!grepl('PROSTATE', colnames(data))]
data <- data[2:ncol(data)]
data <- data[rowSums(data) > 0,] 
rownames(data) <- fix_date_gene(rownames(data))
data <- quant_norm(data)
data <- clean_dat(data)


annot <- read.delim('~/Dropbox/Doug/Prostate_Lung_Comparisons/CCLE_sample_info_file_2012-10-18.txt', stringsAsFactors = F)
lp <- annot[annot$CCLE.name %in% colnames(data),]
lp <- lp[match(colnames(data),lp$CCLE.name),]

type <- factor(as.character(sapply(lp$Hist.Subtype1, function(x) if(grepl("^small_cell",x)) "small_cell" else if(grepl("^adenocarcinoma",x)) "adeno" else if(grepl("^squamous",x)) "squamous" else 'other')),levels=c('small_cell','adeno','squamous','other'))

pvalue <- log_pvalues_from_ttest(data, grepl('^small_cell',lp$Hist.Subtype1), grepl('^adenocarcinoma|^squamous',lp$Hist.Subtype1), remove.zeroes = T)
pvalue <- randomize_ties(pvalue)
write.table(pvalue, 'ttest_CCLE_SC_vs_adenosquamous_genes.rnk', sep = '\t', quote = F, row.names = F)


### TCGA luad + George SCLC
annot <- as.data.frame(fread('~/Dropbox/Doug/Prostate_Lung_Comparisons/gencode_all_transcript_id_names_types.txt', stringsAsFactors = F))
annot2 <- as.data.frame(fread('~/Dropbox/Doug/Prostate_Lung_Comparisons/ensembl_all_transcript_id_names_types.txt', stringsAsFactors = F))

luad <- as.data.frame(fread('~/LUAD_expected_count_transcript.txt', stringsAsFactors = F))
rownames(luad) <- luad$sample
rownames(luad) <- gsub('\\..*', '', rownames(luad))
luad <- luad[,-1]
luad <- 2^(luad) - 1
luad <- quant_norm(luad)
luad <- clean_dat(luad)
luad <- luad[,!grepl('-11$',colnames(luad))]

##transcripts
george <- read.delim('~/Dropbox/Datasets_Small_Cell/George.2015/George_2015_rsem_genes_upper_norm_counts.txt', stringsAsFactors = F, row.names = 1)
txs <- rownames(george)
george <- sapply(george, as.numeric)
rownames(george) <- txs
george <- na.omit(george)
george <- clean_dat(george)


intx <- annot2[annot2$transcript.id %in% rownames(luad),]
intx <- intx[match(rownames(luad), intx$transcript.id),]
rownames(luad) <- intx$transcript.name
luad <- luad[!is.na(rownames(luad)),]

common <- intersect_all(rownames(luad),rownames(lusc),rownames(george))
comb.data <- data.frame(gene = common, luad[common,],lusc[common,],george[common,])

write.table(comb.data, 'george_TCGA_lung.txt', quote = F, sep = '\t', row.names = F)

pvalues <- log_pvalues_from_ttest(comb.data, colnames(george), gsub('-','\\.',colnames(luad)), remove.zeroes = T)
write.table(pvalues, 'ttest_george_SCLC_vs_TCGA_LUAD_transcripts.rnk', quote = F, sep = '\t', row.names = F)


#### Jiang SCLC + TCGA lung


jiang <- read.delim('~/Dropbox/Datasets_Small_Cell/Jiang.2016/Jiang_2016_rsem_genes_upper_norm_counts.txt', stringsAsFactors = F, row.names = 1)
txs <- rownames(jiang)
jiang <- sapply(jiang, as.numeric)
rownames(jiang) <- txs
jiang <- na.omit(jiang)
jiang <- clean_dat(jiang)


common <-intersect_all(rownames(luad),rownames(lusc), rownames(jiang))
comb.data <- data.frame(gene = common, luad[common,],lusc[common,],jiang[common,])

write.table(comb.data, 'jiang_TCGA_lung.txt', quote = F, sep = '\t', row.names = F)
pvalues <- log_pvalues_from_ttest(comb.data, colnames(jiang), gsub('-','\\.',colnames(luad)), remove.zeroes = T)
write.table(pvalues, 'ttest_jiang_SCLC_vs_TCGA_LUAD_transcripts.rnk', quote = F, sep = '\t', row.names = F)



### TCGA LUAD vs LUSC
luad <- as.data.frame(fread('/Volumes/data0/users/dyao/TCGA_genes_raw/LUAD_expected_count_gene.txt', stringsAsFactors = F))
rownames(luad) <- luad$sample
luad <- luad[,-1]
luad <- luad[,!grepl('-11$',colnames(luad))]

lusc <- as.data.frame(fread('/Volumes/data0/users/dyao/TCGA_genes_raw/LUSC_expected_count_gene.txt', stringsAsFactors = F))
rownames(lusc) <- lusc$sample
lusc <- lusc[,-1]
lusc <- lusc[,!grepl('-11$',colnames(lusc))]

prad <- as.data.frame(fread('/Volumes/data0/users/dyao/TCGA_genes_norm/prad_expected_count_gene_norm.txt', stringsAsFactors = F))
rownames(prad) <- prad$sample
prad <- prad[,-1]
prad <- prad[,!grepl('-11$',colnames(prad))]

comb.data <- cbind(lusc, prad)
pvalues <- log_pvalues_from_ttest(comb.data, colnames(prad), colnames(lusc), remove.zeroes = T)
pvalues <- randomize_ties(pvalues)
write.table(pvalues, 'ttest_TCGA_PRAD_vs_LUSC.rnk', quote = F, sep = '\t', row.names = F)


### WCDT 
wcdt <- read.csv('~/Dropbox/Datasets_Small_Cell/WCDT/SU2C WCDT Batch 1-8 Gene Norm Log2.csv', row.names = 1, stringsAsFactors = F)


annot <- read.csv('~/Dropbox/Datasets_Small_Cell/WCDT/WCDT_histology_June2015.csv', stringsAsFactors = F)
sc.annot <- annot[2:7,1]
adeno.annot <- annot[25:46,1]
pvalues <- log_pvalues_from_ttest(as.matrix(wcdt), gsub('-','\\.',sc.annot), gsub('-','\\.',adeno.annot), remove.zeroes = T)
pvalues <- randomize_ties(pvalues)
write.table(pvalues, 'ttest_wcdt_sc_vs_adeno_genes.rnk', quote = F, sep = '\t', row.names = F)
