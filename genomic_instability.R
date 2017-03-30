### genomic_instability.R
# Perform data processing for various datasets related to genomic instability


### processing myc GSE77356
data <- read.delim('~/Downloads/GSE77356_RNAseq_MYC_depletion_logFC.txt')

library(dplyr)
library(org.Hs.eg.db)
library(limma)
s <- AnnotationDbi::select(org.Hs.eg.db, keys= keys(org.Hs.eg.db), columns=c("SYMBOL","ENSEMBL"))
m <- match(data$ensembl_gene_id, s$ENSEMBL)
genes <- s$SYMBOL[m]
data$ensembl_gene_id <- genes
data <- na.omit(data)

data <- as.data.frame(data %>% group_by(ensembl_gene_id) %>% filter(min_rank(PValue) == 1))

pvalues = NULL
for (i in 1:nrow(data)) {
  d <- data[i,]
  pvalue <- -log(d$PValue)
  if (d$logFC < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}

data <- data.frame(gene = data$ensembl_gene_id, pvalue = pvalues)
data <- data[order(-data$pvalue),]
write.table(data, '/Volumes/data0/users/dyao/LORENZIN_myc_KD.rnk', quote = F, sep = '\t', row.names = F)


### titration
data <- read.delim('~/Downloads/GSE77356_RNAseq_titration_logFC.txt')
s <- AnnotationDbi::select(org.Hs.eg.db, keys= keys(org.Hs.eg.db), columns=c("SYMBOL","ENSEMBL"))
m <- match(data$ensembl_gene_id, s$ENSEMBL)
genes <- s$SYMBOL[m]
data$ensembl_gene_id <- genes
data <- na.omit(data)

temp <- data[,c(1,11,12)]
temp <- as.data.frame(temp %>% group_by(ensembl_gene_id) %>% filter(min_rank(Pvalue_DOX_1) == 1))
temp <- temp[!duplicated(temp$ensembl_gene_id),]
pvalues = NULL
for (i in 1:nrow(temp)) {
  d <- temp[i,]
  pvalue <- -log(d[,3])
  if (d[,2] < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = temp$ensembl_gene_id, pvalue = pvalues)
pvalues <- pvalues[order(-pvalues$pvalue),]
write.table(pvalues, '/Volumes/data0/users/dyao/LORENZIN_myc_1.rnk', quote = F, sep = '\t', row.names = F)



### processing myc GSE68219
myc <- read.csv('~/Downloads/1P Metabolism Data.csv')
myc$gene <- gsub('-chr.*','',myc$gene)
pvalues = NULL
for (i in 1:nrow(myc)) {
  d <- myc[i,]
  pvalue <- -log(d$pValueMC_M)
  if (d$fcMC_MK < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}

myc <- data.frame(gene = myc$gene, pvalue = pvalues)
myc <- myc[order(-myc$pvalue),]
write.table(myc, '/Volumes/data0/users/dyao/LY_myc_KD.rnk', quote = F, sep = '\t', row.names = F)



### processing myc GSE66250
cd44 <- read.delim('~/Downloads/GSE66250_cd44high_dox_vs_etoh.txt', stringsAsFactors = F)
cd44 <- cd44[!cd44$hgnc_symbol == '',]
cd44 <- as.data.frame(cd44 %>% group_by(hgnc_symbol) %>% filter(min_rank(PValue) == 1))

pvalues = NULL
for (i in 1:nrow(cd44)) {
  d <- cd44[i,]
  pvalue <- -log(d$PValue)
  if (d$logFC < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = cd44$hgnc_symbol, pvalue = pvalues)
pvalues <- pvalues[order(-pvalues$pvalue),]
write.table(pvalues, '/Volumes/data0/users/dyao/VON_EYSS_CD44_high_myc.rnk', quote = F, sep = '\t', row.names = F)

imecs <- read.delim('~/Downloads/GSE66250_imecs_myc_dox_vs_vector_dox.txt', stringsAsFactors = F)
imecs <- imecs[!imecs$hgnc_symbol == '',]
imecs <- as.data.frame(imecs %>% group_by(hgnc_symbol) %>% filter(min_rank(PValue) == 1))

pvalues = NULL
for (i in 1:nrow(imecs)) {
  d <- imecs[i,]
  pvalue <- -log(d$PValue)
  if (d$logFC < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = imecs$hgnc_symbol, pvalue = pvalues)
pvalues <- pvalues[order(-pvalues$pvalue),]
write.table(pvalues, '/Volumes/data0/users/dyao/VON_EYSS_imecs_myc.rnk', quote = F, sep = '\t', row.names = F)




### rank rank myc vs instability
elkon <- read.delim('/Volumes/data0/users/dyao/ELKON_myc.rnk', stringsAsFactors = F)
walz1 <- read.delim('/Volumes/data0/users/dyao/WALZ_myc_DOX.rnk', stringsAsFactors = F)
walz2 <- read.delim('/Volumes/data0/users/dyao/WALZ_myc_VECT.rnk', stringsAsFactors = F)
ave.myc <- read.delim('/Volumes/data0/users/dyao/ave_myc_rank.rnk', stringsAsFactors = F)
ave.myc$rnk <- -ave.myc$rnk
cell1 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group.rnk', stringsAsFactors = F)
cell2 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group_v3.rnk', stringsAsFactors = F)
tumor1 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group1.rnk', stringsAsFactors = F)
tumor2 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group2.rnk', stringsAsFactors = F)

create_RRHO_file_intersect_then_rank(elkon, cell1, x_axis = 'Elkon_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(elkon, cell2, x_axis = 'Elkon_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v2', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(elkon, tumor1, x_axis = 'Elkon_Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(elkon, tumor2, x_axis = 'Elkon_Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_2', reverse = F, ranked = F, outputdir = getwd(), plot = T)

create_RRHO_file_intersect_then_rank(walz1, cell1, x_axis = 'Walz1_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(walz1, cell2, x_axis = 'Walz1_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v2', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(walz1, tumor1, x_axis = 'Walz1_Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(walz1, tumor2, x_axis = 'Walz1_Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_2', reverse = F, ranked = F, outputdir = getwd(), plot = T)

create_RRHO_file_intersect_then_rank(ave.myc, cell1, x_axis = 'Ave_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(ave.myc, cell2, x_axis = 'Ave_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v2', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(ave.myc, tumor1, x_axis = 'Ave_Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(ave.myc, tumor2, x_axis = 'Ave_Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_2', reverse = F, ranked = F, outputdir = getwd(), plot = T)


### finding different genes between myc and instability
find_diff_genes <- function(ave.myc, cell1, label1, label2) {
  myc.cell1 <- create_RRHO_file_intersect_then_rank(ave.myc, cell1, x_axis = 'Ave_Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = F)
  cell1.threshold.top <- ceiling(nrow(myc.cell1)/5)
  cell1.threshold.bot <- ceiling(nrow(myc.cell1)/5*4)
  myc.threshold1 <- ceiling(nrow(myc.cell1)/3)
  myc.threshold2 <- ceiling(nrow(myc.cell1)/3*2)
  
  notsig.myc.up.cell1 <- myc.cell1[myc.cell1$ave.myc.rank >= myc.threshold1 & myc.cell1$cell1.rank <= cell1.threshold.top,]
  downgenes <- notsig.myc.up.cell1$Unigene
  write.table(downgenes, file = paste0('Up_', label2, '_Not_up_', label1, '_genes.txt'), sep = '\n', quote = F, row.names = F, col.names = F)
              
  notsig.myc.down.cell1 <- myc.cell1[myc.cell1$ave.myc.rank <= myc.threshold2 & myc.cell1$cell1.rank >= cell1.threshold.bot,]
  downgenes <- notsig.myc.down.cell1$Unigene
  write.table(downgenes, file = paste0('Down_', label2, '_Not_down_', label1, '_genes.txt'), sep = '\n', quote = F, row.names = F, col.names = F)
  
  up.myc.notsig.cell1 <- myc.cell1[myc.cell1$ave.myc.rank <= cell1.threshold.top & myc.cell1$cell1.rank >= myc.threshold1,]
  downgenes <- up.myc.notsig.cell1$Unigene
  write.table(downgenes, file = paste0('Up_', label1, '_Not_up_', label2, '_genes.txt'), sep = '\n', quote = F, row.names = F, col.names = F)
  
  down.myc.notsig.cell1 <- myc.cell1[myc.cell1$ave.myc.rank >= cell1.threshold.bot & myc.cell1$cell1.rank <= myc.threshold2,]
  downgenes <- down.myc.notsig.cell1$Unigene
  write.table(downgenes, file = paste0('Down_', label1, '_Not_down_', label2, '_genes.txt'), sep = '\n', quote = F, row.names = F, col.names = F)
}

find_diff_genes(ave.myc, cell1, 'ave_myc', 'cell1')
find_diff_genes(ave.myc, tumor1, 'ave_myc', 'tumor1')

## gene sets genes 
gsea_lookup_gene_v3(list_genes_geneset('SECRETION'),'SECRETION_GENES')
gsea_lookup_gene_v3(list_genes_geneset('EXTRACELLULAR_MATRIX'),'ECM_GENES')

gsea_lookup_gene_v3(list_genes_geneset(list_genesets_keyword('GOLGI','all')), 'ALL_GOLGI_GENESETS_GENES')
gsea_lookup_gene_v3(list_genes_geneset(list_genesets_keyword('LIPID','all')), 'ALL_LIPID_GENESETS_GENES')

gsea_lookup_gene_v3('ERBB2', 'ERBB2')


### GETTING AVERAGE MYC sig
elkon <- read.delim('/Volumes/data0/users/dyao/myc/ELKON_myc.rnk', stringsAsFactors = F, row.names = 1)
sorted.elkon <- elkon
elkon$pvalue <- 1:nrow(elkon)
walz1 <- read.delim('/Volumes/data0/users/dyao/myc/WALZ_myc_DOX.rnk', stringsAsFactors = F, row.names = 1)
sorted.walz1 <- walz1
walz1$pvalue <- 1:nrow(walz1)
walz2 <- read.delim('/Volumes/data0/users/dyao/myc/WALZ_myc_VECT.rnk', stringsAsFactors = F, row.names = 1)
walz2$pvalue <- 1:nrow(walz2)
sorted.walz2 <- walz2
voneyss1 <- read.delim('/Volumes/data0/users/dyao/myc/VON_EYSS_CD44_high_myc.rnk', stringsAsFactors = F, row.names = 1)
voneyss1$pvalue <- 1:nrow(voneyss1)
sorted.voneyss1 <- voneyss1
voneyss2 <- read.delim('/Volumes/data0/users/dyao/myc/VON_EYSS_imecs_myc.rnk', stringsAsFactors = F)
voneyss2 <- voneyss2[!duplicated(voneyss2$gene),]
rownames(voneyss2) <- voneyss2$gene
voneyss2 <- voneyss2[,-1,drop = F]
voneyss2$pvalue <- 1:nrow(voneyss2)
sorted.voneyss2 <- voneyss2
ly.kd <- read.delim('/Volumes/data0/users/dyao/myc/LY_myc_KD.rnk', stringsAsFactors = F, row.names = 1)
ly.kd$pvalue <- 1:nrow(ly.kd)
sorted.ly.kd <- ly.kd
lorenzin.kd <- read.delim('/Volumes/data0/users/dyao/myc/LORENZIN_myc_KD.rnk', stringsAsFactors = F, row.names = 1)
lorenzin.kd$pvalue <- 1:nrow(lorenzin.kd)
sorted.lorenzin.kd <- lorenzin.kd

genes <- Reduce(intersect, list(rownames(elkon), rownames(walz1), rownames(walz2)))

elkon <- elkon[genes,,drop = F]
walz1 <- walz1[genes,,drop = F]
walz2 <- walz2[genes,,drop = F]



ave <- rowMeans(data.frame(elkon$pvalue, walz1$pvalue, walz2$pvalue))
ave <- data.frame(gene = genes, rank = ave)
ave <- ave[order(ave$rank),]
write.table(ave, file = '/Volumes/data0/users/dyao/ave_myc_rank.rnk', row.names = F, quote = F, sep = '\t')

pcs <- read.delim('~/Dropbox/Doug/CNA PCA/BRCA/BRCA.pca.data.R_prcomp_scores.PC1-5_sample_short_names.txt', stringsAsFactors = F)
pcs$Sample <- gsub('-','.',pcs$Sample)
pcs <- pcs[order(pcs$PC1),]

exp <- read.delim('/Volumes/data2/users/nbalanis/Input_Files_and_Standards/TCGA_RNASEQV2_DATA_from_tcga_assembler/BRCA_normalized_results_processed.txt', stringsAsFactors = F)
exp1 <- as.data.frame(collapse_data(exp, s.start = 3, type = 'maxavg', group = 'GeneSymbol'))

samples <- colnames(exp1)
samples <- strtrim(samples, 15)
intx <- intersect(pcs$Sample, samples)

colnames(exp1) <- strtrim(colnames(exp1), 15)
exp2 <- exp1[,intx]
exp2 <- exp2[!(rownames(exp2) == '?'),]
exp2 <- t(scale(t(exp2)))
exp2 <- na.omit(exp2)

temp <- intersect(rownames(sorted.walz1)[(nrow(sorted.walz1)-500):nrow(sorted.walz1)], rownames(exp2))
temp <- as.character(ave$gene[1:500])
toavg <- exp2[temp,]
toavg <- colSums(toavg)

setwd('/Volumes/data0/users/dyao/')
png(filename = 'BRCA_PC1_ave_myc.png', height = 700, width = 700)
plot(toavg,  main = 'BRCA CNA PC1 samples vs Average Myc Score', xlab = 'BRCA TCGA samples sorted by PC1', ylab = 'Average Myc Score', pch = 20)
y_sym2 <- movingAverage(toavg, 200, T)
lines(1:length(toavg), y_sym2, col="red", lwd = 3)
dev.off()

movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}



### MYC SIGNATURE
source('~/Dropbox/Doug/Sanaz/collapse_data.R')

data <- read.delim('~/Downloads/GSE66789_Myc_both_rounds_rnaseq_fpkms_floor.txt',stringsAsFactors = F)
data <- data[,3:7]
data <- as.data.frame(collapse_data(data, type = 'maxavg', group = 'sym'))
data <- log2(data)
data <- data[!abs(((data$Myc.rna.rep1-data$C.rna.rep1)-(data$Myc.rna.rep2-data$C.rna.rep2))) < 0.00001,]

pvalues = NULL
for (i in 1:nrow(data)) {
  ctrl <- data[i,c(1,3)]
  myc <- data[i,c(2,4)]
  tTest <- t.test(as.numeric(myc),as.numeric(ctrl), paired = T)
  pvalue <- -log(tTest$p.value)
  if (tTest$statistic < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues, pvalue)
}

pvalues <- data.frame(gene = rownames(data), pvalue = pvalues)
pvalues <- pvalues[order(-pvalues$pvalue),]
write.table(pvalues, file = '/Volumes/data0/users/dyao/myc_pvalues.rnk', row.names = F, quote = F, sep = '\t')
pvalues <- read.delim('myc_pvalues.rnk')

create_RRHO_file_intersect_then_rank(pvalues, cell1, x_axis = 'Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(pvalues, cell2, x_axis = 'Myc_overexpression_pvalues', y_axis = 'Core_cell_instability_signature_v2', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(pvalues, tumor1, x_axis = 'Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(pvalues, tumor2, x_axis = 'Myc_overexpression_pvalues', y_axis = 'Core_tumor_instability_signature_group_2', reverse = F, ranked = F, outputdir = getwd(), plot = T)


s <- AnnotationDbi::select(org.Hs.eg.db, keys= keys(org.Hs.egENSEMBL, keytype = 'ENSEMBL'), columns=c("SYMBOL","ENSEMBL"))
m2 <- match(data[,1], s$ENSEMBL)

genes <- s$SYMBOL[m2]
data <- data.frame(GENE = genes, data[,3:ncol(data)])
data <- na.omit(data)

library(dplyr)
library(org.Hs.eg.db)
library(limma)

data <- read.delim('~/Downloads/GSE44672_U2OS_shMiz1_+-Dox_RNAseq.txt',stringsAsFactors = F)
data <- read.delim('~/Downloads/GSE44672_U2OS_MycWT_MycV394D_RNAseq.txt')
genes <- intersect(data$hgnc_symbol, tumor1$GENE)
data <- data[data$hgnc_symbol %in% genes,]
data <- as.data.frame(data %>% group_by(hgnc_symbol) %>% filter(min_rank(p.value_MycWT_vs_empty) == 1))
data <- data[!duplicated(data$hgnc_symbol),]

pvalues = NULL
for (i in 1:nrow(data)) {
  pvalue <- -log(data[i,5])
  stat <- data[i,4]
  if (stat < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues, pvalue)
}
pvalues <- data.frame(gene = data$hgnc_symbol, pvalue = pvalues)
pvalues <- pvalues[order(-pvalues$pvalue),]
pvalues[1,2] <- 500

write.table(pvalues, file = '/Volumes/data0/users/dyao/WALZ_myc_VECT.rnk', quote = F, row.names = F, sep = '\t')


### MELANOMA vs instability
cell1 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group.rnk')
cell2 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group_v3.rnk')
tumor1 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group1.rnk')
tumor2 <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group2.rnk')
pc1 <- read.delim('/Volumes/data0/users/dyao/Melanoma_cell_Line_PC1_loadings.rnk')
diff <- read.delim('/Volumes/data0/users/dyao/undiff_vs_diff_melanoma.rnk')

create_RRHO_file_intersect_then_rank(diff, cell1, x_axis = 'Undifferentiated_Neural_crest_like_vs_Transitory_Melanocytic', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(diff, cell2, x_axis = 'Undifferentiated_Neural_crest_like_vs_Transitory_Melanocytic', y_axis = 'Core_cell_instability_signature_v2', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(diff, tumor1, x_axis = 'Undifferentiated_Neural_crest_like_vs_Transitory_Melanocytic', y_axis = 'Core_tumor_instability_signature_group_1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(diff, tumor2, x_axis = 'Undifferentiated_Neural_crest_like_vs_Transitory_Melanocytic', y_axis = 'Core_tumor_instability_signature_group_2', reverse = F, ranked = F, outputdir = getwd(), plot = T)

create_RRHO_file_intersect_then_rank(pc1, cell1, x_axis = 'Melanoma_PC1_loadings', y_axis = 'Core_cell_instability_signature_v1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(pc1, cell2, x_axis = 'Melanoma_PC1_loadings', y_axis = 'Core_cell_instability_signature_v2', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(pc1, tumor1, x_axis = 'Melanoma_PC1_loadings', y_axis = 'Core_tumor_instability_signature_group_1', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(pc1, tumor2, x_axis = 'Melanoma_PC1_loadings', y_axis = 'Core_tumor_instability_signature_group_2', reverse = F, ranked = F, outputdir = getwd(), plot = T)


### COMPARE GSEA genesets
compare_gsea_v2('GSEA_core_c2/1/ave_stem_cell_sig/','GSEA_core_c2/1/core_cell_group/','ave_stem','core_cell', type = 'cpg')
compare_gsea_v2('GSEA_core_c2/1/ave_stem_cell_sig/','GSEA_core_c2/1/core_cell_group/','ave_stem','core_cell', type = 'c2')
compare_gsea_v2('GSEA_core_c5/1/ave_stem_cell_sig/','GSEA_core_c5/1/core_cell_group/','ave_stem','core_cell', type = 'c5')

compare_gsea_v2('GSEA_core_c2/1/ave_stem_cell_sig/','GSEA_core_c2/1/core_tumor_group1/','ave_stem','core_tumor1', type = 'cpg')
compare_gsea_v2('GSEA_core_c2/1/ave_stem_cell_sig/','GSEA_core_c2/1/core_tumor_group1/','ave_stem','core_tumor1', type = 'c2')
compare_gsea_v2('GSEA_core_c5/1/ave_stem_cell_sig/','GSEA_core_c5/1/core_tumor_group1/','ave_stem','core_tumor1', type = 'c5')

compare_gsea_v2('GSEA_core_c2/1/ave_stem_cell_sig/','GSEA_core_c2/1/core_tumor_group2/','ave_stem','core_tumor2', type = 'cpg')
compare_gsea_v2('GSEA_core_c2/1/ave_stem_cell_sig/','GSEA_core_c2/1/core_tumor_group2/','ave_stem','core_tumor2', type = 'c2')
compare_gsea_v2('GSEA_core_c5/1/ave_stem_cell_sig/','GSEA_core_c5/1/core_tumor_group2/','ave_stem','core_tumor2', type = 'c5')

compare_gsea_v2('GSEA_core_c5/1/core_tumor_group1/','GSEA_core_c5/1/core_cell_group/','core_tumor1','core_cell', type = 'c5')

###COMBINED SUBGROUPS

tumorgroup1 <- c('OV', 'CESC', 'HNSC', 'BLCA', 'LUAD', 'LUSC', 'PAAD', 'PRAD', 'LIHC', 'COAD', 'READ', 'STAD', 'ESCA', 'BRCA')
tumorgroup1v2 <- c('HNSC', 'BLCA', 'LUAD', 'LUSC', 'PAAD', 'PRAD', 'LIHC', 'COAD', 'READ', 'STAD', 'ESCA', 'BRCA')
tumorgroup2 <- c('MESO', 'TGCT', 'SKCM', 'GBM')
tumorgroup2v2 <- c('MESO', 'TGCT', 'SKCM')
cellgroup <- c('BREAST', 'SKIN', 'PLEURA', 'UPPER_AERODIGESTIVE_TRACT')
cellgroupv2 <- c('BREAST', 'SKIN', 'PLEURA')

cellgroupv3 <- c('BREAST_cell_line_ICNA', 'BREAST_cell_line_bkpt','SKIN_cell_line_ICNA','BONE_cell_line_ICNA')

tumordir <- '/Volumes/data0/users/dyao/RANK_FILES/tumor_combined_ICNA_bkpt_ranks/'
celldir <- '/Volumes/data0/users/dyao/RANK_FILES/cell_line_combined_ICNA_bkpt_ranks/'

a <- get_ave_tumor(tumorgroup1, tumordir)
write.table(a, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group1.rnk', row.names = F, quote = F, sep = '\t')
b <- get_ave_tumor(tumorgroup2, tumordir)
write.table(b, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group2.rnk', row.names = F, quote = F, sep = '\t')
c <- get_ave_tumor(cellgroup, celldir)
write.table(c, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group.rnk', row.names = F, quote = F, sep = '\t')

d <- get_ave_tumor(tumorgroup1v2, tumordir)
write.table(d, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group1_v2.rnk', row.names = F, quote = F, sep = '\t')
e <- get_ave_tumor(tumorgroup2v2, tumordir)
write.table(e, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_tumor_group2_v2.rnk', row.names = F, quote = F, sep = '\t')
f <- get_ave_tumor(cellgroupv2, celldir)
write.table(f, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group_v2.rnk', row.names = F, quote = F, sep = '\t')

g <- get_ave_tumor(cellgroupv3, celldir)
write.table(g, file = '/Volumes/data0/users/dyao/RANK_FILES/ave_groups/core_cell_group_v3.rnk', row.names = F, quote = F, sep = '\t')

get_ave_tumor <- function(tumorgroup1, tumordir) {
  tumorfiles <- list.files(tumordir, pattern = '.rnk')
  alltumor <- NULL
  for (i in tumorgroup1) {
    files <- tumorfiles[grepl(i,tumorfiles)]
    for (j in files) {
      temp <- read.delim(paste0(tumordir, j), row.names = 1)
      if (is.null(ncol(alltumor))) {
        alltumor <- temp 
      } else {
        gene.names <- rownames(alltumor)
        temp <- temp[gene.names,, drop = F]
        alltumor <- cbind(alltumor, temp)
      }
    }
  }
  avetumor <- rowMeans(alltumor)
  avetumor <- data.frame(GENE = names(avetumor), RVALUE = avetumor)
  avetumor <- avetumor[order(-avetumor$RVALUE),]
  return (avetumor)
}

####GSEA GENES
receptor_genes <- c(c('AATK','ALK','AXL','CSF1R','DDR1','DDR2'), paste0('EPHA',1:8), c('EPHA10'), paste0('EPHB',1:4), c('EPHB6','EGFR','ERBB2','ERBB3','ERBB4','FGFR1','FGFR2','FGFR3','FGFR4','FLT1','FLT3','FLT4','IGF1R','INSR','INSRR','KDR','KIT','LMTK2','LMTK3','LTK','MERTK','MET','MST1R','MUSK','NTRK1','NTRK2','NTRK3','PDGFRA','PDGFRB','PTK7','RET','ROR1','ROR2','ROS1','RYK','STYK1','TEK','TIE1','TYRO3'))
ligand_genes <- c('HGF','FGF2','NRG1','SGCD','EDN1','EDN2','EDN3','FGF4','FGF8','EGF','IL4','FGF1','IGF2','TNF','CSF1','IL10','FGF12','FIBP','PILRB','CTSH','BTC','LTA','MMP13','TNFSF14','TGFA','TNFSF12','FGF9','IGF1','PDGFA','PDGFB','IL2','TGFB1','TGFB2','TGFB3')
gsea_lookup_gene_v2(receptor_genes, 'receptors')
gsea_lookup_gene_v2(ligand_genes, 'ligands')
gsea_lookup_gene_v2(c(paste0('SDC',1:4), paste0('GPC',1:6),'NDST1','HPSE', 'DCN','VCAN', 'PRKCA', 'PRKCB'), 'GAG_genes')
gsea_lookup_gene_v2(c('SRGN','AGRN','HSPG2','COL18A1'), 'GAG_genes_2')

gsea_lookup_gene_v2(c('ACLY','ACACA','ACACB','MLYCD','FASN','ACSS1','ACSS2','ACSL1','SCD'), 'lipid_genes')
gsea_lookup_gene_v2(c('ERN1','EIF2AK3','ATF6'), 'UPR_genes')
gsea_lookup_gene_v2(c('PRKCA','PRKCB'), 'PRKC_genes')
gsea_lookup_gene_v2(c('PTDSS1','PTDSS2'), 'PTDSS_genes')
gsea_lookup_gene_v2(c('PDGFRA','PDGFRB'), 'PDGFR_genes')

data <- read.delim('~/Downloads/er_golgi_kinase_genes.txt', header = F, stringsAsFactors = F)
data <- data$V1
gsea_lookup_gene_v2(data, 'ER_GOLGI_kinase_genes')


###GSEA CLUSTERING

setwd('/Volumes/data0/users/dyao/')
genesets <- list_genesets_keyword(c('GOLGI', 'DAMAGE'), 'all')
genesets <- genesets[!grepl('CYTOKINE', genesets)]
genesets <- genesets[!grepl('NEUROTRANSMITTER', genesets)]
genesets <- genesets[!grepl('BODY_FLUID', genesets)]
genesets <- genesets[!grepl('INTERLEUKIN', genesets)]
genesets <- genesets[!grepl('GLYCAN', genesets)]
genesets <- genesets[!grepl('GHRELIN', genesets)]
genesets <- genesets[!grepl('INACTIVATION_OF_GIP', genesets)]
genesets <- genesets[!grepl('INCRETIN', genesets)]
genesets <- genesets[!grepl('INSULIN', genesets)]
genesets <- genesets[!grepl('GLP1', genesets)]
gsea_heatmap(genesets, 'GOLGI_CORE', 'all', core = 'v1', cluster = F)
gsea_heatmap_v2(genesets, 'DAMAGE_TREATED_GOLGI', '/Volumes/data0/users/dyao/DAMAGE_GSEA/')
gsea_heatmap_v2(list_genesets_keyword('SECRETION', 'all'), 'DAMAGE_TREATED_SECRETION', '/Volumes/data0/users/dyao/DAMAGE_GSEA/')

gsea_heatmap(list_genesets_keyword('RECEPTOR', 'all'), 'RECEPTOR', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('GOLGI', 'all'), 'GOLGI', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('MEMBRANE', 'all'), 'MEMBRANE', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('LIPID', 'all'), 'LIPID', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('APOPT', 'all'), 'APOPTOSIS', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('DEATH', 'all'), 'DEATH', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('FATTY', 'all'), 'FATTY_ACID', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('_ER_', '^ER', 'ENDOPLASMIC'), 'all'), 'ENDOPLASMIC_RETICULUM', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('METABOL'), 'all'), 'METABOLISM', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('AKT','all'), 'AKT', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('TNF', '_NF_'), 'all'), 'NF', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('CHOLESTEROL', 'all'), 'CHOLESTEROL', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('SYNTHESIS', 'all'), 'SYNTHESIS', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('SECRETION', 'all'), 'SECRETION', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('NOTCH', 'all'), 'NOTCH', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('GROWTH', 'all'), 'GROWTH', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('GLYCOLYSIS', 'all'), 'GLYCOLYSIS', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('PROTEIN', 'all'), 'PROTEIN', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('METAL', 'all'), 'METAL', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('WOUND', 'all'), 'WOUND', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('ANGIOGENESIS', 'all'), 'ANGIOGENESIS', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('B_CELL','T_CELL'), 'all'), 'T_CELL_B_CELL', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('EPITHELIAL', 'MESENCHYMAL'), 'all'), 'EPITHELIAL_MESENCHYMAL', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('CYTO', 'all'), 'CYTO', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('TGF', 'all'), 'TGF', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('EGF','FGF'), 'all'), 'EGF_FGF', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('LIPASE','PROTEASE'), 'all'), 'LIPASE_PROTEASE', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('GTP'), 'all'), 'GTP', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('INTEGRIN'), 'all'), 'INTEGRIN', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('MTOR'), 'all'), 'MTOR', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('MATRIX', 'ECM'), 'all'), 'MATRIX', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword('JAK', 'all'), 'JAK', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('UBIQUI','PROTEASOM'), 'all'), 'UBIQUITIN_PROTEASOME', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('MITOCH'), 'all'), 'MITOCHONDRIA', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('MITOT','MITOS'), 'all'), 'MITOSIS', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('DAMAGE'), 'all'), 'DAMAGE', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('CYCLE'), 'all'), 'CYCLE', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('BIND'), 'all'), 'BIND', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('DNA'), 'all'), 'DNA', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('RNA'), 'all'), 'RNA', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('IMMUNE'), 'all'), 'IMMUNE', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('REPAIR'), 'all'), 'REPAIR', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('TRANSCRIPT'), 'all'), 'TRANSCRIPTION', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('TRANSLAT'), 'all'), 'TRANSLATION', 'all', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('KINASE'), 'all'), 'KINASE', 'all', core = 'v1', cluster = F)


gsea_heatmap(list_genesets_keyword(c('RECEPTOR', 'MEMBRANE'), 'all'), 'RECEPTOR_MEMBRANE', 'all')
gsea_heatmap(list_genesets_keyword(c('RECEPTOR', 'MEMBRANE'), 'all'), 'RECEPTOR_MEMBRANE_CORE', 'all', core = T, cluster = F)
gsea_heatmap(list_genesets_keyword(c('STEM'), 'cpg'), 'STEM_CORE', 'c2', core = 'v1', cluster = F)
gsea_heatmap(list_genesets_keyword(c('_MYC'), 'cpg'), 'MYC_CORE', 'c2', core = 'v1', cluster = F)



###RRHO CLUSTERING
data.p <- read.delim('/Volumes/data0/users/dyao/RANK_FILES/combined_all/RRHO_Matrix.txt', row.names = 1)
data.p <- data.p[!grepl('stem',colnames(data.p)),!grepl('stem',colnames(data.p))]
data.p <- data.p[!grepl('avg_tumor',colnames(data.p)),!grepl('avg_tumor',colnames(data.p))]
data.p <- data.p[!grepl('cell_line_cell_line',colnames(data.p)),!grepl('cell_line_cell_line',colnames(data.p))]

names <- colnames(data.p)


finalnames <- gsub('_tumor_bkpt_rvalues','',names)
finalnames <- gsub('_tumor_ICNA_rvalues','',finalnames)
finalnames <- gsub('_cell_line_bkpt','',finalnames)
finalnames <- gsub('_cell_line_ICNA','',finalnames)
finalnames <- unique(finalnames)

alltypes <- list()

nums <- read.delim('~/Dropbox/Doug/ICNA_BKPT/num_samples_per_cell_line.txt', header = F)

alltypes[['kidney']] <- list('KIDNEY', c('KIRC', 'KIRP', 'KICH'))
alltypes[['bile']] <- list('BILIARY_TRACT', 'CHOL')
alltypes[['uterus']] <- list('ENDOMETRIUM', c('UCEC','UCS','CESC'))
alltypes[['skin']] <- list('SKIN', 'SKCM')
alltypes[['liver']] <- list('LIVER', 'LIHC')
alltypes[['lung']] <- list('LUNG', c('LUAD', 'LUSC'))
alltypes[['breast']] <- list('BREAST', 'BRCA')
alltypes[['blood']] <- list('HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', c('LAML','DLBC'))
alltypes[['bone']] <- list('BONE', 'LAML')
alltypes[['tissue']] <- list('SOFT_TISSUE', 'SARC')
alltypes[['prostate']] <- list('PROSTATE', 'PRAD')
alltypes[['brain']] <- list('CENTRAL_NERVOUS_SYSTEM', c('GBM', 'LGG'))
alltypes[['largeint']] <- list('LARGE_INTESTINE', c('COAD', 'READ'))
alltypes[['headneck']] <- list('UPPER_AERODIGESTIVE_TRACT', 'HNSC')
alltypes[['pleura']] <- list('PLEURA', 'MESO')
alltypes[['stomach']] <- list('STOMACH', 'STAD')
alltypes[['thyroid']] <- list('THYROID', 'THCA')
alltypes[['ovary']] <- list('OVARY', 'OV')
alltypes[['pancreas']] <- list('PANCREAS', 'PAAD')
alltypes[['esophagus']] <- list('OESOPHAGUS', 'ESCA')
alltypes[['bladder']] <- list('URINARY_TRACT', 'BLCA')

alltypes2 <- c('KIDNEY (n=33)','BILIARY TRACT (n=8)','ENDOMETRIUM (n=27)','SKIN (n=61)','LIVER (n=27)','LUNG (n=180)','BREAST (n=59)','BLOOD (n=178)','BONE (n=27)','SOFT TISSUE (n=21)','PROSTATE (n=7)','CNS (n=62)','LARGE INTESTINE (n=58)','UPPER AERO TRACT (n=30)','PLEURA (n=10)','STOMACH (n=37)','THYROID (n=10)','OVARY (n=50)','PANCREAS (n=44)','ESOPHAGUS (n=26)','URINARY TRACT (n=21)')



total.combined <- NULL
for (i in 1:length(alltypes)) {
  cell <- alltypes[[i]][[1]]
  tumor <- alltypes[[i]][[2]]
  cell.bkpt <- paste0(cell, '_cell_line_bkpt')
  cell.icna <- paste0(cell, '_cell_line_ICNA')
  tumor.bkpt <- paste0(tumor, '_tumor_bkpt_rvalues')
  tumor.icna <- paste0(tumor, '_tumor_ICNA_rvalues')
  
  if (length(tumor) == 1) {
    icna.data <- data.p[cell.icna, tumor.icna]
    names(icna.data) <- paste0(alltypes2[i], '_vs_', tumor, '_ICNA')
    bkpt.data <- data.p[cell.bkpt, tumor.bkpt]
    names(bkpt.data) <- paste0(alltypes2[i], '_vs_', tumor, '_bkpt')
    combined <- c(icna.data, bkpt.data)
    total.combined <- c(total.combined, combined)
  }
  
  else if (length(tumor) > 1) {
    for (j in 1:length(tumor)) {
      icna.data <- data.p[cell.icna, tumor.icna[j]]
      names(icna.data) <- paste0(alltypes2[i], '_vs_', tumor[j], '_ICNA')
      bkpt.data <- data.p[cell.bkpt, tumor.bkpt[j]]
      names(bkpt.data) <- paste0(alltypes2[i], '_vs_', tumor[j], '_bkpt')
      combined <- c(icna.data, bkpt.data)
      total.combined <- c(total.combined, combined)
    }
  }
}
total.combined <- data.frame(TYPE = names(total.combined), RRHOvalue = total.combined)
total.combined$group <- factor(gsub('_.*','',total.combined$TYPE))

total.combined$group <- reorder(total.combined$group, -total.combined$RRHOvalue, mean, order = T)
total.combined <- total.combined[order(total.combined$group),]
rownames(total.combined) <- gsub('_', ' ', rownames(total.combined))


par(oma = c(8,0,0,0))
barplot(total.combined$RRHOvalue, names.arg = gsub('_',' ',total.combined$TYPE), las = 2, main = 'Cell lines vs cancer RRHO all cancers', cex.names = 0.7, ylab = 'RRHO value', ylim = c(0,800))



i = 1
j = 1
rrhos <- NULL
while (i <= length(names)) {
  group1 <- data.p[names[i],]
  group1 <- group1[,!(grepl(finalnames[j],colnames(group1)))]
  group2 <- data.p[names[i+1],]
  group2 <- group2[,!(grepl(finalnames[j],colnames(group2)))]
  
  pvalue <- -log(t.test(as.numeric(group1),as.numeric(group2), paired = T)$p.value)
  rrhos <- c(rrhos, pvalue)
  i = i+2
  j = j+1
}
rrhos <- data.frame(TYPE = finalnames, PVALUE = rrhos)
rrhos <- rrhos[order(rrhos$PVALUE),]
rrhos$TYPE <- factor(rrhos$TYPE, levels = rrhos$TYPE)


i = 1
rrhos <- NULL
while (i <= length(names)) {
  rrhos <- c(rrhos, data.p[names[i],names[i+1]])
  i = i+2
}



for (i in 1:length(nums$V1)) {
  finalnames[finalnames == nums$V1[i]] <- paste0(finalnames[finalnames == nums$V1[i]], ' (n=', nums$V2[i],')')
}

finalnames <- gsub('_', ' ', finalnames)
finalnames[finalnames == 'CENTRAL NERVOUS SYSTEM (n=62)'] <- 'CNS (n=62)'
finalnames[finalnames == 'HAEMATOPOIETIC AND LYMPHOID TISSUE (n=178)'] <- 'BLOOD (n=178)'
finalnames[finalnames == 'UPPER AERODIGESTIVE TRACT (n=30)'] <- 'UPPER AERO TRACT (n=30)'


rrhos <- data.frame(TYPE = finalnames, RRHOvalue = rrhos)
rrhos <- rrhos[order(-rrhos$RRHOvalue),]
rrhos$TYPE <- factor(rrhos$TYPE, levels = rrhos$TYPE)
write.table(rrhos, file = 'RRHO_ICNA_vs_bkpt_all_cancers.txt', sep = '\t', row.names = F, quote = F)

par(oma = c(5,0,0,0))
barplot(rrhos$RRHOvalue, names.arg = rrhos$TYPE, las = 2, main = 'ICNA vs bkpt RRHO all cancers', cex.names = 0.7, ylab = 'RRHO value', ylim = c(0,6000))

data.p <- as.matrix(data.p)







data <- read.delim('RANK_FILES/stem_cell_rank_files/ave_stem_cell_sig.rnk')
data <- data[order(data$pvalue),]

setwd('/Volumes/data0/users/dyao/')
compare_gsea_all_cross_comparisons('GSEA_stem_sig', 'GSEA_cell_line_bkpt', 'simple')
compare_gsea_all_cross_comparisons('GSEA_stem_sig', 'GSEA_cell_line_ICNA', 'simple')
compare_gsea_all_cross_comparisons('GSEA_stem_sig', 'GSEA_tumor_bkpt', 'simple')
compare_gsea_all_cross_comparisons('GSEA_stem_sig', 'GSEA_tumor_ICNA', 'simple')

### GSEA COMPARISONS
compare_gsea('melanoma_PC1_ICNA', 'melanoma_cell_lines_ICNA', 'PC1', 'cell_line_ICNA')
compare_gsea('melanoma_PC1_ICNA', 'SKCM_tumor_ICNA_rvalues', 'PC1', 'tumor_ICNA')
compare_gsea('SKCM_tumor_ICNA_rvalues', 'melanoma_cell_lines_ICNA', 'tumor_ICNA', 'cell_line_ICNA')
compare_gsea('ca125negvpos', 'OV_tumor_ICNA_rvalues', 'CA125negvpos', 'OV_ICNA')

newdir <- 'tumor_ICNA_individual_genes'
dir.create(newdir)
setwd(newdir)
directory <- '/Volumes/data0/users/dyao/GSEA/GSEA_tumor_ICNA/'
files <- list.files(directory)
for (i in files) {
  get_genes(paste0(directory, '/', i), i, 18990)
}
setwd('..')

newdir <- 'stem_individual_genes'
dir.create(newdir)
setwd(newdir)
directory <- '/Volumes/data0/users/dyao/GSEA/GSEA_stem_sig/'
files <- list.files(directory)
lengths <- c(15518, 24693, 18802, 19104, 18653)
for (i in 1:length(files)) {
  get_genes(paste0(directory, '/', files[i]), files[i], lengths[i])
}
setwd('..')

setwd('/Volumes/data0/users/dyao/GSEA_genes/')
compare_genes_all_cross_comparisons('stem_individual_genes', 'avg_tumor_cell_line_genes')



### small cell PCA copy number


loadings <- read.delim('CCLE_CNV_annotation_by_genes_cleaned_cleaned_prcomp_loadings_PC1-5.gct')
pc1 <- data.frame(gene = loadings$Gene, PC1 = loadings$PC1)
pc2 <- data.frame(gene = loadings$Gene, PC2 = loadings$PC2)
write.table(pc1, file = 'small_cell_PC1_loadings.rnk', quote = F, row.names = F, sep = '\t')
write.table(pc2, file = 'small_cell_PC2_loadings.rnk', quote = F, row.names = F, sep = '\t')

t <- pc2[order(-pc2$PC2),]
create_RRHO_file_intersect_then_rank(ave2, pc1, x_axis = 'Average stem signature pvalue', y_axis = 'Small cell PC1 loadings', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(ave2, pc2, x_axis = 'Average stem signature pvalue', y_axis = 'Small cell PC2 loadings', outputdir = getwd(), reverse = F, ranked = F)

gsea_obtain_top_single('~/gsea_home/output/jul19/my_analysis.GseaPreranked.1468913051834/', 'small_cell_PC2')
gsea_obtain_top_single('~/gsea_home/output/jul19/my_analysis.GseaPreranked.1468913777919/', 'small_cell_PC1')

### small cell PCA expression
loadings <- read.delim('CCLE_RMA_normalized_lung_prostate_prcomp_loadings.txt')
pc1 <- data.frame(gene = loadings$Loading, PC1 = loadings$PC1)
pc2 <- data.frame(gene = loadings$Loading, PC2 = loadings$PC2)
write.table(pc1, file = 'small_cell_PC1_loadings_expression.rnk', quote = F, row.names = F, sep = '\t')
write.table(pc2, file = 'small_cell_PC2_loadings_expression.rnk', quote = F, row.names = F, sep = '\t')

t <- pc2[order(-pc2$PC2),]
create_RRHO_file_intersect_then_rank(ave2, pc1, x_axis = 'Average stem signature pvalue', y_axis = 'Small cell PC1 loadings expression', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(ave2, pc2, x_axis = 'Average stem signature pvalue', y_axis = 'Small cell PC2 loadings expression', outputdir = getwd(), reverse = F, ranked = F)

gsea_obtain_top_single('~/gsea_home/output/jul19/my_analysis.GseaPreranked.1468913051834/', 'small_cell_PC2')
gsea_obtain_top_single('~/gsea_home/output/jul19/my_analysis.GseaPreranked.1468913777919/', 'small_cell_PC1')


### Plotting stem vs instability scatterplots
plots_directory <- 'TCGA_stem_scores_vs_instability_scores_scatterplots'
plots <- list.files(plots_directory)
filtered_plots <- plots[grepl('filtered', plots)]
unfiltered_plots <- plots[!grepl('filtered', plots)]
pdf(paste0('TCGA_stem_scores_vs_instability_scores1.pdf'), width = 8.5, height = 11)
a <- layout(matrix(c(1,1,2,2,1,1,2,2,3,4,5,6,7,8,9,10,11,12,13,14), 5, 4, byrow = TRUE), respect = T)

for (i in 1:14) {
par(mar = c(0,0,0,0))
jj <- readPNG(paste0(plots_directory, '/', filtered_plots[i]))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
rasterImage(jj,0,0,1,1)
}

dev.off()

pdf(paste0('TCGA_stem_scores_vs_instability_scores2.pdf'), width = 8.5, height = 11)
a <- layout(matrix(c(1:20), 5, 4, byrow = TRUE), respect = T)

for (i in 15:22) {
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(plots_directory, '/', filtered_plots[i]))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  rasterImage(jj,0,0,1,1)
}

dev.off()





### more stem v instability stuff
stem <- read.delim('stem_cell_rank_files/ave_stem_cell_sig.rnk')
stem <- stem[order(-stem$pvalue),]
rownames(stem) <- 1:nrow(stem)
tumor1 <- read.delim('tumor_ICNA_ranks/avg_tumor_ICNA_rvalues.rnk')
tumor2 <- read.delim('tumor_bkpt_ranks/avg_tumor_bkpt_rvalues.rnk')
inx <- intersect(stem$gene[1:500], tumor2$gene[1:500])
write.table(inx, 'intersecting_genes.txt', row.names = F, col.names = F, quote = F, sep = '\t')


### small cell ranks 
small <- read.delim('./small_cell_lung/small_cell_lung_RNA_pearson_correlation.txt')
small <- small[!small$X == '',]
non_small <- read.delim('./small_cell_lung/non_small_cell_lung_RNA_pearson_correlation.txt')
non_small <- small[!small$X == '',]
small_bkpt <- data.frame(gene = small$X, bkpt = small$bkpt_vs_RNA_rvalue)
small_icna <- data.frame(gene = small$X, ICNA = small$ICNA_vs_RNA_rvalue)
non_small_bkpt <- data.frame(gene = non_small$X, bkpt = non_small$bkpt_vs_RNA_rvalue)
non_small_icna <- data.frame(gene = non_small$X, bkpt = non_small$ICNA_vs_RNA_rvalue)

write.table(small_bkpt, file = 'small_cell_lung_bkpt.rnk', row.names = F, quote = F, sep = '\t')
write.table(small_icna, file = 'small_cell_lung_ICNA.rnk', row.names = F, quote = F, sep = '\t')
write.table(non_small_bkpt, file = 'non_small_cell_lung_bkpt.rnk', row.names = F, quote = F, sep = '\t')
write.table(non_small_icna, file = 'non_small_cell_lung_ICNA.rnk', row.names = F, quote = F, sep = '\t')



### PLOTTING stuff
library(png)
library(jpeg)

move_ave_to_front <- function(l, name) {
  averages <- l[grepl(name, l)]
  l <- l[!(l %in% averages)]
  l <- c(averages, l)
  return (l)
}

tables_directory <- '/Volumes/data0/users/dyao/combined_stem_vs_cell_line_bkpt_ICNA'
cell_tables <- list.files(tables_directory)
ave_cell_tables <- cell_tables[grepl('^ave',cell_tables)]
ave_cell_tables <- move_ave_to_front(ave_cell_tables, 'average')

tables_gsea_directory <- '/Volumes/data0/users/dyao/combined_GSEA_top_cell_line_bkpt_ICNA/'
gsea_tables <- list.files(tables_gsea_directory)
gsea_tables <- move_ave_to_front(gsea_tables, 'average')

plots_directory <- '/Volumes/data0/users/dyao/RRHO_cell_line/Combined_bkpt_ICNA_average_stem_cell_lines'
ave_cell_scatters_and_heat <- list.files(plots_directory)
ave_cell_scatters <- ave_cell_scatters_and_heat[grepl('^RankScatter', ave_cell_scatters_and_heat)]
ave_cell_scatters <- move_ave_to_front(ave_cell_scatters, 'average')
ave_cell_heat <- ave_cell_scatters_and_heat[grepl('^RRHO', ave_cell_scatters_and_heat)]
ave_cell_heat <- move_ave_to_front(ave_cell_heat, 'average')

cancer_names <- gsub('RankScatterStem_cell_pvalues_VS_','',ave_cell_scatters)
cancer_names <- gsub('_bkpt_score.jpg','',cancer_names)
cancer_names <- gsub('_ICNA_score.jpg','',cancer_names)
cancer_names <- gsub('_',' ',cancer_names)

stem_tables_directory <- '/Volumes/data0/users/dyao/GSEA_top_stem/'
stem_tables <- list.files(stem_tables_directory)
stem_names <- c('Average stem signature', 'Breast stem signature', 'Colon EPHB2 stem signature', 'Colon PTK7 stem signature', 'Prostate CD49F Benign stem signature')

remainder = length(cancer_names)%%3
i = 1
j = 1
while (i < (length(cancer_names) - remainder)) {
  pdf(paste0('Ave_stem_signature_vs_cell_line_instability_signatures_' ,i, '.pdf'), width = 8.5, height = 11)
  a <- layout(matrix(c(1,1,1,2,2,2,3,4,5,6,6,6,7,8,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,17,18,19,19,19,20,21,22), 14, 3, byrow = TRUE), widths = rep(c(1,1,4),9), heights = c(0.2,0.2,1,0.2,1,0.2,1,0.2,1,0.2,0.2,1,0.2,1), respect = T)

  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Ave stem vs ', cancer_names[i]), padj = -1, side = 1)
    
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('Top overlapping gene sets by breakpoints', side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_scatters[i]),native=TRUE)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_heat[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_directory,'/',ave_cell_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.03,0.07,1.03,0.93)
  i = i + 1
    
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('Top overlapping gene sets by ICNA score', side = 1, padj = -1, cex = 0.8)
    
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_scatters[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
    
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_heat[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
    
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_directory,'/',ave_cell_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.03,0.07,1.03,0.93)
  i = i+1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ',cancer_names[j],' by breakpoints'), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_gsea_directory,'/',gsea_tables[j]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
  j = j+1
    
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ',cancer_names[j],' by ICNA'), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_gsea_directory,'/',gsea_tables[j]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
  j = j+1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Ave stem vs ', cancer_names[i]), padj = -1, side = 1)
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('Top overlapping gene sets by breakpoints', side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_scatters[i]),native=TRUE)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_heat[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_directory,'/',ave_cell_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.03,0.07,1.03,0.93)
  i = i + 1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('Top overlapping gene sets by ICNA score', side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_scatters[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_heat[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_directory,'/',ave_cell_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.03,0.07,1.03,0.93)
  i = i+1
  
  
  dev.off()
  
  
  pdf(paste0('Ave_stem_signature_vs_cell_line_instability_signatures_' ,i, '.pdf'), width = 8.5, height = 11)
  a <- layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,8,9,10,10,10,11,12,13,14,14,14,15,15,15,16,16,16,17,17,17), 13, 3, byrow = TRUE), widths = rep(c(1,1,4),9), heights = c(0.2,1,0.2,1,0.2,0.2,1,0.2,1,0.2,1,0.2,1), respect = T)

  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ',cancer_names[j],' by breakpoints'), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_gsea_directory,'/',gsea_tables[j]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
  j = j+1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ',cancer_names[j],' by ICNA'), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_gsea_directory,'/',gsea_tables[j]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
  j = j+1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Ave stem vs ', cancer_names[i]), padj = -1, side = 1)
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('Top overlapping gene sets by breakpoints', side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_scatters[i]),native=TRUE)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_heat[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_directory,'/',ave_cell_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.03,0.07,1.03,0.93)
  i = i + 1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('Top overlapping gene sets by ICNA score', side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_scatters[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readJPEG(paste0(plots_directory, '/', ave_cell_heat[i]),native=TRUE)
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,0,0,1,1)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_directory,'/',ave_cell_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.03,0.07,1.03,0.93)
  i = i+1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ',cancer_names[j],' by breakpoints'), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_gsea_directory,'/',gsea_tables[j]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
  j = j+1
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ',cancer_names[j],' by ICNA'), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(tables_gsea_directory,'/',gsea_tables[j]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
  j = j+1
  dev.off()
}


pdf(paste0('Stem_signatures_',i, '.pdf'), width = 8.5, height = 11)
a <- layout(matrix(1:10, 10, 1, byrow = TRUE), widths = rep(6,10), heights = rep(c(0.2,1),5), respect = T)
for (i in 1:5) {
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext(paste0('Most significant up/downregulated gene sets for ', stem_names[i]), side = 1, padj = -1, cex = 0.8)
  
  par(mar = c(0,0,0,0))
  jj <- readPNG(paste0(stem_tables_directory,'/',stem_tables[i]))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(jj,-0.05,0.05,1.05,0.95)
}

dev.off()






#### Cell line instability data from CCLE
source('~/Dropbox/Doug/Sanaz/collapse_data.R')
data <- read.table('~/Documents/CCLE_Expression_Entrez_2012-09-29.gct', header=TRUE, row.names=NULL, sep='\t', check.names = FALSE,stringsAsFactors = FALSE)
data <- data[,-1]
data <- data.frame(collapse_data(data, type = 'maxavg', group = 'Description'))

pheno.data <- read.delim('~/Documents/CCLE_sample_info_file_2012-10-18.txt')
cell_bkpt <- read.table(file = '~/Dropbox/Doug/ICNA_BKPT/cell_line_bkpts_sorted_by_bkpts_per_chrom.txt',header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
cell_bkpt <- cell_bkpt[!duplicated(cell_bkpt[,1]),]
rownames(cell_bkpt) <- cell_bkpt[,1]
cell_bkpt <- cell_bkpt[,-1]

count = NULL
for (i in 1:length(levels(groups))) {
  temp <- data[,groups == levels(groups)[i]]
  tcount <- length(intersect(colnames(temp), rownames(cell_bkpt)))
  count <- c(count, tcount)
  #compute_icna_bkpt(prefix = as.character(levels(groups)[i]), log2_flag = F, temp = temp, breakpoints = cell_bkpt)
}

names(count) <- as.character(levels(groups))
write.table(count, file = 'num_samples_per_cell_line.txt', quote = F, col.names = F, sep = '\t')



#### Cell line instability data from CCLE w/ small cell info
setwd('~/Dropbox/Doug/ICNA_BKPT/')
data <- read.delim('~/Documents/CCLE_sample_info_file_2012-10-18.txt')
lung <- data[grepl('LUNG', data$CCLE.name),]
small_cell_lung <- lung[grepl('^small_cell', lung$Hist.Subtype1),]
non_small_cell_lung <- lung[!grepl('^small_cell', lung$Hist.Subtype1),]

exp_data <- read.table('~/Documents/CCLE_Expression_Entrez_2012-09-29.gct', header=TRUE, row.names=NULL, sep='\t', check.names = FALSE,stringsAsFactors = FALSE)
non_small_cell_lung <- non_small_cell_lung[non_small_cell_lung$CCLE.name %in% colnames(exp_data),]
exp_data <- exp_data[,-1]
exp_data <- data.frame(collapse_data(exp_data, type = 'maxavg', group = 'Description'))
sc_lung_group <- exp_data[,as.character(small_cell_lung$CCLE.name)]
non_sc_lung_group <- exp_data[,as.character(non_small_cell_lung$CCLE.name)]
cell_bkpt <- read.table(file = '~/Dropbox/Doug/ICNA_BKPT/cell_line_bkpts_sorted_by_bkpts_per_chrom.txt',header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
cell_bkpt <- cell_bkpt[!duplicated(cell_bkpt[,1]),]
rownames(cell_bkpt) <- cell_bkpt[,1]
cell_bkpt <- cell_bkpt[,-1]

compute_icna_bkpt(prefix = 'small_cell_lung', log2_flag = F, temp = sc_lung_group, breakpoints = cell_bkpt)
compute_icna_bkpt(prefix = 'non_small_cell_lung', log2_flag = F, temp = non_sc_lung_group, breakpoints = cell_bkpt)

setwd('~/Dropbox/Doug/ICNA_BKPT')
breast <- read.delim('breast_illumhumanWG-6_v3_MaSC-Luminal_maxavg_Pval.tsv', check.names = F)
prostate <- read.delim('prostate_Quant_CD49fHi.Benign-CD49fLo.Benign_maxavg_Pval.tsv', check.names = F)
colon1 <- read.delim('colon_RMA_EPHB2HighMed-EPHB2LowNeg_maxavg_Pval.tsv', check.names = F)
colon2 <- read.delim('colon_RMA_PTK7HighMed-PTK7LowNeg_maxavg_Pval.tsv', check.names = F)

ave <- Reduce(intersect, list(rownames(breast), rownames(prostate), rownames(colon1), rownames(colon2)))
ave1 <- data.frame(breast = breast[ave,], prostate = prostate[ave,], colon1 = colon1[ave,], colon2 = colon2[ave,], row.names = ave)
ave2 <- data.frame(pvalue = rowMeans(ave1))
ave2 <- data.frame(gene = rownames(ave2), pvalue = ave2$pvalue)

breast <- data.frame(gene = rownames(breast), pvalue = breast$logP.Value)
prostate <- data.frame(gene = rownames(prostate), pvalue = prostate$logP.Value)
colon1 <- data.frame(gene = rownames(colon1), pvalue = colon1$logP.Value)
colon2 <- data.frame(gene = rownames(colon2), pvalue = colon2$logP.Value)

small_cell_data <- read.delim('small_cell_lung_RNA_pearson_correlation.txt')
small_cell_data_icna <- data.frame(gene = small_cell_data$X, ICNA = small_cell_data$ICNA_vs_RNA_rvalue)
small_cell_data_bkpt <- data.frame(gene = small_cell_data$X, bkpt = small_cell_data$bkpt_vs_RNA_rvalue)

non_small_cell_data <- read.delim('non_small_cell_lung_RNA_pearson_correlation.txt')
non_small_cell_data_icna <- data.frame(gene = non_small_cell_data$X, ICNA = non_small_cell_data$ICNA_vs_RNA_rvalue)
non_small_cell_data_bkpt <- data.frame(gene = non_small_cell_data$X, bkpt = non_small_cell_data$bkpt_vs_RNA_rvalue)


create_RRHO_file_intersect_then_rank(ave2, small_cell_data_icna, x_axis = 'Average stem signature pvalue', y_axis = 'Small cell lung ICNA score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(ave2, small_cell_data_bkpt, x_axis = 'Average stem signature pvalue', y_axis = 'Small cell lung breakpoints', outputdir = getwd(), reverse = F, ranked = F)

create_RRHO_file_intersect_then_rank(ave2, non_small_cell_data_icna, x_axis = 'Average stem signature pvalue', y_axis = 'Non small cell lung ICNA score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(ave2, non_small_cell_data_bkpt, x_axis = 'Average stem signature pvalue', y_axis = 'Non small cell lung breakpoints', outputdir = getwd(), reverse = F, ranked = F)

create_RRHO_file_intersect_then_rank(small_cell_data_icna, non_small_cell_data_icna, x_axis = 'Small cell lung ICNA score', y_axis = 'Non small cell lung ICNA score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(small_cell_data_bkpt, non_small_cell_data_bkpt, x_axis = 'Small cell lung breakpoints', y_axis = 'Non small cell lung breakpoints', outputdir = getwd(), reverse = F, ranked = F)


### Boxplots non small cell vs small cell
setwd('~/Dropbox/Doug/ICNA_BKPT/')
data <- read.delim('~/Documents/CCLE_sample_info_file_2012-10-18.txt', stringsAsFactors = F)
lung <- data[grepl('LUNG', data$CCLE.name),]
small_cell_lung <- lung[grepl('^small_cell', lung$Hist.Subtype1),]$CCLE.name
squamous <- lung[grepl('squamous', lung$Hist.Subtype1),]$CCLE.name
adeno <- lung[grepl('adenocarcinoma', lung$Hist.Subtype1),]$CCLE.name
non_small_cell_lung <- lung[!grepl('^small_cell', lung$Hist.Subtype1),]$CCLE.name

cell_bkpt <- read.table(file = '~/Dropbox/Doug/ICNA_BKPT/cell_line_bkpts_sorted_by_bkpts_per_chrom.txt',header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
cell_bkpt <- cell_bkpt[!duplicated(cell_bkpt[,1]),]
rownames(cell_bkpt) <- cell_bkpt[,1]
cell_bkpt <- cell_bkpt[,-1]

small_cell_icna <- na.omit(cell_bkpt[small_cell_lung,'ICNA_score'])
small_cell_bkpt <- na.omit(cell_bkpt[small_cell_lung,'bkpt_samp'])
non_small_cell_icna <- na.omit(cell_bkpt[non_small_cell_lung,'ICNA_score'])
non_small_cell_bkpt <- na.omit(cell_bkpt[non_small_cell_lung,'bkpt_samp'])
squamous_icna <- na.omit(cell_bkpt[squamous,'ICNA_score'])
squamous_bkpt <- na.omit(cell_bkpt[squamous,'bkpt_samp'])
adeno_icna <- na.omit(cell_bkpt[adeno,'ICNA_score'])
adeno_bkpt <- na.omit(cell_bkpt[adeno,'bkpt_samp'])

icna <- list(small_cell_icna, non_small_cell_icna)
names(icna) <- c('Small cell', 'Non small cell')

bkpt <- list(small_cell_bkpt, non_small_cell_bkpt)
names(bkpt) <- c('Small cell', 'Non small cell')

icna2 <- list(small_cell_icna, squamous_icna, adeno_icna)
names(icna2) <- c('Small cell', 'Squamous', 'Adenocarcinoma')

bkpt2 <- list(small_cell_bkpt, squamous_bkpt, adeno_bkpt)
names(bkpt2) <- c('Small cell', 'Squamous', 'Adenocarcinoma')

boxplot(icna, main = 'Range of ICNA scores for lung cancers', ylab = 'ICNA score', outline = F)
boxplot(bkpt, main = 'Range of breakpoints for lung cancers', ylab = 'Breakpoints', outline = F)

boxplot(icna2, main = 'Range of ICNA scores for lung cancers', ylab = 'ICNA score', outline = F)
boxplot(bkpt2, main = 'Range of breakpoints for lung cancers', ylab = 'Breakpoints', outline = F)



### COADREAD 
bkpt_coad <- read.table(file = '~/Dropbox/Doug/COADREAD_normalized_results_processed.txt', header = T, row.names = NULL, sep = '\t', check.names = F, stringsAsFactors = F)
bkpt_coad <- bkpt_coad[!duplicated(bkpt_coad[,1]),]
rownames(bkpt_coad) <- bkpt_coad[,1]
bkpt_coad <- bkpt_coad[,-1]
bkpt_coad <- bkpt_coad[,-1]

bkpts <- read.table(file = '~/Dropbox/Doug/COAD_bkpts_sorted_by_bkpts_per_chrom.txt',header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
bkpts <- bkpts[!duplicated(bkpts[,1]),]
rownames(bkpts) <- bkpts[,1]
bkpts <- bkpts[,-1]

a <- strtrim(rownames(bkpts), nchar(rownames(bkpts))-13)
bkpts <- data.frame(bkpts[!(duplicated(a)),], row.names = a[!(duplicated(a))])
compute_icna_bkpt(prefix = 'COADREAD', log2_flag = T, temp = bkpt_coad, breakpoints = bkpts)

coadread <- read.delim('COADREAD_log2_RNA_pearson_correlation.txt')
coadread_icna <- data.frame(GENE = coadread$X, ICNA = coadread$ICNA_vs_RNA_rvalue)
coadread_bkpt <- data.frame(GENE = coadread$X, bkpt = coadread$bkpt_vs_RNA_rvalue)
write.table(coadread_icna, file = 'COADREAD_tumor_ICNA_rvalues.rnk', sep = '\t', quote = F, row.names = F)
write.table(coadread_bkpt, file = 'COADREAD_tumor_bkpt_rvalues.rnk', sep = '\t', quote = F, row.names = F)


### PLOTTING TCGA INSTABILITY VS STEMNESS
setwd('~/Dropbox/Doug/ICNA_BKPT/TCGA_stem_signature/')
file_list <- list.files()
names <- sub('_.*','',file_list)
dataset_stem = list() 
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i], row.names = 1)
  dataset_stem[[names[i]]] <- temp
}

setwd('~/Dropbox/Doug/ICNA_BKPT/bkpt.icna files')
file_list2 <- list.files()
file_list2 <- file_list2[!(file_list2 == 'others')]
names2 <- sub('_.*','',file_list2)
dataset_ICNA = list()
for (i in 1:length(file_list2)) {
  temp <- read.delim(file_list2[i], row.names = 1)
  temp <- data.frame(ICNA_score = temp$ICNA_score, row.names = rownames(temp))
  dataset_ICNA[[names2[i]]] <- temp
}

dataset_bkpt = list()
for (i in 1:length(file_list2)) {
  temp <- read.delim(file_list2[i], row.names = 1)
  temp <- data.frame(bkpt = temp$bkpt_samp, row.names = rownames(temp))
  dataset_bkpt[[names2[i]]] <- temp
}

setwd('~/Dropbox/Doug/ICNA_BKPT/TCGA_stem_scores_vs_instability_scores_scatterplots/')

total_stem = NULL
total_ICNA = NULL
total_bkpt = NULL
intx <- intersect(names(dataset_stem), names(dataset_ICNA))

for (i in intx) {
  total_stem <- rbind(total_stem, dataset_stem[[i]])
  total_ICNA <- rbind(total_ICNA, dataset_ICNA[[i]])
  total_bkpt <- rbind(total_bkpt, dataset_bkpt[[i]])
}

dataset_stem[['all']] <- total_stem
dataset_ICNA[['all']] <- total_ICNA
dataset_bkpt[['all']] <- total_bkpt
intx <- intersect(names(dataset_stem), names(dataset_ICNA))


for (i in intx) {
  temp_stem <- dataset_stem[[i]]
  temp_icna <- dataset_ICNA[[i]]
  temp_bkpt <- dataset_bkpt[[i]]
  a <- strtrim(rownames(temp_stem), nchar(rownames(temp_stem))-9)
  b <- strtrim(rownames(temp_icna), nchar(rownames(temp_icna))-9)
  c <- strtrim(rownames(temp_bkpt), nchar(rownames(temp_bkpt))-9)
  temp_stem <- data.frame(scores = temp_stem[!(duplicated(a)),], row.names = a[!(duplicated(a))])
  temp_icna <- data.frame(scores = temp_icna[!(duplicated(b)),], row.names = b[!(duplicated(b))])
  temp_bkpt <- data.frame(scores = temp_bkpt[!(duplicated(c)),], row.names = c[!(duplicated(c))])
  
  temp_stem_filtered <- data.frame(scores = temp_stem[!temp_stem$scores %in% boxplot.stats(temp_stem$scores)$out,], row.names = rownames(temp_stem)[!temp_stem$scores %in% boxplot.stats(temp_stem$scores)$out])
  temp_icna_filtered <- data.frame(scores = temp_icna[!temp_icna$scores %in% boxplot.stats(temp_icna$scores)$out,], row.names = rownames(temp_icna)[!temp_icna$scores %in% boxplot.stats(temp_icna$scores)$out])
  temp_bkpt_filtered <- data.frame(scores = temp_bkpt[!temp_bkpt$scores %in% boxplot.stats(temp_bkpt$scores)$out,], row.names = rownames(temp_bkpt)[!temp_bkpt$scores %in% boxplot.stats(temp_bkpt$scores)$out])
  
  temp_intersect_icna <- intersect(rownames(temp_stem), rownames(temp_icna))
  temp_combined_icna <- data.frame(scores = temp_stem[temp_intersect_icna,], icna = temp_icna[temp_intersect_icna,], row.names = temp_intersect_icna)
  png(file = paste0(i, '_TCGA_stem_score_vs_ICNA_score.png'))
  plot(temp_combined_icna$scores, temp_combined_icna$icna, main = paste0(i,' TCGA stem score vs ICNA score'), xlab = 'Stem score', ylab = 'ICNA score', pch = 20)
  abline(fit <- lm(temp_combined_icna$icna~temp_combined_icna$scores), col="black")
  legend("topright", bty="n", legend=paste("R =", format(sqrt(summary(fit)$adj.r.squared), digits=4)))
  dev.off()
  
  temp_intersect_icna <- intersect(rownames(temp_stem_filtered), rownames(temp_icna_filtered))
  temp_combined_icna <- data.frame(scores = temp_stem_filtered[temp_intersect_icna,], icna = temp_icna_filtered[temp_intersect_icna,], row.names = temp_intersect_icna)
  png(file = paste0(i, '_TCGA_stem_score_vs_ICNA_score_filtered.png'))
  plot(temp_combined_icna$scores, temp_combined_icna$icna, main = paste0(i,' TCGA stem score vs ICNA score filtered'), xlab = 'Stem score', ylab = 'ICNA score', pch = 20)
  abline(fit <- lm(temp_combined_icna$icna~temp_combined_icna$scores), col="black")
  legend("topright", bty="n", legend=paste("R =", format(sqrt(summary(fit)$adj.r.squared), digits=4)))
  dev.off()
  
  temp_intersect_bkpt <- intersect(rownames(temp_stem), rownames(temp_bkpt))
  temp_combined_bkpt <- data.frame(scores = temp_stem[temp_intersect_bkpt,], bkpt = temp_bkpt[temp_intersect_bkpt,], row.names = temp_intersect_bkpt)
  png(file = paste0(i, '_TCGA_stem_score_vs_bkpt_score.png'))
  plot(temp_combined_bkpt$scores, temp_combined_bkpt$bkpt, main = paste0(i,' TCGA stem score vs Breakpoints'), xlab = 'Stem score', ylab = 'Breakpoints', pch = 20)
  abline(fit <- lm(temp_combined_bkpt$bkpt~temp_combined_bkpt$scores), col="black")
  legend("topright", bty="n", legend=paste("R =", format(sqrt(summary(fit)$adj.r.squared), digits=4)))
  dev.off()
  
  temp_intersect_bkpt <- intersect(rownames(temp_stem_filtered), rownames(temp_bkpt_filtered))
  temp_combined_bkpt <- data.frame(scores = temp_stem_filtered[temp_intersect_bkpt,], bkpt = temp_bkpt_filtered[temp_intersect_bkpt,], row.names = temp_intersect_bkpt)
  png(file = paste0(i, '_TCGA_stem_score_vs_bkpt_score_filtered.png'))
  plot(temp_combined_bkpt$scores, temp_combined_bkpt$bkpt, main = paste0(i,' TCGA stem score vs Breakpoints filtered'), xlab = 'Stem score', ylab = 'Breakpoints', pch = 20)
  abline(fit <- lm(temp_combined_bkpt$bkpt~temp_combined_bkpt$scores), col="black")
  legend("topright", bty="n", legend=paste("R =", format(sqrt(summary(fit)$adj.r.squared), digits=4)))
  dev.off()
}





### CELL LINE INSTABILITY vs CELL LINE STEM
setwd('~/Dropbox/Doug/ICNA_BKPT/CCLE_cell_line_ICNA/')
file_list <- list.files()
names <- strtrim(file_list, nchar(file_list) - 28)
dataset_cell_icna = list()
dataset_cell_bkpt = list()
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i], row.names = 1)
  temp <- temp[!(rownames(temp) == ''),]
  temp1 <- data.frame(ICNA = temp$ICNA_vs_RNA_rvalue, row.names = rownames(temp))
  temp2 <- data.frame(bkpt = temp$bkpt_vs_RNA_rvalue, row.names = rownames(temp))
  dataset_cell_icna[[names[i]]] <- temp1
  dataset_cell_bkpt[[names[i]]] <- temp2
}

for (i in 1:length(dataset_cell_icna)) {
  temp <- data.frame(GENE = rownames(dataset_cell_icna[[i]]), ICNA = dataset_cell_icna[[i]]$ICNA)
  write.table(temp, file = paste0(names(dataset_cell_icna)[i], '_cell_line_ICNA.rnk'), quote = F, sep = '\t', row.names = F)
}

for (i in 1:length(dataset_cell_bkpt)) {
  temp <- data.frame(GENE = rownames(dataset_cell_bkpt[[i]]), bkpt = dataset_cell_bkpt[[i]]$bkpt)
  write.table(temp, file = paste0(names(dataset_cell_bkpt)[i], '_cell_line_bkpt.rnk'), quote = F, sep = '\t', row.names = F)
}

setwd('~/Dropbox/Doug/ICNA_BKPT')
breast <- read.delim('breast_illumhumanWG-6_v3_MaSC-Luminal_maxavg_Pval.tsv', check.names = F)
prostate <- read.delim('prostate_Quant_CD49fHi.Benign-CD49fLo.Benign_maxavg_Pval.tsv', check.names = F)
colon1 <- read.delim('colon_RMA_EPHB2HighMed-EPHB2LowNeg_maxavg_Pval.tsv', check.names = F)
colon2 <- read.delim('colon_RMA_PTK7HighMed-PTK7LowNeg_maxavg_Pval.tsv', check.names = F)

ave <- Reduce(intersect, list(rownames(breast), rownames(prostate), rownames(colon1), rownames(colon2)))
ave1 <- data.frame(breast = breast[ave,], prostate = prostate[ave,], colon1 = colon1[ave,], colon2 = colon2[ave,], row.names = ave)
ave2 <- data.frame(pvalue = rowMeans(ave1))
ave2 <- data.frame(gene = rownames(ave2), pvalue = ave2$pvalue)
write.table(ave2, file = 'ave_stem_cell_sig.rnk', quote = F, row.names = F, sep = '\t')

breast <- data.frame(gene = rownames(breast), pvalue = breast$logP.Value)
prostate <- data.frame(gene = rownames(prostate), pvalue = prostate$logP.Value)
colon1 <- data.frame(gene = rownames(colon1), pvalue = colon1$logP.Value)
colon2 <- data.frame(gene = rownames(colon2), pvalue = colon2$logP.Value)


do_RRHO_single(breast, c('BREAST','PROSTATE','LARGE_INTESTINE'), 'pvalues', outfile.prefix = 'breast_stem_cell_sig_cell_line', type = 'cell')
do_RRHO_single(prostate, c('BREAST','PROSTATE','LARGE_INTESTINE'), 'pvalues', outfile.prefix = 'prostate_stem_cell_sig_cell_line', type = 'cell')
do_RRHO_single(colon1, c('BREAST','PROSTATE','LARGE_INTESTINE'), 'pvalues', outfile.prefix = 'colon_EPHB2_stem_cell_sig_cell_line', type = 'cell')
do_RRHO_single(colon2, c('BREAST','PROSTATE','LARGE_INTESTINE'), 'pvalues', outfile.prefix = 'colon_PTK7_stem_cell_sig_cell_line', type = 'cell')



### INSTABILITY RANK FILES CORRECTED PVALUES
setwd('~/Dropbox/Doug/ICNA_BKPT/CCLE_cell_line_ICNA/')
file_list <- list.files()
names <- strtrim(file_list, nchar(file_list) - 28)
dataset_cell_icna = list()
dataset_cell_bkpt = list()
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i], row.names = 1)
  temp <- temp[!(rownames(temp) == ''),]
  tempicna <- p.adjust(temp$ICNA_vs_RNA_p_val, method = 'BH')
  tempbkpt <- p.adjust(temp$bkpt_vs_RNA_p_val, method = 'BH')
  for (j in 1:nrow(temp)) {
    if (temp$ICNA_vs_RNA_rvalue[j] < 0) {
      tempicna[j] <- tempicna[j] * -1
    }
    if (temp$bkpt_vs_RNA_rvalue[j] < 0) {
      tempbkpt[j] <- tempbkpt[j] * -1
    }
  }
  temp1 <- data.frame(ICNA = tempicna, row.names = rownames(temp))
  temp1 <- temp1[order(abs(temp1$ICNA)),,drop = F]
  temp2 <- data.frame(bkpt = tempbkpt, row.names = rownames(temp))
  temp2 <- temp2[order(abs(temp2$bkpt)),,drop = F]
  dataset_cell_icna[[names[i]]] <- temp1
  dataset_cell_bkpt[[names[i]]] <- temp2
}

for (i in 1:length(dataset_cell_icna)) {
  temp <- data.frame(GENE = rownames(dataset_cell_icna[[i]]), ICNA = dataset_cell_icna[[i]]$ICNA)
  write.table(temp, file = paste0(names(dataset_cell_icna)[i], '_pvalue_corrected_cell_line_ICNA.rnk'), quote = F, sep = '\t', row.names = F)
}

for (i in 1:length(dataset_cell_bkpt)) {
  temp <- data.frame(GENE = rownames(dataset_cell_bkpt[[i]]), bkpt = dataset_cell_bkpt[[i]]$bkpt)
  write.table(temp, file = paste0(names(dataset_cell_bkpt)[i], '_pvalue_corrected_cell_line_bkpt.rnk'), quote = F, sep = '\t', row.names = F)
}


setwd('~/Dropbox/Doug/ICNA_BKPT/bkpt.icna files')
file_list <- list.files()
file_list <- file_list[!(file_list == 'others')]
names <- sub('_.*','',file_list)
dataset_ICNA = list()
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i], row.names = 1)
  temp <- temp$ICNA_score
  dataset_ICNA[[names[i]]] <- temp
}

dataset_bkpt = list()
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i], row.names = 1)
  temp <- temp$bkpt_samp
  dataset_bkpt[[names[i]]] <- temp
}

nsets = NULL
for (i in 1:length(dataset_ICNA)) {
  nsets <- c(nsets, length(dataset_ICNA[[i]]))
}
names(nsets) <- names
write.table(nsets, file = 'num_samples_per_tumor.txt', sep = '\t', quote = F, col.names = F)


icnaexpdata <- read.delim('~/Dropbox/Doug/pearson_ICNA_vs_RNA_rvalue_avg.txt')
bkptexpdata <- read.delim('~/Dropbox/Doug/pearson_bkpt_vs_RNA_rvalue_avg.txt')

for (i in 1:(ncol(icnaexpdata)-2)) {
  rvalues <- icnaexpdata[,i+1]
  negs <- rvalues < 0
  tvalues <- rvalues / (sqrt((1-(rvalues)^2)/(length(rvalues)-2)))
  pvalues = 2*pt(-abs(tvalues), df=(nsets[[i]]-1))
  adjpvalues <- p.adjust(pvalues, method = 'BH')
  adjpvalues[negs] <- adjpvalues[negs] * -1
  adjpvalues <- data.frame(GENE = icnaexpdata$Gene, ICNA = adjpvalues)
  adjpvalues <- adjpvalues[order(abs(adjpvalues$ICNA)),,drop = F]
  write.table(adjpvalues, file = paste0(names(nsets)[i], '_pvalue_corrected_tumor_ICNA.rnk'), row.names = F, sep = '\t', quote = F)
}

### Jennifer data
setwd("~/Dropbox/Doug/Jennifer")
library(DESeq2)
library(RRHO)

cell_line <- read.delim("~/Dropbox/Doug/ICNA_BKPT/melanoma_RNA_pearson_correlation.txt")
icna <- data.frame(gene = cell_line$X, ICNA = cell_line$ICNA_vs_RNA_rvalue)
breakpt <- data.frame(gene = cell_line$X, bkpt = cell_line$bkpt_vs_RNA_rvalue)

GIdata_icna <- read.delim(file = '~/Dropbox/Doug/pearson_ICNA_vs_RNA_rvalue_avg.txt')
#GIdata_icna <- data.frame(gene = GIdata_icna$Gene, rvalue = GIdata_icna$SKCM)
GIdata_bkpt <- read.delim(file = '~/Dropbox/Doug/pearson_bkpt_vs_RNA_rvalue_avg.txt')
#GIdata_bkpt <- data.frame(gene = GIdata_bkpt$Gene, rvalue = GIdata_bkpt$SKCM)

source('~/Dropbox/Doug/generate_RRHO_plots_all_cancers.R')
create_RRHO_file_intersect_then_rank(icna,GIdata_icna, x_axis = 'cell_line_ICNA', y_axis = 'tumor_ICNA', reverse = F, ranked = F, outputdir = getwd(), plot = T)
create_RRHO_file_intersect_then_rank(breakpt,GIdata_bkpt, x_axis = 'cell_line_bkpt', y_axis = 'tumor_bkpt', reverse = F, ranked = F, outputdir = getwd())


### Ranked lists for GSEA
dir.create('~/GraeberLab/GSEA')
setwd('~/GraeberLab/GSEA')
for (i in 2:ncol(GIdata_icna)) {
  GIdata <- data.frame(gene = GIdata_icna$Gene, rvalue = GIdata_icna[,i])
  GIdata <- GIdata[order(-GIdata$rvalue),]
  write.table(GIdata, file = paste0(colnames(GIdata_icna)[i], '_tumor_ICNA_rvalues.rnk'), row.names = F, quote = F, sep = '\t')
}

for (i in 2:ncol(GIdata_bkpt)) {
  GIdata <- data.frame(gene = GIdata_bkpt$Gene, rvalue = GIdata_bkpt[,i])
  GIdata <- GIdata[order(-GIdata$rvalue),]
  write.table(GIdata, file = paste0(colnames(GIdata_bkpt)[i], '_tumor_bkpt_rvalues.rnk'), row.names = F, quote = F, sep = '\t')
}

for (i in 2:ncol(GIdata_bkpt)) {
  GIdata <- data.frame(gene = GIdata_bkpt$Gene, blank = GIdata_bkpt[,i])
  colnames(GIdata)[2] <- colnames(GIdata_bkpt)[i]
  create_RRHO_file_intersect_then_rank(pvalues,GIdata, x_axis = x_axis, y_axis = colnames(GIdata_bkpt)[i], outputdir = paste0('./RRHO_bkpt_', outfile.prefix), reverse = reverse, ranked = T)
}

write.table(icna, file = 'melanoma_cell_line_ICNA.rnk', row.names = F, quote = F, sep = '\t')


setwd('/Volumes/data0/users/dyao')
file_list <- list.files('GSEA_cell_line_bkpt')
file_list2 <- list.files('GSEA_cell_line_ICNA')
file_list3 <- list.files('GSEA_tumor_bkpt')
file_list4 <- list.files('GSEA_tumor_ICNA')
file_list5 <- list.files('GSEA_stem_sig')

names <- gsub('_tumor_bkpt_rvalues', '', file_list3)
dir.create('GSEA_top_stem')
setwd('GSEA_top_stem')
for (i in 1:length(file_list5)) {
  gsea_obtain_top_single(paste0('../GSEA_stem_sig/', '/', file_list5[i]), file_list5[i])
}

directory1 = '../GSEA_cell_line_bkpt/BONE_cell_line_bkpt'
label = 'BONE'

#### CELL LINE BKPT VS ICNA
setwd('~/Dropbox/Doug/ICNA_BKPT')
create_RRHO_file_intersect_then_rank(icna, breakpt, x_axis = 'cell_line_ICNA', y_axis = 'cell_line_bkpt', reverse = F, ranked = F, outputdir = getwd())
raw <- read.delim('melanoma_bkpts_sorted_by_bkpts_per_chrom.txt')
plot(raw$ICNA_score, raw$bkpt_samp, main = 'ICNA vs Breakpoints melanoma cell lines', xlab = 'ICNA score', ylab = 'Breakpoints', pch = 19)
pc_data <- read.delim('~/Dropbox/Doug/Jennifer/PCA_all/all_not_scaled_scores.txt', row.names = 1)
rownames(pc_data) <- strtrim(rownames(pc_data), nchar(rownames(pc_data))-5)
rownames(pc_data) <- toupper(rownames(pc_data))
raw$Sample <- toupper(raw$Sample)
intersect <- intersect(rownames(pc_data),raw$Sample)
pc_data <- pc_data[intersect,]
rownames(raw) <- raw$Sample
raw <- raw[intersect,]
plot(pc_data$PC1, raw$ICNA_score, main = 'PC1 score vs ICNA melanoma cell lines', xlab = 'PC1 score', ylab = 'ICNA', pch = 19)
abline(lm(raw$ICNA_score~pc_data$PC1), col="black")
plot(pc_data$PC1, raw$bkpt_samp, main = 'PC1 score vs breakpoints melanoma cell lines', xlab = 'PC1 score', ylab = 'breakpoints', pch = 19)
abline(lm(raw$bkpt_samp~pc_data$PC1), col="black")


### BOXPLOTS RANGE OF ICNA AND BKPTS 
### TUMOR
setwd('~/Dropbox/Doug/ICNA_BKPT/bkpt.icna files')
file_list <- list.files()
file_list <- file_list[!(file_list == 'others')]
names <- sub('_.*','',file_list)

dataset_ICNA = list()
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i])
  temp <- temp$ICNA_score
  dataset_ICNA[[names[i]]] <- temp
}

dataset_bkpt = list()
for (i in 1:length(file_list)) {
  temp <- read.delim(file_list[i])
  temp <- temp$bkpt_samp
  dataset_bkpt[[names[i]]] <- temp
}

medians.ICNA <- sapply(dataset_ICNA, median)
medians.bkpt <- sapply(dataset_bkpt, median)
medians <- cbind(medians.ICNA, medians.bkpt)

png(filename = 'ICNA_vs_bkpts_all_TCGA_tumors.png', height = 750, width = 750)
plot(medians, main = 'ICNA vs number of breakpoints for all TCGA tumors', xlab = 'ICNA score', ylab = 'Number of breakpoints', xlim = c(-10000000, 1e9))
text(medians, labels = rownames(medians), adj = c(1.25,1), cex = 0.8)
dev.off()

sort_list_by_median <- function(list) {
  medians <- sapply(list, median)
  sorted.list <- list[order(medians)]
  return (sorted.list)
}

dataset_ICNA <- sort_list_by_median(dataset_ICNA)
dataset_bkpt <- sort_list_by_median(dataset_bkpt)

png(filename = 'Range_of_ICNA_scores_for_all_TCGA_tumors.png', width = 1000, height = 750)
boxplot(dataset_ICNA, main = 'Range of ICNA scores for all TCGA tumors', xlab = '', ylab = 'ICNA score', las = 3, outline = F)
dev.off()
png(filename = 'Range_of_breakpoints_for_all_TCGA_tumors.png', width = 1000, height = 750)
boxplot(dataset_bkpt, main = 'Range of number of breakpoints for all TCGA tumors', xlab = '', ylab = 'Breakpoints', las = 3, outline = F)
dev.off()


### CELL
setwd('~/Dropbox/Doug/ICNA_BKPT/CCLE_cell_line_ICNA/')
nums <- read.delim('~/Dropbox/Doug/ICNA_BKPT/num_samples_per_cell_line.txt', header = F)
nums <- nums[nums$V2 > 2,]
file_list <- list.files()
file_list <- file_list[!grepl('cell_line', file_list)]
names <- sub('_RNA_pearson_correlation.txt','',file_list)
names <- names[!names == 'cell_line']
names[names == 'CENTRAL_NERVOUS_SYSTEM'] <- 'CNS'
names[names == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"] <- 'BLOOD'
names[names == "UPPER_AERODIGESTIVE_TRACT"] <- 'UPPER_AERO_TRACT'
names <- sub('_',' ',names)
names <- sub('_',' ',names)
names <- paste0(names, ' (n=', nums$V2, ')')
data <- read.delim('~/Dropbox/Doug/ICNA_BKPT/cell_line_bkpts_sorted_by_bkpts_per_chrom.txt')

dataset_ICNA = list()
for (i in 1:length(file_list)) {
  temp <- data[grepl(nums$V1[i], data$Sample),]
  temp <- temp$ICNA_score
  dataset_ICNA[[names[i]]] <- temp
}

dataset_bkpt = list()
for (i in 1:length(file_list)) {
  temp <- data[grepl(nums$V1[i], data$Sample),]
  temp <- temp$bkpt_samp
  dataset_bkpt[[names[i]]] <- temp
}

medians.ICNA <- sapply(dataset_ICNA, median)
medians.bkpt <- sapply(dataset_bkpt, median)
medians <- cbind(medians.ICNA, medians.bkpt)

png(filename = 'ICNA_vs_bkpts_all_CCLE_cancer_cell_lines.png', height = 750, width = 750)
plot(medians, main = 'ICNA vs number of breakpoints for all CCLE cancer cell lines', xlab = 'ICNA score', ylab = 'Number of breakpoints', xlim = c(5e8,1.2e9))
text(medians, labels = rownames(medians), adj = c(NA, 1.25), cex = 0.8)
dev.off()

sort_list_by_median <- function(list) {
  medians <- sapply(list, median)
  sorted.list <- list[order(medians)]
  return (sorted.list)
}

dataset_ICNA <- sort_list_by_median(dataset_ICNA)
dataset_bkpt <- sort_list_by_median(dataset_bkpt)

png(filename = 'Range_of_ICNA_scores_for_all_CCLE_cell_lines.png', width = 1000, height = 750)
par(oma = c(6,0,0,0))
boxplot(dataset_ICNA, main = 'Range of ICNA scores for all CCLE cell lines', xlab = '', ylab = 'ICNA score', las = 3, outline = F, cex = 0.8, cex.axis = 0.8)
dev.off()
png(filename = 'Range_of_breakpoints_for_all_CCLE_cell_lines.png', width = 1000, height = 750)
par(oma = c(6,0,0,0))
boxplot(dataset_bkpt, main = 'Range of number of breakpoints for all CCLE cell lines', xlab = '', ylab = 'Breakpoints', las = 3, outline = F, cex = 0.8, cex.axis = 0.8)
dev.off()



#### PCA JENNIFER
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1, check.names = F) 
data<- round(data)
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
column_data <- data.frame(data = colnames(data))
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~data)
rlog_data <- rlog(dds)

source('pcaplotting.R')
dir.create('./PCA_all')
setwd('./PCA_all')
rlog_data_PCA <- t(assay(rlog_data))
rlog_data_PCA <- rlog_data_PCA[,apply(rlog_data_PCA, 2, var, na.rm=TRUE) != 0]
output_both_PCA_files_and_graphs(dataset = rlog_data_PCA, outfile_prefix = 'all', colors = factor(c(rep(1,ncol(data)))), shapes = factor(c(rep(1,ncol(data)))), pcs = 4)  


####INTRINSIC RESISTANCE
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1) 
data<-round(data)
a <- paste0('M', c(244,296,381,410,233,243,229,397,375,395,230,285,408,417,249), '_DMSO')
data <- data[,a]
data <- data[ rowSums(data) > 1, ]
sample_condition <- c('sensitive', 'sensitive', 'resistant', 'resistant', 'resistant', 'sensitive', 'sensitive', 'sensitive', 'sensitive', 'sensitive', 'sensitive', 'sensitive', 'sensitive', 'resistant', 'sensitive')
column_data <- data.frame(sample_condition)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition')
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~condition)
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c('condition', 'resistant', 'sensitive'), independentFiltering = F, cooksCutoff = F))
res <- res[order(res$pvalue),]
res <- na.omit(res)

values <- log_pvalues_from_deseq(res)
#create_RRHO_file_intersect_then_rank(pvalues,icna, x_axis = 'all_intrinsic_resistant_vs_sensitive', y_axis = 'cell_line_icna', reverse = F, ranked = F, outputdir = getwd())
dir.create('./all_intrinsic_resistant_vs_sensitive_VS_cell_line_bkpt')
create_RRHO_file_intersect_then_rank(pvalues,breakpt, x_axis = 'all_intrinsic_resistant_vs_sensitive', y_axis = 'cell_line_bkpt', reverse = F, ranked = F, outputdir = './all_intrinsic_resistant_vs_sensitive_VS_cell_line_bkpt')



####ACQUIRED RESISTANCE
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1) 
data<-round(data)
data <- data[,c(3,4,7,8,12,13,29,30,31,32,40,41)]
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
sample_condition <- (rep(c('nonAR','AR'),6))
sample_type <- rep(c('M229','M238','M249','M395','M397','M409'), each = 2)
column_data <- data.frame(sample_condition,sample_type)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition','type')
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~type + condition)
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c('condition', 'AR', 'nonAR'), independentFiltering = F, cooksCutoff = F))
res <- res[order(res$pvalue),]
res <- na.omit(res)
pvalues = NULL
for (i in 1:nrow(res)) {
  d <- res[i,]
  pvalue <- -log(d$pvalue)
  if (d$stat < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = rownames(res), pvalue = pvalues)

#create_RRHO_file_intersect_then_rank(pvalues,icna, x_axis = 'all_acquired_resistant_vs_sensitive', y_axis = 'cell_line_icna', reverse = F, ranked = F, outputdir = getwd())

dir.create('./all_acquired_resistant_vs_sensitive_VS_cell_line_bkpt')
create_RRHO_file_intersect_then_rank(pvalues,breakpt, x_axis = 'all_acquired_resistant_vs_sensitive', y_axis = 'cell_line_bkpt', reverse = F, ranked = F, outputdir = './all_acquired_resistant_vs_sensitive_VS_cell_line_bkpt')


####ACQUIRED RESISTANCE W/DIFFERENTIATION
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1) 
data<-round(data)
data <- data[,c(3,4,7,8,12,13,29,30,31,32,40,41)]
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
data <- data[,1:4]
sample_condition <- (rep(c('nonAR','AR'),2))
sample_type <- rep(c('M229','M238'), each = 2)
column_data <- data.frame(sample_condition,sample_type)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition','type')
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~type + condition)
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c('condition', 'AR', 'nonAR'), independentFiltering = F, cooksCutoff = F))
res <- res[order(res$pvalue),]
res <- na.omit(res)
pvalues = NULL
for (i in 1:nrow(res)) {
  d <- res[i,]
  pvalue <- -log(d$pvalue)
  if (d$stat < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = rownames(res), pvalue = pvalues)

#create_RRHO_file_intersect_then_rank(pvalues,icna, x_axis = 'dedifferentiated_acquired_resistant_vs_sensitive', y_axis = 'cell_line_icna', reverse = F, ranked = F, outputdir = getwd())

dir.create('./dedifferentiated_acquired_resistant_vs_sensitive_VS_cell_line_bkpt')
create_RRHO_file_intersect_then_rank(pvalues,breakpt, x_axis = 'dedifferentiated_acquired_resistant_vs_sensitive', y_axis = 'cell_line_bkpt', reverse = F, ranked = F, outputdir = './dedifferentiated_acquired_resistant_vs_sensitive_VS_cell_line_bkpt')

####ACQUIRED RESISTANCE W/O DIFFERENTIATION
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1) 
data<-round(data)
data <- data[,c(3,4,7,8,12,13,29,30,31,32,40,41)]
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
data <- data[,5:10]
sample_condition <- (rep(c('nonAR','AR'),3))
sample_type <- rep(c('M249','M395','M397'), each = 2)
column_data <- data.frame(sample_condition,sample_type)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition','type')
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~type + condition)
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c('condition', 'AR', 'nonAR'), independentFiltering = F, cooksCutoff = F))
res <- res[order(res$pvalue),]
res <- na.omit(res)
pvalues = NULL
for (i in 1:nrow(res)) {
  d <- res[i,]
  pvalue <- -log(d$pvalue)
  if (d$stat < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = rownames(res), pvalue = pvalues)

#create_RRHO_file_intersect_then_rank(pvalues,icna, x_axis = 'not_dedifferentiated_acquired_resistant_vs_sensitive', y_axis = 'cell_line_icna', reverse = F, ranked = F, outputdir = getwd())
dir.create('./not_dedifferentiated_acquired_resistant_vs_sensitive_VS_cell_line_bkpt')
create_RRHO_file_intersect_then_rank(pvalues,breakpt, x_axis = 'not_dedifferentiated_acquired_resistant_vs_sensitive', y_axis = 'cell_line_bkpt', reverse = F, ranked = F, outputdir = './not_dedifferentiated_acquired_resistant_vs_sensitive_VS_cell_line_bkpt')


###DIFFERENTIATION TRAJECTORY
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1) 
data<-round(data)
a <- paste0('M', c(244,296,381,410,233,243,255,409,418,229,376,397,375,395,202,230,285,368,408,417,249), '_DMSO')
data <- data[,a]
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
sample_condition <- c(rep('Undifferentiated',4), rep('Neural_crest_like',5), rep('Transitory',5), rep('Melanocytic',7))
column_data <- data.frame(sample_condition)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition')
conditions <- c('Undifferentiated', 'Neural_crest_like','Transitory','Melanocytic')
conditions <- combn(conditions,2)
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~condition)
dds <- DESeq(dds)

for (x in 1:ncol(conditions)) {
  res <- data.frame(results(dds, contrast = c('condition', conditions[1,x], conditions[2,x]), independentFiltering = F, cooksCutoff = F))
  res <- res[order(res$pvalue),]
  res <- na.omit(res)
  pvalues = NULL
  for (i in 1:nrow(res)) {
    d <- res[i,]
    pvalue <- -log(d$pvalue)
    if (d$stat < 0) {
      pvalue = (pvalue * -1)
    }
    pvalues <- c(pvalues,pvalue)
  }
  pvalues <- data.frame(gene = rownames(res), pvalue = pvalues)

  #create_RRHO_file_intersect_then_rank(pvalues,icna, x_axis = paste0(conditions[1,x],'_vs_', conditions[2,x]), y_axis = 'cell_line_icna', reverse = F, ranked = F, outputdir = getwd())
  dir.create(paste0('./', conditions[1,x],'_vs_', conditions[2,x], '_VS_cell_line_bkpt'))
  create_RRHO_file_intersect_then_rank(pvalues,breakpt, x_axis = paste0(conditions[1,x],'_vs_', conditions[2,x]), y_axis = 'cell_line_bkpt', reverse = F, ranked = F, outputdir = paste0('./', conditions[1,x],'_vs_', conditions[2,x], '_VS_cell_line_bkpt'))
}


### DIFFERENTIATION W/GROUPS
setwd('~/Dropbox/Doug/Jennifer/')
data<-read.csv(file = "Merck_RNASeq_expected_count.genes.results.csv", row.names = 1) 
data<-round(data)
a <- paste0('M', c(244,296,381,410,233,243,255,409,418,229,376,397,375,395,202,230,285,368,408,417,249), '_DMSO')
data <- data[,a]
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
sample_condition <- c(rep('Undifferentiated_Neural_crest_like',9), rep('Transitory_Melanocytic',12))
column_data <- data.frame(sample_condition)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition')
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~condition)
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c('condition', 'Undifferentiated_Neural_crest_like', 'Transitory_Melanocytic'), independentFiltering = F, cooksCutoff = F))
res <- res[order(res$pvalue),]
res <- na.omit(res)
pvalues = NULL
for (i in 1:nrow(res)) {
  d <- res[i,]
  pvalue <- -log(d$pvalue)
  if (d$stat < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = rownames(res), pvalue = pvalues)
write.table(pvalues,file = '/Volumes/data0/users/dyao/undiff_vs_diff_melanoma.rnk', quote = F, row.names = F, sep = '\t')
create_RRHO_file_intersect_then_rank(pvalues,GIdata_icna, x_axis = 'Undifferentiated_Neural_crest_like_vs_Transitory_Melanocytic', y_axis = 'tumor_icna', reverse = F, ranked = F, outputdir = getwd())


#### PCs
pc1 <- read.delim('Melanoma_cell_Line_PC1_loadings.rnk')
pc2 <- read.delim('Melanoma_cell_Line_PC2_loadings.rnk')
pc1_reverse <- data.frame(Gene = pc1$Gene, pc1 = pc1$PC1 * (-1))
write.table(pc1_reverse, file = 'melanoma_PC1.rnk', row.names = F, quote = F, sep = '\t')
create_RRHO_file_output_table(pc1,GIdata_icna,'PC1','TUMOR_SKCM_ICNA', rank1 = 'asc')
create_RRHO_file_output_table(pc1,icna,'PC1','CELL_SKCM_ICNA', rank1 = 'asc')

a <- create_RRHO_file_intersect_then_rank(pc1,icna, x_axis = 'PC1', y_axis = 'tumor_icna', reverse = T, ranked = F, outputdir = getwd(), plot = F)
create_RRHO_file_intersect_then_rank(pc1,icna, x_axis = 'PC1', y_axis = 'cell_line_icna', reverse = T, ranked = F, outputdir = getwd())
create_RRHO_file_intersect_then_rank(pc2,GIdata_icna, x_axis = 'PC2', y_axis = 'tumor_icna', reverse = T, ranked = F, outputdir = getwd())
create_RRHO_file_intersect_then_rank(pc2,icna, x_axis = 'PC2', y_axis = 'cell_line_icna', reverse = T, ranked = F, outputdir = getwd())


#### STEM CELL SIGNATURES
setwd('~/Dropbox/Doug/ICNA_BKPT')
breast <- read.delim('breast_illumhumanWG-6_v3_MaSC-Luminal_maxavg_Pval.tsv', check.names = F)
breast <- data.frame(gene = rownames(breast), pvalue = breast$logP.Value)
prostate <- read.delim('prostate_Quant_CD49fHi.Benign-CD49fLo.Benign_maxavg_Pval.tsv', check.names = F)
prostate <- data.frame(gene = rownames(prostate), pvalue = prostate$logP.Value)
colon1 <- read.delim('colon_RMA_EPHB2HighMed-EPHB2LowNeg_maxavg_Pval.tsv', check.names = F)
colon1 <- data.frame(gene = rownames(colon1), pvalue = colon1$logP.Value)
colon2 <- read.delim('colon_RMA_PTK7HighMed-PTK7LowNeg_maxavg_Pval.tsv', check.names = F)
colon2 <- data.frame(gene = rownames(colon2), pvalue = colon2$logP.Value)

write.table(breast, file = 'breast_stem_sig.rnk', row.names = F, quote = F, sep = '\t')
write.table(prostate, file = 'prostate_stem_sig.rnk', row.names = F, quote = F, sep = '\t')
write.table(colon1, file = 'colon_EPHB2_stem_sig.rnk', row.names = F, quote = F, sep = '\t')
write.table(colon2, file = 'colon_PTK7_stem_sig.rnk', row.names = F, quote = F, sep = '\t')


do_RRHO_single(breast, c('BRCA','PRAD','COAD','READ'), 'pvalues', outfile.prefix = 'breast_stem_cell_sig')
do_RRHO_single(prostate, c('BRCA','PRAD','COAD','READ'), 'pvalues', outfile.prefix = 'prostate_stem_cell_sig')
do_RRHO_single(colon1, c('BRCA','PRAD','COAD','READ'), 'pvalues', outfile.prefix = 'colon_EPHB2_COAD_stem_cell_sig')
do_RRHO_single(colon2, c('BRCA','PRAD','COAD','READ'), 'pvalues', outfile.prefix = 'colon_PTK7_COAD_stem_cell_sig')

icna <- read.delim('C:/Users/DWYao/Dropbox/Doug/ICNA_BKPT/tumor_ICNA_rank_files/COADREAD_tumor_ICNA_rvalues.rnk')
bkpt <- read.delim('C:/Users/DWYao/Dropbox/Doug/ICNA_BKPT/tumor_bkpt_ranks/COADREAD_tumor_bkpt_rvalues.rnk')
create_RRHO_file_intersect_then_rank(breast, icna, x_axis = 'breast_pvalues', y_axis = 'COADREAD_ICNA_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(breast, bkpt, x_axis = 'breast_pvalues', y_axis = 'COADREAD_bkpt_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(prostate, icna, x_axis = 'prostate_pvalues', y_axis = 'COADREAD_ICNA_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(prostate, bkpt, x_axis = 'prostate_pvalues', y_axis = 'COADREAD_bkpt_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(colon1, icna, x_axis = 'colon_EPHB2_pvalues', y_axis = 'COADREAD_ICNA_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(colon1, bkpt, x_axis = 'colon_EPHB2_pvalues', y_axis = 'COADREAD_bkpt_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(colon2, icna, x_axis = 'colon_PTK7_pvalues', y_axis = 'COADREAD_ICNA_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(colon2, bkpt, x_axis = 'colon_PTK7_pvalues', y_axis = 'COADREAD_bkpt_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(ave2, icna, x_axis = 'ave_pvalues', y_axis = 'COADREAD_ICNA_score', outputdir = getwd(), reverse = F, ranked = F)
create_RRHO_file_intersect_then_rank(ave2, bkpt, x_axis = 'ave_pvalues', y_axis = 'COADREAD_bkpt_score', outputdir = getwd(), reverse = F, ranked = F)


#### Mellinghoff data
setwd('~/Dropbox/Doug/Mellinghoff')
data<-read.table(file = '~/Dropbox/Doug/Mellinghoff/Mellinghoff_matrices/Mellinghoff_raw_gene_count_matrix.txt',row.names=1,header=T,sep="\t",stringsAsFactors = F) 
data<-round(data)
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]

cell<-factor(as.character(sapply(colnames(data), function(x) ifelse(grepl("AM38",x), "AM38", "NMCG1"))))
drug<-factor(as.character(sapply(colnames(data), function(x) if(grepl(".Vem.",x)) "Vem" else if(grepl(".Tram.",x)) "Tram" else if(grepl(".LY.",x)) "LY")),levels=c("Vem","Tram","LY") )
time<-factor(as.character(sapply(colnames(data), function(x) if(grepl("\\.48_",x)) "48" else if(grepl("\\.24_",x)) "24" else if(grepl("\\.12_",x)) "12" else if(grepl("\\.6_",x)) "6"  else if(grepl("\\.1_",x)) "1" else if(grepl("\\.0_",x)) "0")),levels=c(0,1,6,12,24,48))
patient<-factor(as.character(sapply(colnames(data), function(x) if(grepl("PITT_0052",x)) "PITT_0052" else if(grepl("KIM_0453",x)) "KIM_0453" else if(grepl("KIM_0452",x)) "KIM_0452")),levels=c("PITT_0052","KIM_0453","KIM_0452"))
pooledtime<-factor(sapply(time,function(x) if(x==24 || x == 48) 'Late' else if (x==12) 'Middle' else if (x == 6 || x == 1) 'Early' else if (x == 0) 'DMSO'),levels=c('Late', 'Middle', 'Early','DMSO'))
cell_drug_time<-paste0(cell,drug,time)
drug_time<-paste0(drug,time)
cell_drug<-paste0(cell,drug)
drug_pooledtime <- paste0(drug,pooledtime)


column_data <- data.frame(cell,drug,time,cell_drug_time,drug_time,cell_drug)
rownames(column_data)<-colnames(data)

dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~cell + drug_time)
dds <- DESeq(dds) 
res <- results(dds, contrast = list(c("drug_timeVem24","drug_timeVem48"), c("drug_timeVem0")), listValues = c(0.5,-1), independentFiltering = F, cooksCutoff = F)
res <- na.omit(res)
pvalues <- log_pvalues_from_deseq(res)

res_LY <- results(dds, contrast = list(c("drug_timeLY24","drug_timeLY48"), c("drug_timeLY0")), listValues = c(0.5,-1), independentFiltering = F, cooksCutoff = F)
res_LY <- na.omit(res_LY)
pvalues_LY <- log_pvalues_from_deseq(res_LY)
do_RRHO(pvalues_LY, outfile.prefix = 'LY_Late_vs_DMSO')


#### Sanaz data
setwd('~/Dropbox/Doug/Sanaz')
source('collapse_data.R')
data <- read.delim('GSE70072_HumanSerousCancer_rawCounts.txt', row.names = 1, stringsAsFactors = F)
rownames(data) <- NULL
instabgenes <- read.delim('~/Dropbox/Doug/ICNA_BKPT/RANK_FILES/tumor_bkpt_ranks/ACC_tumor_bkpt_rvalues.rnk')
instabgenes <- instabgenes$gene
data <- data[data$gname %in% instabgenes,]
data <- data.frame(collapse_data(data, s.start = 2, type = 'maxavg', group = 'gname'))
data <- data[ rowSums(data[1:ncol(data)]) > 1, ]
sample_condition <- (rep(c('CA125negative','CA125positive'),10))
sample_type <- rep(c(LETTERS[1:10]), each = 2)
column_data <- data.frame(sample_condition,sample_type)
rownames(column_data) <- colnames(data)
colnames(column_data) <- c('condition','type')
dds <- DESeqDataSetFromMatrix(countData = data, colData = column_data, design = ~type + condition)
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c('condition','CA125negative','CA125positive'), cooksCutoff = F, independentFiltering = F))
res <- res[order(res$pvalue),]
res <- na.omit(res)
pvalues = NULL
for (i in 1:nrow(res)) {
  d <- res[i,]
  pvalue <- -log(d$pvalue)
  if (d$stat < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = rownames(res), pvalue = pvalues)
do_RRHO_single(pvalues, cancer = 'OV', x_axis = 'CA125neg_vs_CA125pos', outfile.prefix = 'CA125')
pvalues <- pvalues[order(pvalues$pvalue),]
write.table(pvalues, file = 'ca125negvpos_short.rnk', quote = F, sep = '\t', row.names = F)

# OTHER DATA
setwd('~/GraeberLab/Scripts')
data <- read.delim('RMA_resistant-control_maxavg_Pval.tsv')
data <- data.frame(gene = rownames(data), pvalues = data$logP.Value)
do_RRHO(data, outfile.prefix = 'GSE24460_resistant_vs_control')

data2 <- read.delim('RMA_G1-G0_maxavg_Pval.tsv')
data2 <- data.frame(gene = rownames(data2), pvalues = data2$logP.Value)
do_RRHO(data2, outfile.prefix = 'GSE51373_resistant_vs_control')

setwd('~/Dropbox/Doug/TGFB')
data <- read.delim('~/Dropbox/Doug/TGFB/mcf10a_tgf_ttest.txt')
data <- na.omit(data)
pvalues = NULL
for (i in 1:nrow(data)) {
  d <- data[i,]
  pvalue <- -log(d$p.value)
  if (d$statistic < 0) {
    pvalue = (pvalue * -1)
  }
  pvalues <- c(pvalues,pvalue)
}
pvalues <- data.frame(gene = data$X, pvalue = pvalues)
do_RRHO_single(pvalues, 'BRCA', 'pvalues', outfile.prefix = 'TGF')




