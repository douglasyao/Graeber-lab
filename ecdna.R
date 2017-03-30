### ecdna.R
# Perform data processing for specific ECDNA data sets



# data <- read.delim('C:/Users/DWYao/Dropbox/Doug/ECDNA/GBM_just_cancer.seg', stringsAsFactors = F)
# borders <- c(55084725, 55277031)
# data <- data[data$Chromosome == 7,]
# data <- data[data$Start<borders[2],]
# data <- data[data$End>borders[1],]
# data <- data[order(-data$Segment_Mean),]
# data$Length <- data$End - data$Start


# TCGA GBM
library(DESeq2)
library(org.Hs.eg.db)
source('~/Dropbox/Doug/Sanaz/collapse_data.R')

data <- read.delim('~/Dropbox/Doug/ECDNA/trackNames.tab', stringsAsFactors = F)
data <- data[!duplicated(data$Sample),]
data <- data$Name
cutoff = 'TCGA-12-0776-01A-01D-0333-01'
amp <- data[1:match(cutoff,data)]
notamp <- data[(match(cutoff,data)+1):length(data)]
amp <- strtrim(amp, 15)
amp <- gsub('-','.',amp)
notamp <- strtrim(notamp, 15)
notamp <- gsub('-','.',notamp)


tcga <- read.table('/Volumes/data0-1/users/dyao/TCGA_genes_raw/GBM_expected_count_gene.txt', stringsAsFactors = F, header = 1)
tcga[,2:ncol(tcga)] <- round(exp(tcga[,2:ncol(tcga)])-1)


amp <- intersect(amp,colnames(tcga))
notamp <- intersect(notamp,colnames(tcga))

s <- AnnotationDbi::select(org.Hs.eg.db, keys= keys(org.Hs.eg.db), columns=c("SYMBOL","ENSEMBL"))
m <- match(strtrim(tcga$sample, 15), s$ENSEMBL)
genes <- s$SYMBOL[m]

tcga$sample <- genes
tcga <- na.omit(tcga)
tcga <- collapse_data(tcga, group = 'sample')

tcga <- tcga[,c(amp,notamp)]
tcga <- tcga[rowSums(tcga) > 0,]

column.data <- data.frame(c(rep('amp',length(amp)),rep('notamp',length(notamp))))
rownames(column.data) <- c(amp,notamp)
colnames(column.data) <- 'condition'

dds <- DESeqDataSetFromMatrix(countData = tcga, colData = column.data, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition','amp','notamp'))
res <- data.frame(res[order(res$pvalue),])


### CCLE
source('~/Dropbox/Doug/Sanaz/collapse_data.R')
source('~/GitHub/Graeber-lab/perform_ttests.R')
library(stringr)
ccle <- read.table('~/Documents/CCLE_Expression_Entrez_2012-09-29.gct', header=TRUE, row.names=NULL, sep='\t', check.names = FALSE,stringsAsFactors = FALSE)
ccle <- ccle[,-1]
ccle <- data.frame(collapse_data(ccle, type = 'maxavg', group = 'Description'))
ccle <- ccle[!(rownames(ccle) == ''),]
name_only <- str_match(colnames(ccle), '^[^_]*')
type <- str_match(colnames(ccle), '^[^_]*_(.*$)')[,2]
colnames(ccle) <- name_only
ccle <- ccle[,intersect(colnames(ccle), c(ecpos, ecneg))]
lookup <- data.frame(name = name_only, type = type, stringsAsFactors = F)
lookup <- lookup[lookup$name %in% c(ecpos, ecneg),]
lookup$ec <- sapply(lookup$name, function(x) if (x %in% ecpos) 'yes' else 'no')
rownames(lookup) <- lookup$name
lookup$name <- NULL

annot <- read.csv('~/Dropbox/Doug/ECDNA/12 10 16 Source Data Table 2.csv', stringsAsFactors = F)
annot$Sample.name <- gsub('-', '', toupper(annot$Sample.name))
annot <- annot[annot$Sample.name %in% intersect(annot$Sample.name, colnames(ccle)),]
ecpos <- annot$Sample.name[which(annot$EC.positive == 'yes')]
ecneg <- annot$Sample.name[which(annot$EC.positive == 'no')]


model <- model.matrix(~type + ec, lookup)

log_pvalues_from_ttest(ccle, ecpos, ecneg, output = '~/Dropbox/Doug/ECDNA/CCLE_ECpos_vs_ECneg.rnk')
log_pvalues_from_glm(ccle, model, 'modelecyes', output = '~/Dropbox/Doug/ECDNA/CCLE_ECpos_vs_ECneg_glm.rnk')

write.table(lookup[order(lookup$type),], '~/Dropbox/Doug/ECDNA/types.txt', sep = '\t', col.names = NA, quote = F)
