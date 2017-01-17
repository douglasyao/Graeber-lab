#### RRHO.R
# Performs Rank-Rank Hypergeometric Overlap analysis (Plaisier et al., 2010) from command line
# Requires two input metric files with two columns each: left column containing gene names and right column containing metric
# Other arguments specify title and labels for output plots

library('optparse')
source('ExpressionAnalysis_3.0_v2.R')
dyn.load('libgsl.so')

optionlist = list(make_option(c('--metric_file1','-a'), type = 'character', help = 'metric file 1', metavar = 'character'), 
                  make_option(c('--metric_file2','-b'), type = 'character', help = 'metric file 2', metavar = 'character'), 
                  make_option(c('--x_axis','-x'), type = 'character', help = 'label for x axis', metavar = 'character'),
                  make_option(c('--y_axis','-y'), type = 'character', help = 'label for y axis', metavar = 'character'),
                  make_option(c('--x_sub_left','-c'), type = 'character', default = '', help = 'label for left subtitle on x axis', metavar = 'character'),
                  make_option(c('--x_sub_right','-d'), type = 'character', default = '', help = 'label for right subtitle on x axis', metavar = 'character'),
                  make_option(c('--y_sub_left','-e'), type = 'character', default = '', help = 'label for left subtitle on y axis', metavar = 'character'),
                  make_option(c('--y_sub_right','-f'), type = 'character', default = '', help = 'label for right subtitle on y axis', metavar = 'character'),
                  make_option(c('--outputdir','-o'), type = 'character', default = getwd(), help = 'output directory', metavar = 'character'),
                  make_option(c('--stepsize','-s'), type = 'integer', help = 'stepsize', metavar = 'integer'))

opt_parser = OptionParser(option_list = optionlist)
opt = parse_args(opt_parser)

if (is.null(opt$metric_file1) | is.null(opt$metric_file2)){
  print_help(opt_parser)
  stop('Must supply two input metric files')
}

if (is.null(opt$x_axis)) {
  opt$x_axis <- opt$metric_file1
}

if (is.null(opt$y_axis)) {
  opt$y_axis <- opt$metric_file2
}

metric_file1 = opt$metric_file1
metric_file2 = opt$metric_file2
x_axis = opt$x_axis
y_axis = opt$y_axis
xsubs = c(opt$x_sub_left,opt$x_sub_right)
ysubs = c(opt$y_sub_left,opt$y_sub_right)
outputdir= opt$outputdir
stepsize = opt$stepsize

metric1 <- read.delim(metric_file1)
metric2 <- read.delim(metric_file2)
metric1 <- metric1[!duplicated(metric1[,1]),]
metric2 <- metric2[!duplicated(metric2[,1]),]

metric1[,2] <- 1:nrow(metric1)
metric2[,2] <- 1:nrow(metric2)
rownames(metric1)<-make.names(metric1[,1])
rownames(metric2)<-make.names(metric2[,1])
metric1 <- metric1[,-1,drop=F]
metric2 <- metric2[,-1,drop=F]

common_names<-intersect(rownames(metric1), rownames(metric2))
print(paste0("File1 genes=",length(rownames(metric1))," File2 genes=",length(rownames(metric2))," Intersect genes=",length(common_names)))

metric1 <- metric1[common_names,,drop = F]
metric2 <- metric2[common_names,,drop = F]

if (is.null(stepsize)) RRHO(data.frame(gene = common_names, value = -metric1), data.frame(gene = common_names, value = -metric2), xsubs = xsubs, ysubs = ysubs, alternative = 'two.sided', plots = T, labels = c(x_axis, y_axis), outputdir = outputdir, log10.ind = T) else 
RRHO(data.frame(gene = common_names, value = -metric1), stepsize = stepsize, data.frame(gene = common_names, value = -metric2), xsubs = xsubs, ysubs = ysubs, alternative = 'two.sided', plots = T, labels = c(x_axis, y_axis), outputdir = outputdir, log10.ind = T)





