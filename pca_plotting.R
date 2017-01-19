#### pca_plotting.R
#### Functions for performing principal component analysis (PCA)


library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(GGally)
library(Hmisc)
library(limma)

#Input dataset should have nothign but actual #'s.(you can have header and rownames of course)
#INPUT: |dataset| is a dataset where columns are genes, and rows are samples, because we are using the pr.comp function.
#|outfile_prefix| is the output name without the extension for example "this_is_the_dataset_name"
#|colors| is a list of the colors, its length is equal to the # of samples e.g. factor(rep(c(1,2),each=12)) .. 
#factor has 2 levels i.e. first 12 samples one color 2nd 12 samples another color
#|shapes| is a list of the shapes, its length is equal to # of samples e.g. factor(rep(c(1,2,3,4,5,6,7,8,9,10,11,12) , each =2)) 
#|pcs| is the number of pcs to show can vary default value is 4, there is a page for each pc comparison so if pcs=7 there are 7 choose 2 = 21 pages. Look at the var explained
#plot to see how high you should go. 4 is a good place to start.
#factor has 12 levels, first two samples one shape, 2nd two next shape...etc
#OUTPUT: output both the scaled/centered and unscaled/centered  scores,loadings, and varexplained files, as well ouputs
#A pdf that has by default has the 4 principal component comparisons, first page is scree plot, 2nd page is grid view, with a density plot as well, 
#each subsequent page is a comparison of the PCs, Can take up to 9 colors, basically unlimited shapes
output_both_PCA_files_and_graphs<-function(dataset,outfile_prefix,colors,shapes,pcs=4){
  pr_data_scaled <- prcomp(dataset,scale=T,center=T) 
  pr_data_scaled_scores <- pr_data_scaled$x
  pr_data_scaled_var <- summary(pr_data_scaled)[[6]][2,]
  scores_file1 <- paste0(outfile_prefix, "_scaled_scores.txt")
  variance_file1 <- paste0(outfile_prefix, "_scaled_varexplained.txt")
  
  write.table(pr_data_scaled$rotation, file = paste0(outfile_prefix, "_scaled_loadings.txt"), col.names=NA, row.names=T, sep="\t", quote=F)
  write.table(pr_data_scaled$x, file = paste0(outfile_prefix, "_scaled_scores.txt"), col.names=NA, row.names=T, sep="\t", quote=F)
  write.table(summary(pr_data_scaled)[[6]][2,], file=paste0(outfile_prefix, "_scaled_varexplained.txt"), sep="\t", quote=F, col.names=F, row.names=T)
  
  pr_data_not_scaled <- prcomp(dataset,scale=F,center=T) 
  pr_data_not_scaled_scores <- pr_data_not_scaled$x
  pr_data_not_scaled_var <- summary(pr_data_not_scaled)[[6]][2,]
  scores_file2 <- paste0(outfile_prefix, "_not_scaled_scores.txt")
  variance_file2 <- paste0(outfile_prefix, "_not_scaled_varexplained.txt")
  
  write.table(pr_data_not_scaled$rotation, file = paste0(outfile_prefix, "_not_scaled_loadings.txt"), col.names=NA, row.names=T, sep="\t", quote=F)
  write.table(pr_data_not_scaled$x, file = paste0(outfile_prefix, "_not_scaled_scores.txt"), col.names=NA, row.names=T, sep="\t", quote=F)
  write.table(summary(pr_data_not_scaled)[[6]][2,], file=paste0(outfile_prefix, "_not_scaled_varexplained.txt"), sep="\t", quote=F, col.names=F, row.names=T)
  
  #scaled graph
  draw_qplot_PCA(scores_file1,variance_file1,colors=colors,shapes=shapes,pcs=pcs,data=t(dataset),scale_flag=T)
  #unscaled graph
  draw_qplot_PCA(scores_file2,variance_file2,colors=colors,shapes=shapes,pcs=pcs,data=t(dataset),scale_flag=F)
  
}


draw_ggpairs_PCA <-function(dataset,gg_title,pc=4) {
  print(ggpairs(dataset, columns = c(1:pc),
                lower = list(continuous="points"),
                diag = list(continuous ="blank"),
                upper = list(continuous = "density"),
                title=gg_title))
}



#INPUT: A |scores_file| and |variance file| which are the pca scores file after write.table of the scores output from pr.comp, same for var explained, 
#|pcs| input can be from 2-5, outputs comparisons of all those pcs, by default it is 4. |colors| is a list of the colors, its length is equal to 
# of samples e.g. factor(rep(c(1,2),each=12)) .. factor has 2 levels i.e. first 12 samples one color 2nd 12 samples another color
#|shapes| is a list of the shapes, its length is equal to # of samples e.g. factor(rep(c(1,2,3,4,5,6,7,8,9,10,11,12) , each =2)) factor has 12 levels, 
#first two samples one shape, 2nd two next shape...etc
#|pcs| is the number of pcs to show can vary 
#OUTPUT: A pdf that has by default has the 4 principal component comparisons, first page is scree plot, 2nd page is grid view, with a density plot as well, 
#each subsequent page is a comparison of the PCs, Can take up to 9 colors, basically unlimited shapes


draw_qplot_PCA <-function(scores_file,variance_file,colors,shapes,pcs=4,data=NULL,scale_flag=T){
  dat <- read.table(scores_file, header=TRUE, row.names= 1, check.names=FALSE, sep="\t",stringsAsFactors=FALSE )
  var_dat <- read.table(variance_file,header=FALSE,row.names=1,check.names=FALSE,sep="\t",stringsAsFactors=FALSE)
  PC_pairs= list()
  
  combos<-combn(seq(1,pcs),2)
  PC_pairs<- split(combos,col(combos))
  
  pdf(paste0(removeExt(scores_file),"_PC_PLOTS", ".pdf"))
  #par(mar=c(5,4,4,5))
  plot(var_dat[1:10,], type='b', pch=20, col='black', ylab= 'Variance Percentage', xlab='PC', main=removeExt(scores_file),xaxp=c(1,10,9))
  draw_ggpairs_PCA(dat,paste0(removeExt(scores_file),"_PCs_1_to_",pcs),pc=pcs) 
  if(!is.null(data)) { 
    
    if(scale_flag==T) {
    centered.data<-t(scale(t(data),center = T,scale=T))
    quantile.range <- quantile(centered.data, probs = seq(0, 1, 0.01))
    bks<-seq(quantile.range["1%"],quantile.range["99%"],length.out =20)
    cols <- colorRampPalette(c("blue","white","red"))(length(bks)+1)    
    
 #   do.call("pheatmap",list(centered.data,color = cols, kmeans_k = 20, breaks=bks,border_color="black",name = paste0(gsub("_scores","",removeExt(scores_file)), "_complete_pearson"),clustering_method="complete",clustering_distance_rows="correlation",clustering_distance_cols="correlation",cluster_method = "complete",cluster_rows=T,cluster_cols=T,show_rows = F,show_cols=F, dpi=300))
#    do.call("pheatmap",list(centered.data,color = cols, kmeans_k = 20, breaks=bks,border_color="black",name = paste0(gsub("_scores","",removeExt(scores_file)), "_centroid_pearson"),clustering_method="centroid",clustering_distance_rows="correlation",clustering_distance_cols="correlation",cluster_method = "centroid",cluster_rows=T,cluster_cols=T,show_rows = F,show_cols=F, dpi=300)) 
  }
  if(scale_flag==F){
    centered.data<-t(scale(t(data),center = T,scale=F))
  quantile.range <- quantile(centered.data, probs = seq(0, 1, 0.01))
  bks<-seq(quantile.range["1%"],quantile.range["99%"],length.out =20)
  cols <- colorRampPalette(c("blue","white","red"))(length(bks)+1)    
  
 # do.call("pheatmap",list(centered.data,color = cols, breaks=bks,border_color="black",name = paste0(gsub("_scores","",removeExt(scores_file)), "_complete_pearson"),clustering_method="complete",clustering_distance_rows="correlation",clustering_distance_cols="correlation",cluster_method = "complete",cluster_rows=T,cluster_cols=T,show_rows = F,show_cols=F, dpi=300))
#  do.call("pheatmap",list(centered.data,color = cols, breaks=bks,border_color="black",name = paste0(gsub("_scores","",removeExt(scores_file)), "_centroid_pearson"),clustering_method="centroid",clustering_distance_rows="correlation",clustering_distance_cols="correlation",cluster_method = "centroid",cluster_rows=T,cluster_cols=T,show_rows = F,show_cols=F, dpi=300)) 
  }
  }
  
  
  for (PC in PC_pairs) {     
    qp<- qplot(data=dat, x=dat[, PC[[1]] ],y=dat[, PC[[2]] ], color=colors,shape=shapes)+
      scale_colour_brewer(palette = "Set1")+  
      geom_text(aes(label=rownames(dat)), size =2.5, vjust =-0.4,hjust=-0.04)+
      geom_point(size=2.3)+
      theme_bw()+
      scale_shape_manual(values=1:nlevels(shapes))+
      theme(legend.position="none",plot.title=element_text(size=13,hjust=1))+
      scale_x_continuous(expand = c(0.25,0))+
      labs(x= paste0("PC",PC[[1]]),
           y= paste0("PC",PC[[2]]),
           title= paste0(removeExt(scores_file), "_", paste0("PC",PC[[1]]), "_vs_",paste0("PC",PC[[2]])))
    print(qp)      
  }
  dev.off()
}



###################################################
##########filtering genes by their coefficient of variation
#####################################################



#returns coefficient of variation value for a gene, i.e. row of a dataset, na.rm if you want to remove NAs before calculation

cv_val_gene<-function(vector,na.rm){
  sdx <- sd(vector, na.rm = na.rm)
  val <- sdx/abs(mean(vector, na.rm = na.rm))
  return(val)
}

#returns a coefficient of variation vector of length equal to # rows, each value is the cv value of that row.
cv_val_dataset<-function(dataset,na.rm) {
  cv_vector <- numeric(length=nrow(dataset))
  for(row in 1:nrow(dataset)) {
    cv_vector[row] <- cv_val_gene(dataset[row,],na.rm)
  }
  return(cv_vector)
}

#returns a dataset filtered by a CV cutoff. 
#dataset is the dataset: genes rows, samples columns
#filter argument tells how to perform filtering, 
#"rm_below" which removes those genes with cvs below cutoff,this is the default
#"rm_above" which removes those genes with cvs above cutoff,
#"bounded" which removes those genes not within cutoff bounds
#default na.rm=TRUE, removes NAs before calculation of CVs
#cutoff is the cutoff value, if "bounded" is used the this is a vector of form c(lowvalue,highvalue), e.g.(0.3,0.7)
#which keeps genes with CV values within these values 

coeff_of_v_filter <- function(dataset,filter="rm_below",cutoff,na.rm=T) {
  cv_vector<-cv_val_dataset(dataset,na.rm)
  cv_logical<- logical(length(cv_vector))
  switch(filter,
         rm_below={rows_filtered<- vapply(cv_vector,function(x) {if (x<=cutoff){FALSE} else{TRUE}},logical(1))},
         rm_above={rows_filtered<- vapply(cv_vector,function(x) {if (x>=cutoff){FALSE} else{TRUE}},logical(1))},
         bounded={rows_filtered<- vapply(cv_vector,function(x) {if (x>=cutoff[[2]] || x<=cutoff[[1]]){FALSE} else{TRUE}},logical(1))},
         stop("Invalid filter name")
  )
  return(dataset[rows_filtered,])
}



# quantile.range <- quantile(data.p, probs = seq(0, 1, 0.01))
# bks<-seq(quantile.range2["2%"],quantile.range2["98%"],0.2)
# cols <- colorRampPalette(c("blue","white","red"))(length(bks)+1)
# params <- list(data.p,name=paste0(gsub("_scores","",removeExt(scores_file)),color = cols,breaks=bks,border_color="black",cellwidth=3,cellheight=3,fontsize=3)
# 
# 
# do.call("pheatmap",c(params,clustering_method="complete",filename=paste0(savename,"_completelink_pearson_heatmap.png"),clustering_distance_rows="correlation",clustering_distance_cols="correlation",cluster_method = "complete",show_rows = F,show_cols=F, dpi=300))
# do.call("pheatmap",c(params,clustering_method="centroid",filename=paste0(savename,"_centroidlink_pearson_heatmap.png"),clustering_distance_rows="correlation",clustering_distance_cols="correlation",cluster_method = "centroid",show_rows = F,show_cols=F, dpi=300))
# 



do_heatmap <- function (dataset,name,show_rows,show_cols,cluster_metric,cluster_method,cluster_r,cluster_c ) {
  col.pal <- brewer.pal(9,"Blues")
  hm.parameters <- list(dataset, 
                        color = col.pal,
                        cellwidth = NA, cellheight = NA, scale = "none",
                        treeheight_row = 0,
                        kmeans_k = NA,
                        show_rownames = show_rows, show_colnames = show_cols,
                        main = name,
                        fontsize_row= 5,
                        fontsize_column=5,
                        clustering_method = cluster_method,
                        cluster_rows = cluster_r, cluster_cols = cluster_c,
                        clustering_distance_rows = cluster_metric, 
                        clustering_distance_cols = cluster_metric)
  
  do.call("pheatmap", hm.parameters)
}







