source('../../exac_constants.R')

load_R_libraries( 'vioplot' )

if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}

gtex_path = 'data/constraint/gtex/'
tissues <- colnames(read.table(sprintf('%sGTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm_medians.txt.gz',gtex_path),header=T))[-1]

## Function for plotting violin plots
myvioplot <- function(list1,list2,list3,color1,color2,color3,gap){
  
  nviolins <- 3
  maxval <- max(c(max(unlist(list1)),max(unlist(list2)),max(unlist(list3))))
  minval <- min(c(min(unlist(list1)),min(unlist(list2)),min(unlist(list3))))
  
  positions <- c(1:nviolins,(nviolins+gap):(nviolins+gap+nviolins-1),(nviolins+nviolins+gap+gap-1):(nviolins+nviolins+gap+gap+nviolins-2))
  
  par(bty='n')
  plot(positions,c(minval-0.1*minval,rep(maxval+0.1*maxval,length(positions)-1)),col='white',axes=F,xlab='',ylab='',xlim=c(0.5,max(positions+0.5)))
  
  h <- vioplot(list1[[1]],list1[[2]],list1[[3]],at=positions[1:nviolins],add=T,col=alpha(color1, 0.4),border=color1,lwd=2)
  points(positions[1:nviolins],h$median[1:3],bg='white',col='black',pch=21,cex=2)
  h <- vioplot(list2[[1]],list2[[2]],list2[[3]],at=positions[(nviolins+1):(nviolins+nviolins)],add=T,col=alpha(color2, 0.4),border=color2,lwd=2)
  points(positions[(nviolins+1):(nviolins+nviolins)],h$median[1:3],bg='white',col='black',pch=21,cex=2)
  h <- vioplot(list3[[1]],list3[[2]],list3[[3]],at=positions[(nviolins+1+nviolins):(nviolins+nviolins+nviolins)],add=T,col=alpha(color3, 0.4),border=color3,lwd=2)
  points(positions[(nviolins+1+nviolins):(nviolins+nviolins+nviolins)],h$median[1:3],bg='white',col='black',pch=21,cex=2)
  
  shift <- 0
  mtext(c('synonymous','missense','protein-truncating'),side=3,line=-0.9,at=positions[seq(ceiling(nviolins/2),nviolins*3,nviolins)]+shift,col=c(color1,color2,color3), cex=1.2)
  mtext(c('low','medium','high'),side=1,line=-0.25,at=positions,cex=1)
  
}

prepare_gtex_data = function() {
  # Load expression data
  transcripts <- read.table(sprintf('%scanonical_transcripts_v75.txt.gz',gtex_path),header=F)
  colnames(transcripts) <- c('GeneID','TranscriptID')
  constraint_gene <- merge(transcripts,constraint,by.x=2,by.y="feature")
  expr <- read.table(sprintf('%sGTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm_medians.txt.gz',gtex_path),header=T)
  
  # Combine brain expression excluding two duplicates: "Brain_Cortex" and "Brain_Cerebellum"
  brain_tissues <- c("Brain_Amygdala","Brain_Anterior_cingulate_cortex","Brain_Caudate","Brain_Cerebellar_Hemisphere","Brain_Frontal_Cortex","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens","Brain_Putamen","Brain_Spinal_cord","Brain_Substantia_nigra")
  expr$Brain <- apply(expr[,which(!is.na(match(colnames(expr),brain_tissues)))],1,median)
  
  # List of all tissues
  tissues <- colnames(expr)[-1]
  
  # Combine expression and constraint data
  gtex <- merge(constraint_gene,expr,by.x="GeneID",by.y=1)
  return(list(gtex, tissues))
}

metrics = c("syn_z","mis_z","pLI")

all_tissue_expression_constraint = function(gtex, tissues) {
  ## Association between expression and constraint
  # Association statistics
  res <- NULL
  for (i in metrics){
    res <- rbind(res,cbind(i,t(apply(gtex[,which(!is.na(match(colnames(gtex),tissues)))],2,function(x) summary(lm(log2(x[which(x>0)])~gtex[which(x>0),i]))$coefficients[2,]))))
  }
  # Prepare plot
  pdf('figures/GTEx_expression_constraint_all_tissues.pdf',height=9.50,width=9.50)
  nwide <- 2; nhigh <- 4; k <- 0
  par(mfrow=c(nhigh,nwide))
  par(oma=c(3.5,3.5,1,1))
  par(mar=c(3,2,4,1))
  # Loop through tissues
  for (tissue in tissues){
    # Choose genes with RPKM > 0
    ind <- which(gtex[,tissue]>0)
    # Select data
    list1 <- split(log2(gtex[ind,tissue]),gtex$syn_cut[ind])
    list2 <- split(log2(gtex[ind,tissue]),gtex$mis_cut[ind])
    list3 <- split(log2(gtex[ind,tissue]),gtex$lof_phi_cut[ind])
    # Plot data
    myvioplot(list1,list2,list3,color_syn,color_mis,color_lof,1.25)
    # Add axes
    axis(2,at=seq(-10,20,2),las=2,lwd=0,lwd.tick=1,line=0.25)
    mtext(gsub("_"," ",tissue),side=3,line=1,cex=1.2)
    k <- k+1
    if (k %% nwide == 1) {mtext('log2(RPKM)',side=2,line=3,cex=1)}
    if (k>(nwide*(nhigh-1))&k<=(nhigh*nwide)) {mtext('Gene constraint',side=1,line=3,cex=1)}
    if (k<=(nwide*(nhigh-1))&!is.na(match(tissue,tissues[(length(tissues)-(length(tissues) %% nwide)+1):length(tissues)]))) {mtext('Gene constraint',side=1,line=3,cex=1)}
    if (k==(nhigh*nwide)){k <- 0}
  }
  dev.off()
  
  ## Number of tissues gene is expressed in, excluding individual brain subregions
  # Prepare plot
  pdf('figures/GTEx_NTissues_constraint.pdf',height=9.50,width=9.50)
  nwide <- 2; nhigh <- 4; k <- 0
  par(mfrow=c(nhigh,nwide))
  par(oma=c(3.5,3.5,1,1))
  par(mar=c(3,2,4,1))
  
  # Loop through RPKM cutoffs
  rpkms <- c(0.1,1,5,10)
  for (i in rpkms){
    # Select data
    data <- apply(gtex[,which(!is.na(match(colnames(gtex),tissues[grep("Brain_",tissues,invert=T)])))],1,function(x) length(which(x>i)))
    list1 <- split(data,gtex$syn_cut)
    list2 <- split(data,gtex$mis_cut)
    list3 <- split(data,gtex$lof_phi_cut)
    # Plot data
    myvioplot(list1,list2,list3,color_syn,color_mis,color_lof,1.25)
    # Add axes
    axis(2,at=seq(0,40,5),las=2,lwd=0,lwd.tick=1,line=0.25)
    mtext(sprintf('RPKM > %.1f',i),side=3,line=1,cex=1.2)
    k <- k+1
    if (k %% nwide == 1) {mtext('N(Tissues)',side=2,line=3,cex=1)}
    if (k>(nwide*(nhigh-1))&k<=(nhigh*nwide)) {mtext('Gene constraint',side=1,line=3,cex=1)}
    if (k<=(nwide*(nhigh-1))&!is.na(match(i,rpkms[(length(rpkms)-(length(rpkms) %% nwide)+1):length(rpkms)]))) {mtext('Gene constraint',side=1,line=3,cex=1)}
    if (k==(nhigh*nwide)){k <- 0}
    # Association statistics
    for (j in metrics){
      res <- rbind(res,c(j,summary(lm(data~gtex[,j]))$coefficients[2,]))
      rownames(res)[dim(res)[1]] <- sprintf('RPKM>%.1f',i)
    }
  }
  
  dev.off()
}

plot_expression_constraint = function(gtex, tissues, title=F, save_plot=T) {
  ## Plot RPKM > 1 separately
  # Initialize plot
  if (save_plot) {
    pdf('figures/GTEx_NTissues_constraint_RPKM>1.pdf',height=4,width=9.50)
    par(oma=c(2,3,1,1))
    par(mar=c(2.5,1.5,2.5,1))
  }
  # RPKM cutoff
  i <- 1
  # Select data
  data <- apply(gtex[,which(!is.na(match(colnames(gtex),tissues[grep("Brain_",tissues,invert=T)])))],1,function(x) length(which(x>i)))
  list1 <- split(data,gtex$syn_cut)
  list2 <- split(data,gtex$mis_cut)
  list3 <- split(data,gtex$lof_phi_cut)
  # Plot data
  myvioplot(list1,list2,list3,color_syn,color_mis,color_lof,1.25) 
  # Add axes
  axis(2,at=seq(0,40,5),las=2,lwd=0,lwd.tick=1,line=0.25)
  if (title) mtext(sprintf('RPKM > %.1f',i),side=3,line=1,cex=1.2)
  mtext('expressed in N tissues',side=2,line=3,cex=1.2)
  mtext('constraint bins',side=1,line=1,cex=1.2)
  if (save_plot) dev.off()
}

