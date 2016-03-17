source('../../exac_constants.R')

if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}

gtex_path = 'data/constraint/gtex/'
transcripts <- read.table(sprintf('%scanonical_transcripts_v75.txt.gz',gtex_path),header=F)
colnames(transcripts) <- c('GeneID','TranscriptID')
constraint_gene <- merge(transcripts,constraint,by.x=2,by.y="feature")

### Function for eQTL plots
eqtlplot <- function(eqtl,se,colors,maxval,minval){
  par(lwd=2.5)
  n <- length(eqtl[1,])
  plot(1:n,eqtl[1,],type='l',col=colors[1],ylim=c(minval,maxval),axes=F,ylab='',xlab='',xlim=c(1,n))
  lines(1:n,eqtl[2,],col=colors[2])
  lines(1:n,eqtl[3,],col=colors[3])
  
  polygon(c(1:n, rev(1:n)),c(eqtl[1,]-se[1,]*1.96, rev(eqtl[1,]+se[1,]*1.96)), col = alpha(colors[1], 0.2), border = NA)
  polygon(c(1:n, rev(1:n)),c(eqtl[2,]-se[2,]*1.96, rev(eqtl[2,]+se[2,]*1.96)), col = alpha(colors[2], 0.2), border = NA)
  polygon(c(1:n, rev(1:n)),c(eqtl[3,]-se[3,]*1.96, rev(eqtl[3,]+se[3,]*1.96)), col = alpha(colors[3], 0.2), border = NA)
}

load_tissue = function(tissue) {
  # Genes included in eQTL analysis
  expr <- as.matrix(read.table(sprintf('%s%s.expr.genes.txt',gtex_path,tissue),header=T,as.is=T))
  # Genes with significant eQTL
  eqtl <- as.matrix(read.table(sprintf('%s%s.eqtl.genes.txt',gtex_path,tissue),header=T,as.is=T))
  
  gtex <- matrix(0,dim(constraint_gene)[1],2)
  rownames(gtex) <- as.character(constraint_gene[,2])
  gtex[which(!is.na(match(rownames(gtex),expr))),1] <- 1
  gtex[which(!is.na(match(rownames(gtex),eqtl))),2] <- 1
  
  peqtl <- NULL
  
  for (zcut in c('syn_cut','mis_cut','lof_phi_cut')){
    # Number of eQTL genes vs. constraint bin
    ntable <- table(gtex[gtex[,1]==1,2],constraint_gene[gtex[,1]==1,zcut])
    peqtl <- rbind(peqtl,as.numeric(ntable[2,]/colSums(ntable)))
    # Chisq test comparing low and high constraint bins
    pval <- c(tissue,zcut,chisq.test(ntable[,-2])$p.value)
  }
  
  # Standard errors for proportion
  peqtl_se <- sqrt((peqtl*(1-peqtl))/colSums(ntable))
  
  # % eQTLs in the most constrained bin compared to average
  tmp <- peqtl[,3]/as.numeric(table(gtex[gtex[,1]==1,2])[2]/sum(table(gtex[gtex[,1]==1,2])))
  print(sprintf('%s, synonymous: %.2f%%',tissue,tmp[1]*100))
  print(sprintf('%s, missense: %.2f%%',tissue,tmp[2]*100))
  print(sprintf('%s, loss-of-function: %.2f%%',tissue,tmp[3]*100))
  
  # Scale eQTL proportions by the ratio of total number genes and the total number of eQTLs
  neqtl <- sum(gtex[,1])/sum(gtex[,2])*peqtl
  neqtl_se <- sum(gtex[,1])/sum(gtex[,2])*peqtl_se
  return(list(neqtl, neqtl_se, pval))
}

eqtl_all_tissues = function() {
  # Tissues to be included
  # GTEx V4 eQTL tissues
  tissues <- c("Adipose_Subcutaneous", "Artery_Aorta", "Artery_Tibial", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Left_Ventricle", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Skin_Sun_Exposed_Lower_leg", "Stomach", "Thyroid", "Whole_Blood")
  
  # Initialize plot
  pdf('figures/GTEx_eQTL_constraint_all_tissues.pdf',height=9.50,width=9.50)
  
  # Plot layout
  nwide <- 3; nhigh <- 3; k <- 0
  par(mfrow=c(nwide,nhigh))
  par(oma=c(3.5,3.5,1,1))
  par(mar=c(3,2,4,1))
  
  pval <- NULL
  
  # Loop through tissues
  for (tissue in tissues){
  
    eqtl_data = load_tissue(tissue)
    neqtl = eqtl_data[[1]]
    neqtl_se = eqtl_data[[2]]
    pval = rbind(pval, eqtl_data[[3]])
    # Plot results
    eqtlplot(neqtl,neqtl_se,c(color_syn,color_mis,color_lof),1.65,0.35)
    
    # Add axes
    axis(2,at=seq(0.5,1.5,0.5),labels=seq(0.5,1.5,0.5),las=2,lwd=0,lwd.tick=1,line=0.25,cex.axis=1)
    mtext(c('low','medium','high'),side=1,line=0.5,at=c(1:3),cex=0.8)
    mtext(gsub("_"," ",tissue),side=3,line=1,cex=1.2)
    k <- k+1
    if (k %% nwide == 1) {mtext('eQTL enrichment',side=2,line=3,cex=1)}
    if (k>(nwide*(nhigh-1))&k<=(nhigh*nwide)) {mtext('Gene constraint',side=1,line=3,cex=1)}
    if (k<=(nwide*(nhigh-1))&!is.na(match(tissue,tissues[(length(tissues)-(length(tissues) %% nwide)+1):length(tissues)]))) {mtext('Gene constraint',side=1,line=3,cex=1)}
    if (k==(nhigh*nwide)){k <- 0}
  }
  
  dev.off()

  # Save P-values
  colnames(pval) <- c('Tissue','Constraint','P-value')
  write.table(pval,'data/GTEx_eQTL_constraint_pval.txt',col.names=T,row.names=F,quote=F)
}

# Plot whole blood
plot_eqtl_tissue = function(tissue='Whole_Blood', title=F, save_plot=T) {
  # Initialize plot
  if (save_plot) {
    pdf(paste0('figures/GTEx_eQTL_constraint_', tissue, '.pdf'),height=5,width=5)
  }
  eqtl_data = load_tissue(tissue)
  neqtl = eqtl_data[[1]]
  neqtl_se = eqtl_data[[2]]
  # Plot results and add axes
  eqtlplot(neqtl,neqtl_se,c(color_syn,color_mis,color_lof),1.5,0.5)
  axis(2,at=seq(0.5,1.5,0.5),labels=paste(c('0.5','1.0','1.5'),'',sep=''),las=2,lwd=0,lwd.tick=1,line=0.25,cex.axis=1.2)
  mtext(c('low','medium','high'),side=1,line=0.5,at=c(1:3),cex=1.2)
  if (title) mtext(gsub("_"," ",tissue),side=3,line=1)
  mtext('eQTL enrichment',side=2,line=3, cex=1.2)
  mtext('constraint bins',side=1,line=3, cex=1.2)
  if (save_plot) dev.off() 
}

