#GWAS constrained plots
library(reshape)
gwas<-read.table("data/constraint/gwas.genes.parent.uniq",header=T,sep="\t")
gwas$inGWAS=1

source("exac_constants.R")
if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}

#To use 1 -obs/exp uncomment below
#constraint$corrected_syn = 1 - constraint$n_syn/constraint$adj_exp_syn
#constraint$corrected_mis = 1 - constraint$n_mis/constraint$adj_exp_mis
#constraint$corrected_lof = 1 - constraint$n_lof/constraint$adj_exp_lof

prepare_gwas_data = function() {
  #Cuts to use for figure
  categories= c("syn_cut","mis_cut","lof_phi_cut")

  cuts = unique(sort(constraint[,'lof_phi_cut']))
  data=data.frame(matrix(NA, nrow = length(cuts), ncol = 3))
  names(data)<-categories
  sds=data.frame(matrix(NA, nrow = length(cuts), ncol = 3))
  names(sds)<-categories 
  pvals=data.frame(matrix(NA, nrow = 1, ncol = 3))
  names(pvals) <- categories 
  
  lowEnrichment=data.frame(matrix(NA, nrow = length(cuts), ncol = 3))
  names(lowEnrichment) <- categories
  highEnrichment=data.frame(matrix(NA, nrow = length(cuts), ncol = 3))
  names(highEnrichment) <- categories
  
  for (zset in categories) {
  	#constraint$cut=cut(constraint[,zset],breaks=c(quantile(na.rm=T,constraint[,zset],probs = cuts)),include.lowest=TRUE)
  	constraint$cut = constraint[,zset]
  
  	b<-merge(constraint,gwas,by.x="gene",by.y="Mapped_gene",all.x=T,all.y=F)
  
  	b[is.na(b$inGWAS),"inGWAS"]=0
  	b$PARENT<-as.character(b$PARENT)
  	b[is.na(b$PARENT),"PARENT"]<-""
  
  	complete<-b[complete.cases(b),]
  	cnondup<-complete[!duplicated(complete$gene),]	
  	stats<-data.frame(score=unique(sort(cnondup$cut)))
  	
  	#Number of GWAS hits per category 
  	stats$ALL<-tapply(cnondup$inGWAS,cnondup$cut,sum)
  
  	#Proportion of GWAS hits
  	stats$mean<-tapply(cnondup$inGWAS,cnondup$cut,mean,na.rm=T)
  	
  	#Total number of evaluated genes
  	stats$n<-tapply(cnondup$inGWAS,cnondup$cut,length)
  
  	#Genes not in GWAS catalog
  	stats$notGWAS = stats$n - stats$ALL
  
  	#scaling factor to compute enrichment
  	nf=sum(stats$n)/sum(stats$ALL)
  	
  	#Proportion of GWAS genes in bin scaled by enrichment factor 
  	data[,zset] = stats$mean*nf
  	
  	#SE of proportion
  	SE=sqrt(stats$mean*(1-stats$mean)/stats$n)
  	#multiply by scaling factor
  	SE_scaled=SE*nf
  
  	#Get 95% CIs
  	lowEnrichment[,zset]=data[,zset] - (1.96*SE_scaled)
  	highEnrichment[,zset]=data[,zset] + (1.96*SE_scaled)
  
  	
  	#Get high vs Low
  	cmp<-stats[c(3,1),c("ALL","notGWAS")]
  	pvals[,zset]=chisq.test(cmp)$p.value
  }
  write.table(file="data/GWAS_Constraint_pval.tsv",pvals,quote=F,row.names=F)
  return(list(data, lowEnrichment, highEnrichment))
}

plot_gwas_figure = function(gwas_data, save_plot=T) {
  data = gwas_data[[1]]
  lowEnrichment = gwas_data[[2]]
  highEnrichment = gwas_data[[3]]
  if (save_plot) pdf("figures/Figure_constraint_GWAS.pdf",height=5,width=5)
	plot(NA,NA,xlim=c(0,2.1),ylim=c(0.5,1.5),yaxs='i',xaxs='i',axes=FALSE,xlab='constraint bins',ylab='GWAS enrichment', cex.lab=1.2)
	
	colors = c(k_syn, k_mis, k_lof)
	#comparison line
	#abline(h=1, lwd=.5, lty=2)
	axis(side=2,at=c(0.5,1.0,1.5),labels=paste(c('0.5','1.0','1.5'),'',sep=''),lwd=0,lwd.ticks=1,las=2,cex.axis=1.2)
	for (yval in seq(1,3,1)) {
  	#yvals = (length(data[,yval]):1)/length(data[,yval])
    yvals =seq(0.1,2.1,by=1)
    lines(yvals, data[,yval], col=rep(colors[yval],3), lwd=2.5)
    polygon(c(yvals, rev(yvals)), c(lowEnrichment[,yval], rev(highEnrichment[,yval])), col = alpha(rep(colors[yval],3), 0.2), border = NA)
	}
	lablist=c("low", "medium", "high")
	axis(1, at=seq(0.1, 2.1, by=1), labels = FALSE,tick=FALSE,lwd=0,lwd.ticks=1)
	text(seq(0.1, 2.1, by=1), par("usr")[3] - 0.02, labels = lablist, pos = 1, xpd = TRUE,col="black",las=2,cex=1.2)
	if (save_plot) dev.off()
}


