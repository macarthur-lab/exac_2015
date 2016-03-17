### FIGURE 1C - INDEL DESCRIPTIVE STATS
library(binom)

coding_indels = c("inframe_indel", "frameshift_variant")

load_indel_data = function(exac) {
  indels = ddply(subset(exac, use & indel & category %in% coding_indels, select=c("indel","singleton","bases_inserted","af_global")), 'bases_inserted', function(x) {
    data.frame(n_variants=nrow(x), n_singletons=sum(x$singleton), mean_af=mean(x$af_global))
  })
  # add confidence intervals
  indels$lower95 = binom.confint(x=indels$n_singletons,n=indels$n_variants,method='asymptotic')$lower
  indels$upper95 = binom.confint(x=indels$n_singletons,n=indels$n_variants,method='asymptotic')$upper
  subset(indels, bases_inserted %in% c(-6:-1,1:6))
}

# Some indels non-coding:
# table(subset(exac, indel)$category)

# map values mod 3 to colors for plotting
indel_color = function(indel_size) {
  colors = character(length(indel_size))
  colors[indel_size %% 3 == 0] = k_mis
  colors[indel_size %% 3 != 0] = k_lof
  return (colors)
}

# what proportion variants are in this range
# sum(exac$bases_inserted %in% c(-6:-1,1:6))/sum(exac$bases_inserted != 0)

syn_prop_singletons = sum(exac$singleton & exac$use & exac$consequence=='synonymous_variant', na.rm=TRUE) / sum(exac$use & exac$consequence=='synonymous_variant', na.rm=TRUE)
lof_prop_singletons = sum(exac$singleton & exac$use & (!exac$indel) & exac$lof_use, na.rm=TRUE) / sum(exac$use & (!exac$indel) & exac$lof_use, na.rm=TRUE)

fs_prop_singletons = sum(exac$ac_adj==1 & exac$use & exac$indel & exac$category %in% coding_indels & exac$bases_inserted%%3 != 0, na.rm=TRUE)/sum(exac$use & exac$indel & exac$category %in% coding_indels & exac$bases_inserted%%3 != 0, na.rm=TRUE)
if_prop_singletons = sum(exac$ac_adj==1 & exac$use & exac$indel & exac$category %in% coding_indels & exac$bases_inserted%%3 == 0, na.rm=TRUE)/sum(exac$use & exac$indel & exac$category %in% coding_indels & exac$bases_inserted%%3 == 0, na.rm=TRUE)

plot_indel_hist = function(indels) {
  plot(indels$bases_inserted, indels$n_variants,type='h',col=indel_color(indels$bases_inserted),lwd=25,lend=1,
       axes=FALSE,xlab='',ylab='',xlim=c(-6,6),
       main='')
  abline(v=0,lwd=2,col='black')
  axis(side=1,at=-6:6,labels=c(-6:0,paste('+',1:6,sep='')),lwd=0,lwd.ticks=0)
  abline(h=0,lwd=2,col='black')
  mtext(side=1,adj=0,text='deletions')
  mtext(side=1,adj=1,text='insertions')
  par(xpd=NA)
  mtext('# of variants', side=2, line=4)
  par(xpd=F)
  axis(side=2, at=(0:4)*1e4, labels=(0:4)*1e4, las=1, cex.axis=.8, lwd=0, lwd.ticks=1, cex.axis=.8)
}

plot_indel_af = function(indels) {
  plot(indels$bases_inserted, indels$n_singletons/indels$n_variants,
       col=indel_color(indels$bases_inserted), pch=19, axes=FALSE,
       ylim=c(0.45,lof_prop_singletons+.2),xlim=c(-6,6), yaxs='i', xlab='', ylab='', yaxs='i')
  abline(v=0,lwd=2,col='black')
  abline(h=c(fs_prop_singletons,if_prop_singletons),lwd=2,col=c(k_lof,k_mis),lty=3)
  points(indels$bases_inserted, indels$n_singletons/indels$n_variants,col=indel_color(indels$bases_inserted), pch=19, cex=1.0)
  for (i in c(-6:-1,1:6)) {
    row = which(indels$bases_inserted==i)
    points(x=rep(i,2),y=c(indels$lower95[row],indels$upper95[row]),col=indel_color(indels$bases_inserted[row]),type='l',lwd=3)
  }
  axis(side=2, at=(0:10)/10, labels=paste((0:10)*10,"%",sep=""), las=1, col='black', lwd=0, lwd.ticks=1, cex.axis=.8)
  text(x=c(1,1),y=c(if_prop_singletons-.04,fs_prop_singletons-.04),labels=c('mean in-frame','mean frameshift'),col=c(k_mis,k_lof),cex=.8,pos=4)
  mtext(side=2, text='proportion', line=4)
  mtext(side=2, text='singletons', line=3)
}

plot_both_indels = function(indels, save_plot=F) {
  if (save_plot) {
    pdf('figures/indel_statistics.pdf',width=8,height=4)
    par(mfrow=c(2,1), mar=c(2,4,1,2))
  }
  plot_indel_hist(indels)
  plot_indel_af(indels)
  if (save_plot) dev.off()
}
