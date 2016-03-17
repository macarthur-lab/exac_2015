library(gdata)
library(ggplot2)
library(plotrix)

source('exac_constants.R')
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}
if (!("use_data" %in% ls(globalenv()))) {
  use_data = subset(exac, use)
}
if (!("raw_mutation_data" %in% ls(globalenv()))) {
  raw_mutation_data = load_raw_mutation_data()
}
if (!("mutation_data" %in% ls(globalenv()))) {
  mutation_data = load_mutation_data(use_data)
}

plot_mu_vs_singleton = function(save_plot=T, leg=F) {
  if (save_plot) {
    pdf('figures/s_recur_mu_vs_singleton.pdf', width=6, height=4)
    par(mar=c(4,4,2,1), cex=1.25)
  }
  par(bty='n', pch=20, las=1, xpd=F)
  plot(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n, ylab='proportion singleton', xlab='mutability', col=mutation_data$color)
  #cor.test(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n)
  abline(lm(singletons/n ~ mu_snp, mutation_data, weights=mutation_data$n), lwd=3, lty=2)
  if (leg) legend('topright', variant_types, bty='n', pt.cex=1.5, col=variant_type_colors, text.col=variant_type_colors, pch=20)
  if (save_plot) dev.off()
}

plot_recurrence_sfs = function(save_plot=T, leg=F) {
  mutation_data$type = ifelse(mutation_data$transition, 'non-CpG transition', 'transversion')
  mutation_data$type[mutation_data$cpg] = 'CpG'
  mutation_data$type <- reorder.factor(mutation_data$type, new.order=variant_types)
  collapsed_data = ddply(mutation_data, c('type'), function(x) {
    y = rbind(
      c('singletons', sum(x$singletons)/sum(x$n)),
      c('doubletons', sum(x$doubletons)/sum(x$n)),
      c('tripletons', sum(x$tripletons)/sum(x$n)),
      c('AC=4', sum(x$quad)/sum(x$n)),
      c('AC=5', sum(x$quint)/sum(x$n)),
      c('AC>5', sum(x$ac_gt_five)/sum(x$n))
    )
    data.frame(y)
  })
  names(collapsed_data) = c('type', 'freq', 'count')
  collapsed_data$count = as.numeric(collapsed_data$count)
  collapsed_data$freq <- reorder.factor(collapsed_data$freq, new.order=c('singletons', 'doubletons', 'tripletons', 'AC=4', 'AC=5', 'AC>5'))
  collapsed_data$type <- reorder.factor(collapsed_data$type, new.order=variant_types)
  font_size = 1.2
  if (save_plot) pdf('figures/2b_mutation_sfs.pdf', width=6, height=4)
  a = barplot(t(cast(collapsed_data, freq ~ type, value='count')), beside=T, col=c(color_tv, color_ti, color_cpg), ylab='proportion', xaxt='n', yaxt='n')
  labels = colnames(t(cast(collapsed_data, freq ~ type, value='count')))
  staxlab(side=1, at=a[2,], labels=labels, srt=45)
  axis(side=2, at=(0:6)/10, labels=percent((0:6)/10), lwd=0, lwd.ticks=1, las=1)
  if (leg) legend('topright', variant_types, bty='n', pch=15, col=variant_type_colors, text.col=variant_type_colors)
  if (save_plot) dev.off()
#   ggplot(collapsed_data) + geom_bar(aes(x=freq, fill=type, y=count), stat='identity', position='dodge', col='black') +
#     ylab('Count') + theme(axis.text=element_text(size = rel(font_size)), axis.text.x = element_text(angle = -45, hjust=0), axis.title=element_text(size = rel(font_size)), axis.title.x=element_blank(), legend.text=element_text(size = rel(font_size)), legend.title=element_text(size = rel(font_size)), panel.background = element_blank(), legend.position = 'bottom') +
#     scale_fill_manual(values=c('CpG'=color_cpg, 'non-CpG transition'=color_ti, 'transversion'=color_tv)) +
#     guides(fill = guide_legend(override.aes = list(colour = NULL)))
#   if (save_plot) {
#     pdf('figures/2b_mutation_sfs.pdf', width=6, height=4)
#     print(p)
#     dev.off()
#   } else {
#     return(p)
#   }
}

load_doubleton_data = function() {
  syn_doubletons = subset(use_data, ac_adj == 2 & category == 'synonymous_variant')
  d_syn = merge(syn_doubletons, subset(raw_mutation_data, select=c(context, ref, alt, mu_snp)), all.x=T)
  d_syn_filt = subset(d_syn, !((ac_nfe == 1 & ac_fin == 1) | (ac_nfe == 1 & ac_amr == 1) | (ac_nfe == 1 & ac_afr == 1)))
  doubleton_data = ddply(d_syn_filt, c('mu_snp', 'cpg', 'transition'), function(x) {
    data.frame(dd = mean(as.numeric(x$doubleton_dist), na.rm=T), cross_pop=sum(x$ac_popmax == 1, na.rm=T), total=nrow(x))
  })
  return(doubleton_data)
}

doubleton_cross_pop = function(doubleton_data, save_plot=F, leg=F, distance=F) {
  doubleton_data$color = ifelse(doubleton_data$transition, ifelse(doubleton_data$cpg, color_cpg, color_ti), color_tv)
  if (save_plot) {
    if (distance) {
      pdf('figures/s_recur_doubleton_distance.pdf', width=6, height=4)
    } else {
      pdf('figures/2c_recurrence_doubletons.pdf', width=6, height=4)
    }
    par(mar=c(4,5,1,1), cex=1.25)
  }
  par(bty='n', pch=20, las=1)
  if (distance) {
    plot(doubleton_data$mu_snp, doubleton_data$dd, ylab='', bty='n', pch=20, log='x', xlab='mutability', col=doubleton_data$color)
    par(las=0)
    mtext('mean Euclidean distance', side=2, line=4)
    cor.test(doubleton_data$mu_snp, doubleton_data$dd)
  } else {
    plot(doubleton_data$mu_snp, doubleton_data$cross_pop/doubleton_data$total, bty='n', pch=20, log='x', xlab='mutability', ylab='proportion cross-population', col=doubleton_data$color, yaxt='n', ylim=c(0, 0.5))
    axis(side=2, at=(0:5)/10, labels=percent((0:5)/10), lwd=0, lwd.ticks=1, las=1)
  }
  if (leg) legend('topleft', c('transversion', 'non-CpG transition', 'CpG'), bty='n', pt.cex=1.5, col=c(color_tv, color_ti, color_cpg), text.col=c(color_tv, color_ti, color_cpg), pch=20)
  if (save_plot) dev.off()
}

cross_pop_controls = function() {
  singletons = which(exac$singleton & exac$use & exac$category == 'synonymous_variant')
  cross_pops_singleton = replicate(10000, length(unique(exac$popmax[sample(singletons, 2)])) == 2) #0.7205
  sum(cross_pops_singleton)/10000 # 0.7205
  doubletons = which(exac$ac_adj == 2 & exac$ac_hom == 0 & exac$use & exac$category == 'synonymous_variant')
  ac_pops = paste0('ac_', pops)
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  check_doubleton = function() {
    dat = exac[sample(doubletons, 2),ac_pops]
    resample(which.max(dat[1,])[[1]]) != resample(which.max(dat[2,])[[1]])
  }
  cross_pops = replicate(10000, check_doubleton()) 
  sum(cross_pops)/10000 # 0.7407
}


