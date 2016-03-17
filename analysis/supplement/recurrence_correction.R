source('../../exac_constants.R')
load_R_libraries( "plyr" )

if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}
if (!("use_data" %in% ls(globalenv()))) {
  use_data = subset(exac, use)
}
if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}

get_titv_singleton_by_category = function() {
  titv_categories = ddply(use_data, 'category', function(x) {
    cpgs = subset(x, cpg)
    tis = subset(x, transition & !cpg)
    tvs = subset(x, !transition)
    data.frame(ti_singletons=sum(tis$singleton), ti=nrow(tis),
               tv_singletons=sum(tvs$singleton), tv=nrow(tvs),
               cpg_singletons=sum(cpgs$singleton), cpg=nrow(cpgs))
  })
  titv_categories[titv_categories < 10] = 0
  titv_categories$singletons = titv_categories$ti_singletons + titv_categories$tv_singletons + titv_categories$cpg_singletons
  titv_categories$n = titv_categories$ti + titv_categories$tv + titv_categories$cpg
  titv_categories$prop_sing = titv_categories$singletons/titv_categories$n
  titv_categories$ti_prop_sing = titv_categories$ti_singletons/titv_categories$ti
  titv_categories$tv_prop_sing = titv_categories$tv_singletons/titv_categories$tv
  titv_categories$cpg_prop_sing = titv_categories$cpg_singletons/titv_categories$cpg
  
  titv_categories$sem = sqrt(titv_categories$prop_sing*(1-titv_categories$prop_sing)/titv_categories$n)
  titv_categories$ti_sem = sqrt(titv_categories$ti_prop_sing*(1-titv_categories$ti_prop_sing)/titv_categories$ti)
  titv_categories$tv_sem = sqrt(titv_categories$tv_prop_sing*(1-titv_categories$tv_prop_sing)/titv_categories$tv)
  titv_categories$cpg_sem = sqrt(titv_categories$cpg_prop_sing*(1-titv_categories$cpg_prop_sing)/titv_categories$cpg)
  
  titv_categories$lower95 = titv_categories$prop_sing - titv_categories$sem*1.96
  titv_categories$upper95 = titv_categories$prop_sing + titv_categories$sem*1.96
  titv_categories$ti_lower95 = titv_categories$ti_prop_sing - titv_categories$ti_sem*1.96
  titv_categories$ti_upper95 = titv_categories$ti_prop_sing + titv_categories$ti_sem*1.96
  titv_categories$tv_lower95 = titv_categories$tv_prop_sing - titv_categories$tv_sem*1.96
  titv_categories$tv_upper95 = titv_categories$tv_prop_sing + titv_categories$tv_sem*1.96
  titv_categories$cpg_lower95 = titv_categories$cpg_prop_sing - titv_categories$cpg_sem*1.96
  titv_categories$cpg_upper95 = titv_categories$cpg_prop_sing + titv_categories$cpg_sem*1.96
  return(titv_categories)
}

plot_titv_singleton_by_category = function(titv_categories, save_plot=T, leg=F, n_cutoff=200) {
  titv_categories = subset(titv_categories, !is.na(category) & ti > n_cutoff & tv > n_cutoff)
  
  titv_categories$k = '#000000'
  titv_categories$k[titv_categories$category %in% lof_like] = k_lof
  titv_categories$k[titv_categories$category %in% mis_like] = k_mis
  titv_categories$k[titv_categories$category %in% syn_like] = k_syn
  titv_categories$k[titv_categories$subpanel=='polyphen'] = k_mis
  titv_categories$k[titv_categories$subpanel=='sift'] = k_mis
  titv_categories = subset(titv_categories, k != '#000000')
  titv_categories$display_text = format_vep_category(titv_categories$category)
  
  plot_order = order(titv_categories$prop_sing)
  titv_categories$xval = NA
  titv_categories$xval[plot_order] = 1:length(plot_order)
  
  xrange = c(min(titv_categories$xval,na.rm=TRUE)-.5, max(titv_categories$xval,na.rm=TRUE)+.5)
  # yrange = c(min(titv_categories$cpg_lower95, na.rm=T) - 0.02, max(titv_categories$tv_upper95, na.rm=T)+.02)
  yrange = c(0, 0.9)
  
  text_size = 1
  pt_size = 1.8
  
  if (save_plot) {
    pdf('figures/s_recur_proportion_singletons_uncorrected_titv.pdf', width=7, height=5)
    par(mar=c(7,6,3,3), cex=1.2)
  }
  par(las=1)
  plot(NA,NA,xlim=xrange,ylim=yrange,yaxs='i',xaxs='i',xaxt='n',xlab='',ylab='proportion singleton')
#   points(titv_categories$xval, titv_categories$prop_sing, pch=20, col=titv_categories$k, cex=pt_size)
#   for (i in titv_categories$xval) {
#     lines(x=rep(titv_categories$xval[i], 2),y=c(titv_categories$lower95[i], titv_categories$upper95[i]), lwd=3, col=titv_categories$k[i])
#   }
  points(titv_categories$xval, titv_categories$ti_prop_sing, pch=15, col=color_ti, cex=pt_size)
  for (i in titv_categories$xval) {
    lines(x=rep(titv_categories$xval[i], 2),y=c(titv_categories$ti_lower95[i], titv_categories$ti_upper95[i]), lwd=3, col=color_ti)
  }
  points(titv_categories$xval, titv_categories$tv_prop_sing, pch=17, col=color_tv, cex=pt_size)
  for (i in titv_categories$xval) {
    lines(x=rep(titv_categories$xval[i], 2),y=c(titv_categories$tv_lower95[i], titv_categories$tv_upper95[i]), lwd=3, col=color_tv)
  }
  points(titv_categories$xval, titv_categories$cpg_prop_sing, pch=15, col=color_cpg, cex=pt_size)
  for (i in titv_categories$xval) {
    lines(x=rep(titv_categories$xval[i], 2),y=c(titv_categories$cpg_lower95[i], titv_categories$cpg_upper95[i]), lwd=3, col=color_cpg)
  }
  # axis(side=2,at=(0:10)/10,labels=paste((0:10)/10*100,'%',sep=''),lwd=0,lwd.ticks=1,las=2,cex.axis=par('cex'))
  # mtext(side=2, text='percent singleton', line=4)
  load_R_libraries( "plotrix" )
  staxlab(side=1, at=titv_categories$xval, labels=titv_categories$display_text, las=1, cex=text_size, col=titv_categories$k, srt=45)
  par(xpd=NA)
  if (leg) legend(0.5, 1, c('transversion', 'non-CpG transition', 'CpG transition'), pch=c(17, 15, 15), bty='n', col=c(color_tv, color_ti, color_cpg), text.col=c(color_tv, color_ti, color_cpg))
  if (save_plot) dev.off()
  
}

plot_mu_v_percent_singleton = function() {
  mutation_data = load_mutation_data(exac)
  pdf('figures/proportion_in_ac.pdf')
  par(bty='n', pch=20, mar=c(5,4,4,2), cex=1.5)
  log_type = ''
  plot(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n, ylab='Proportion singleton', xlab='Mu SNP', log=log_type, col=mutation_data$color)
  plot(mutation_data$mu_snp, mutation_data$doubletons/mutation_data$n, ylab='Proportion doubleton', xlab='Mu SNP', log=log_type, col=mutation_data$color)
  plot(mutation_data$mu_snp, mutation_data$tripletons/mutation_data$n, ylab='Proportion tripleton', xlab='Mu SNP', log=log_type, col=mutation_data$color)
  plot(mutation_data$mu_snp, mutation_data$quad/mutation_data$n, ylab='Proportion of AC == 4', xlab='Mu SNP', log=log_type, col=mutation_data$color)
  plot(mutation_data$mu_snp, mutation_data$quint/mutation_data$n, ylab='Proportion of AC == 5', xlab='Mu SNP', log=log_type, col=mutation_data$color)
  plot(mutation_data$mu_snp, mutation_data$ac_gt_five/mutation_data$n, ylab='Proportion of AC > 5', xlab='Mu SNP', log=log_type, col=mutation_data$color)
  dev.off()
  
  pdf('figures/proportion_ac_poster.pdf')
  par(bty='n', pch=20, mar=c(5,4,4,2), cex=1.5)
  plot(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n, ylab='Proportion singleton', xlab='Mu SNP', col=mutation_data$color, cex=1.5)
  abline(lm(singletons/n ~ mu_snp, mutation_data), lwd=3, lty=2)
  legend('topright', c('non-CpG transitions', 'transversions', 'CpG transition'), bty='n', pt.cex=1.5, col=c(color_ti, color_tv, color_cpg), pch=20)
  plot(mutation_data$mu_snp, mutation_data$not_rare/mutation_data$n, ylab='Proportion of AC >= 20', xlab='Mu SNP', col=mutation_data$color, cex=1.5)
  dev.off()
  
  # Another way
  all_ac_data = count(subset(synonymous_snps, select=c('context', 'ref', 'alt', 'ac_adj')))
  raw_numbers = count(subset(synonymous_snps, select=c('context', 'ref', 'alt')))
  colnames(raw_numbers) = c('context', 'ref', 'alt', 'total')
  
  all_ac_mutations = merge(merge(all_ac_data, mutations), raw_numbers)
  all_ac_mutations$frequency = all_ac_mutations$freq/all_ac_mutations$total
  
  #   a = subset(all_ac_mutations, ac_adj >= 20)
  #   b = ddply(a, c('context', 'ref', 'alt'), function(x) {
  #     data.frame(freq=sum(x$freq), total=unique(x$total), mu_snp=unique(x$mu_snp))
  #   })
  intercepts = ddply(all_ac_mutations, 'ac_adj', function(x) {
    data.frame(intercept=lm(x$freq/x$total ~ x$mu_snp)$coefficients[[1]],
               cpg=mean(subset(x, cpg)$frequency),
               n=sum(x$freq))
  })
  barplot(intercepts[1:10,]$intercept, space=0.1, col='black', names.arg=intercepts[1:10,]$ac_adj, xlab='AC_Adj', ylab='Intercept')
  
  # Difference in frequencies
  sum_diff = ddply(all_ac_mutations, 'ac_adj', function(x) {
    cpg = subset(x, cpg)
    noncpg = subset(x, !cpg)
    data.frame(cpg=sum(cpg$freq)/sum(cpg$total), noncpg=sum(noncpg$freq)/sum(noncpg$total))
  })
  plot(sum_diff$ac_adj, sum_diff$cpg - sum_diff$noncpg, pch=20, bty='n', xlim=c(1, 500), xlab='AC Adj', ylab='CpG - non_CpG')
  # Binarized
  plot(sum_diff$ac_adj, sum_diff$cpg - sum_diff$noncpg, pch=20, bty='n', xlim=c(1, 500), xlab='AC Adj', ylab='CpG > non_CpG', yaxt='n')
  axis(2, c(0, 1), c('FALSE', 'TRUE'))
  
  # sum_diff = ddply(all_ac_mutations, 'ac_adj', function(x) {
  #  lm(x$freq/x$total ~ x$mu_snp)
  #})
  
  return(mutation_data)
}

percent_singleton_regression = function() {
  mutation_data = correct_percent_singleton()
  pdf('figures/singleton_by_category.pdf')
  
  par(bty='n', pch=20, mar=c(5,4,4,2))
  # Linear
  plot(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n, ylab='Proportion singleton', xlab='Mu SNP', col=mutation_data$color, main='Linear')
  abline(lm(singletons/n ~ mu_snp, mutation_data))
  
  # Linear weighted
  plot(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n, ylab='Proportion singleton', xlab='Mu SNP', col=mutation_data$color, main='Weighted')
  abline(lm(singletons/n ~ mu_snp, mutation_data, weights = mutation_data$n))
  
  # LOESS
  plot(mutation_data$mu_snp, mutation_data$singletons/mutation_data$n, ylab='Proportion singleton', xlab='Mu SNP', col=mutation_data$color, main='LOESS')
  mu_spline = loess(singletons/n ~ mu_snp, mutation_data[order(mutation_data$mu_snp),], weights=mutation_data$n)
  mu_range = 1:150*8/1e10
  lines(mu_range, predict(mu_spline, data.frame(mu_snp=mu_range)), col = "darkblue")
  
  dev.off()
}