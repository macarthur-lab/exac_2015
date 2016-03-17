library(gdata)
library(Hmisc)
library(reshape)

prepare_all_observed = function(type='use') {
  all_possible_observed = load_all_observed(type)
  return(subset_observed(all_possible_observed))
}

load_all_observed = function(type='use') {
  all_possible_cutoffs = read.table('data/summary/all_possible_cutoffs.txt.gz', header=T)
  observed_cutoffs = read.table(paste0('data/summary/observed_', ifelse(type != '', paste0(type, '_'), ''), 'cutoffs.txt.gz'), header=T)
  
  all_possible_observed = merge(all_possible_cutoffs, observed_cutoffs, by=c('consequence', 'cpg', 'transition', 'cutoff', 'cov'), suffixes=c('.poss', '.obs'))
  all_possible_observed$type = ifelse(all_possible_observed$transition, 'non-CpG transition', 'transversion')
  all_possible_observed$type[all_possible_observed$cpg] = 'CpG transition'
  return(all_possible_observed)
}

subset_observed = function(all_possible_observed) {
  possible_categories = subset(all_possible_observed, cutoff == 0.8 & cov == 10 & consequence %in% c('synonymous_variant', 'stop_gained', 'missense_variant'))
  possible_categories$consequence = gsub('stop_gained', 'nonsense', gsub('_variant', '', possible_categories$consequence))
  
  possible_categories$category <- reorder.factor(possible_categories$consequence, new.order=c('synonymous', 'missense', 'nonsense'))
  possible_categories$type <- reorder.factor(possible_categories$type, new.order=variant_types)
  possible_categories = arrange(possible_categories, category, type)
  possible_categories$proportion_observed = possible_categories$freq.obs/possible_categories$freq.poss
  possible_categories$sem = sqrt(possible_categories$proportion_observed*(1-possible_categories$proportion_observed)/possible_categories$freq.poss)
  possible_categories$upper = possible_categories$proportion_observed + 1.96*possible_categories$sem
  possible_categories$lower = possible_categories$proportion_observed - 1.96*possible_categories$sem
  return(possible_categories)
}

plot_proportion_observed = function(possible_categories, leg=F, save_plot=F) {
  plot_data = t(cast(possible_categories, category ~ type, value='proportion_observed'))
  uppers = t(cast(possible_categories, category ~ type, value='upper'))
  lowers = t(cast(possible_categories, category ~ type, value='lower'))
  errcol = 'black'
  spacing = c(0, 1)
  if (save_plot) pdf('figures/1d_proportion_observed.pdf', width=6, height=4)
  a = barplot(plot_data, beside=T, col=variant_type_colors, ylab='proportion possible found in ExAC', space=spacing, ylim=c(0,1), yaxt='n', cex.names=0.9)
  axis(side=2, at=(0:5)/5, labels=percent((0:5)/5), lwd=0, lwd.ticks=1, las=1)
  errbar(a, plot_data, lowers, uppers, add=T, col=errcol, errbar.col=errcol, lwd=2, pch=NA)
  if (leg) legend('topright', variant_types, col=variant_type_colors, text.col=variant_type_colors, pch=15, bty='n')
  if (save_plot) dev.off()
}