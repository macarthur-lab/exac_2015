source('exac_constants.R')
# source('analysis/properties/category_stats.R')
if (!("exac_canonical" %in% ls(globalenv()))) {
  exac_canonical = load_exac_data('canonical')
}
if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}
if (!("cadd_data" %in% ls(globalenv()))) {
  cadd_data = load_cadd_exac()
}

get_exac_with_mutation_data = function() {
  exac_constraint = merge(subset(exac_canonical, use), constraint, by='feature')
  mutation_data = load_mutation_data(exac_constraint)
  exac_mut = merge(exac_constraint, subset(mutation_data, select=c(context, ref, alt, mu_snp)))
  exac_mut$ps_predicted_weighted = predict(lm(singletons/n ~ mu_snp, mutation_data, weights = mutation_data$n), data.frame(mu_snp = exac_mut$mu_snp))
  return(exac_mut)
}

prepare_category_constraint_data = function(exac_mut) {
  print('Merging with CADD...')
  exac_cadd = merge(exac_mut, cadd_data)
  
  print('Merged. Getting subsets...')
  missense = subset(exac_cadd, category == 'missense_variant')
  missense$cadd_break = cut(missense$rawscore, quantile(missense$rawscore, probs=seq(0, 1, 1/3)))
  missense$cadd = match(missense$cadd_break, levels(missense$cadd_break))
  stops = subset(exac_cadd, category == 'stop_gained')
  stops$cadd_break = cut(stops$rawscore, breaks=quantile(stops$rawscore, probs=seq(0, 1, 1/3)))
  stops$cadd = match(stops$cadd_break, levels(stops$cadd_break))
  
  # Calculating %singleton across categories
  
  print('Done. Calculating categories...')
  mis_metric = 'mis_cut'
  lof_metric = 'lof_phi_cut'
  cat_detail = ddply(exac_mut, 'category', function(category) {
    constraint_metric = 'syn_cut'
    if (unique(category$category) %in% lof_like) {
      constraint_metric = lof_metric
    } else if (unique(category$category) %in% mis_like) {
      constraint_metric = mis_metric
    }
    category$constraint = match(category[,constraint_metric], levels(category[,constraint_metric]))
    ddply(subset(category, !is.na(constraint)), 'constraint', function(category_constraint) {
      data.frame(subpanel='main', n=nrow(category_constraint), singletons=sum(category_constraint$singleton), correction=sum(category_constraint$ps_predicted_weighted))
    })
  })
  
  print('Done. Calculating PolyPhen categories...')
  pphen_detail = ddply(exac_mut, 'polyphen_word', function(category) {
    category$constraint = match(category[,mis_metric], levels(category[,mis_metric]))
    ddply(subset(category, !is.na(mis_cut)), mis_metric, function(category_constraint) {
      data.frame(subpanel='polyphen', n=nrow(category_constraint), singletons=sum(category_constraint$singleton), correction=sum(category_constraint$ps_predicted_weighted))
    })
  })
  print('Done. Calculating SIFT categories...')
  sift_detail = ddply(exac_mut, 'sift_word', function(category) {
    category$constraint = match(category[,mis_metric], levels(category[,mis_metric]))
    ddply(subset(category, !is.na(mis_cut)), mis_metric, function(category_constraint) {
      data.frame(subpanel='sift', n=nrow(category_constraint), singletons=sum(category_constraint$singleton), correction=sum(category_constraint$ps_predicted_weighted))
    })
  })
  print('Done. Calculating CADD categories...')
  missense_cadd_cats_detail = ddply(missense, 'cadd', function(category) {
    category$constraint = match(category[,mis_metric], levels(category[,mis_metric]))
    ddply(subset(category, !is.na(mis_cut)), mis_metric, function(category_constraint) {
      data.frame(subpanel='missense_cadd', n=nrow(category_constraint), singletons=sum(category_constraint$singleton), correction=sum(category_constraint$ps_predicted_weighted))
    })
  })
  stop_cadd_cats_detail = ddply(stops, 'cadd', function(category) {
    category$constraint = match(category[,lof_metric], levels(category[,lof_metric]))
    ddply(subset(category, !is.na(mis_cut)), lof_metric, function(category_constraint) {
      data.frame(subpanel='stop_cadd', n=nrow(category_constraint), singletons=sum(category_constraint$singleton), correction=sum(category_constraint$ps_predicted_weighted))
    })
  })
  
  names(pphen_detail)[1:2] = c('category', 'constraint')
  names(sift_detail)[1:2] = c('category', 'constraint')
  names(missense_cadd_cats_detail)[1:2] = c('category', 'constraint')
  names(stop_cadd_cats_detail)[1:2] = c('category', 'constraint')
  stop_cadd_cats_detail$category = paste0('stop_', c(rep('low', 3), rep('medium', 3), rep('high', 3), NA))
  missense_cadd_cats_detail$category = paste0('missense_', c(rep('low', 3), rep('medium', 3), rep('high', 3), NA))
  categories_detail = subset(rbind(cat_detail, pphen_detail, missense_cadd_cats_detail, stop_cadd_cats_detail), !is.na(category) & n > 200)
  categories_detail$proportion_singletons_raw = categories_detail$singletons/categories_detail$n
  categories_detail$proportion_singletons = (categories_detail$singletons-categories_detail$correction)/categories_detail$n
  categories_detail$sem_ps = sqrt(categories_detail$proportion_singletons_raw*(1-categories_detail$proportion_singletons_raw)/categories_detail$n)
  categories_detail$lower95 = categories_detail$proportion_singletons - categories_detail$sem_ps*1.96
  categories_detail$upper95 = categories_detail$proportion_singletons + categories_detail$sem_ps*1.96
  return(categories_detail)
}

plot_category_stats_constraint = function(categories, categories_detail, save_plot=T) {
  if (!any(grepl('missense_', categories$category[categories$subpanel == 'missense_cadd']))) {
    categories$category[categories$subpanel == 'missense_cadd'] = paste0('missense_', categories$category[categories$subpanel == 'missense_cadd'])
    categories$category[categories$subpanel == 'stop_cadd'] = paste0('stop_', categories$category[categories$subpanel == 'stop_cadd'])
  }
  
  # Plotting
  panel_line = c(max(categories$xval[categories$subpanel=='main']),
                 max(categories$xval[categories$subpanel=='polyphen']),
                 max(categories$xval[categories$subpanel=='missense_cadd'])) + 1
  panel_mids = c(mean(categories$xval[categories$subpanel=='main'],na.rm=TRUE),
                 mean(categories$xval[categories$subpanel=='polyphen'],na.rm=TRUE),
                 mean(categories$xval[categories$subpanel=='missense_cadd'],na.rm=TRUE),
                 mean(categories$xval[categories$subpanel=='stop_cadd'],na.rm=TRUE))
  text_bottom = max(panel_line)
  xrange = c(min(categories$xval,na.rm=TRUE)-.5, max(categories$xval,na.rm=TRUE)+.5)
  yrange = c(min(categories$lower95, categories_detail$lower95) - 0.05, max(categories$upper95, categories_detail$upper95)+.05)
  
  if (save_plot) { 
    pdf('figures/3e_maps_constraint.pdf',width=9,height=6)
    par(mar=c(13,5,1,1), xpd=NA)
  }
  delta = ifelse(save_plot, 0.6, 1.2)
  s = -1; m = 7; l = 15
  plot(NA,NA,xlim=xrange,ylim=yrange,yaxs='i',xaxs='i',axes=FALSE,xlab='',ylab='MAPS')
  library(plotrix)
  staxlab(side=1, at=categories$xval, labels=categories$display_text, las=1, col=categories$k, srt=45)
  abline(v=panel_line,lwd=1)
  abline(h=0, lwd=.5, lty=2)
  text_size = 1
  mtext(side=1,at=panel_mids,text=c('All variants', 'PolyPhen','Missense', 'Nonsense'), line=5, cex=text_size)
  mtext(side=1,at=panel_mids,text=c('', '','CADD', 'CADD'), line=6, cex=text_size)
  axis(side=2,at=(0:10)/10,labels=(0:10)/10,lwd=0,lwd.ticks=1,las=2,cex.axis=1.2)
  
  d_ply(categories_detail, 'category', function(x) {
    this_cat = categories[categories$category == unique(x$category),]
    xval = this_cat$xval
    if (length(xval)) {
      xvals = xval - 0.6 + 0.9*(1:length(x$proportion_singletons))/length(x$proportion_singletons)
      points(xvals, x$proportion_singletons, col=this_cat$k, pch=20)
      lines(xvals, x$proportion_singletons, col=this_cat$k, lwd=2)
      polygon(c(xvals, rev(xvals)), c(x$lower95, rev(x$upper95)), col = alpha(this_cat$k, 0.2), border = NA)
    }
  })
  
  # Additional legends
  par(xpd=NA)
  delta_delta = ifelse(save_plot, 0.04, 0.1)
  legend_data = subset(categories_detail, category == '5_prime_UTR_variant')
  xvals = s:(s+2)
  points(xvals, legend_data$proportion_singletons - delta, col=color_syn, pch=20, cex=2)
  lines(xvals, legend_data$proportion_singletons - delta, col=color_syn, lwd=2)
  text(xvals, max(legend_data$proportion_singletons) - delta + delta_delta, col=color_syn, labels=c('<25%', '25-75%', '>75%'), cex=0.7)
  polygon(c(xvals, rev(xvals)), c(legend_data$lower95, rev(legend_data$upper95)) - delta, col = alpha(color_syn, 0.2), border = NA)
  xvals = m:(m+2)
  points(xvals, legend_data$proportion_singletons - delta, col=color_mis, pch=20, cex=2)
  lines(xvals, legend_data$proportion_singletons - delta, col=color_mis, lwd=2)
  text(xvals, max(legend_data$proportion_singletons) - delta + delta_delta, col=color_mis, labels=c('<25%', '25-75%', '>75%'), cex=0.7)
  polygon(c(xvals, rev(xvals)), c(legend_data$lower95, rev(legend_data$upper95)) - delta, col = alpha(color_mis, 0.2), border = NA)
  xvals = l:(l+2)
  points(xvals, legend_data$proportion_singletons - delta, col=color_lof, pch=20, cex=2)
  lines(xvals, legend_data$proportion_singletons - delta, col=color_lof, lwd=2)
  text(xvals, max(legend_data$proportion_singletons) - delta + delta_delta, col=color_lof, labels=c('<0.1', '0.1-0.9', '>0.9'), cex=0.7)
  polygon(c(xvals, rev(xvals)), c(legend_data$lower95, rev(legend_data$upper95)) - delta, col = alpha(color_lof, 0.2), border = NA)
  
  text(x=c(s+3, m+3, l+3), y=rep(-delta, 3), labels=c('Synonymous Z', 'Missense Z', 'LoF pLI'), adj=0)
  
  if (save_plot) dev.off()
}