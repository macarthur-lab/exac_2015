source('exac_constants.R')
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}
if (!("use_data" %in% ls(globalenv()))) {
  use_data = subset(exac, use)
}
if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}

if (!("cadd_data" %in% ls(globalenv()))) {
  cadd_data = load_cadd_exac()
}

category_comparison = function(exac, correction = TRUE) {
  print('Preparing data...')
  
  if (correction) {
    if ('ps_predicted_weighted' %in% colnames(exac)) {
      exac_mut = exac
    } else {
      mutation_data = load_mutation_data(exac)
      exac_mut = merge(exac, subset(mutation_data, select=c(context, ref, alt, mu_snp)))
      expected_proportion_singleton_lm = lm(singletons/n ~ mu_snp, mutation_data, weights = mutation_data$n)
      exac_mut$ps_predicted_weighted = predict(expected_proportion_singleton_lm, data.frame(mu_snp = exac_mut$mu_snp))
    }
    use_data_mut = subset(exac_mut, use & !is.na(mu_snp))
  } else {
    use_data_mut = subset(exac, use)
  }
  print('Loaded. Merging with CADD...')
  exac_cadd = merge(use_data_mut, cadd_data)
  
  print('Merged. Getting subsets...')
  missense = subset(exac_cadd, category == 'missense_variant')
  missense$cadd_break = cut(missense$rawscore, quantile(missense$rawscore, probs=seq(0, 1, 1/3)))
  missense$cadd = match(missense$cadd_break, levels(missense$cadd_break))
  
  stops = subset(exac_cadd, category == 'stop_gained')
  stops$cadd_break = cut(stops$rawscore, breaks=quantile(stops$rawscore, probs=seq(0, 1, 1/3)))
  stops$cadd = match(stops$cadd_break, levels(stops$cadd_break))
  
  print('Done. Calculating categories...')
  raw_cats = ddply(use_data_mut, 'category', function(category) {
    data.frame(subpanel='main', n=nrow(category), singletons=sum(category$singleton), correction=sum(category$ps_predicted_weighted))
  })
  pphen_cats = ddply(missense, 'polyphen_word', function(category) {
    data.frame(subpanel='polyphen', n=nrow(category), singletons=sum(category$singleton), correction=sum(category$ps_predicted_weighted))
  })
  sift_cats = ddply(missense, 'sift_word', function(category) {
    data.frame(subpanel='sift', n=nrow(category), singletons=sum(category$singleton), correction=sum(category$ps_predicted_weighted))
  })
  missense_cadd_cats = ddply(missense, 'cadd', function(category_cadd) {
    data.frame(subpanel='missense_cadd', n=nrow(category_cadd), singletons=sum(category_cadd$singleton), correction=sum(category_cadd$ps_predicted_weighted))
  })
  stop_cadd_cats = ddply(stops, 'cadd', function(category_cadd) {
    data.frame(subpanel='stop_cadd', n=nrow(category_cadd), singletons=sum(category_cadd$singleton), correction=sum(category_cadd$ps_predicted_weighted))
  })
  names(pphen_cats)[1] = 'category'
  names(sift_cats)[1] = 'category'
  names(missense_cadd_cats)[1] = 'category'
  names(stop_cadd_cats)[1] = 'category'
  stop_cadd_cats$category = c('low', 'medium', 'high', NA)
  missense_cadd_cats$category = c('low', 'medium', 'high', NA)
  
  categories = rbind(raw_cats, pphen_cats, missense_cadd_cats, stop_cadd_cats)
  return(categories)
}

filter_categories = function(categories, correction = TRUE, min_n_category = 500) {
  categories = subset(categories, !is.na(category) & n > min_n_category)
  # Preparing metrics and error bars
  if (correction) {
    categories$proportion_singletons_raw = categories$singletons/categories$n
    categories$proportion_singletons = (categories$singletons-categories$correction)/categories$n
    categories$sem_ps = sqrt(categories$proportion_singletons_raw*(1-categories$proportion_singletons_raw)/categories$n)
  } else {
    categories$proportion_singletons = categories$singletons/categories$n
    categories$sem_ps = sqrt(categories$proportion_singletons*(1-categories$proportion_singletons)/categories$n)
  }
  categories$lower95 = categories$proportion_singletons - categories$sem_ps*1.96
  categories$upper95 = categories$proportion_singletons + categories$sem_ps*1.96
  categories
  
  # Assigning colors
  categories$k = '#000000'
  categories$k[categories$category %in% lof_like] = k_lof
  categories$k[categories$category %in% mis_like] = k_mis
  categories$k[categories$category %in% syn_like] = k_syn
  categories$k[categories$subpanel=='polyphen'] = k_mis
  # categories$k[categories$subpanel=='sift'] = k_mis
  categories$k[categories$subpanel=='missense_cadd'] = k_mis
  categories$k[categories$subpanel=='stop_cadd'] = k_lof
  # categories = subset(categories, k != '#000000')
  
  # Sorting and setting positions
  plot_order = order(categories$proportion_singletons[categories$subpanel=='main'])
  plot_order = c(plot_order, max(plot_order) + order(categories$proportion_singletons[categories$subpanel=='polyphen']))
  plot_order = c(plot_order, max(plot_order) + order(categories$proportion_singletons[categories$subpanel=='missense_cadd']))
  plot_order = c(plot_order, max(plot_order) + order(categories$proportion_singletons[categories$subpanel=='stop_cadd']))
  
  categories$xval = NA
  categories$xval[plot_order] = 1:length(plot_order)
  
  # leave space between panels
  categories$xval[categories$subpanel=='polyphen'] = categories$xval[categories$subpanel=='polyphen'] + 1 
  # categories$xval[categories$subpanel=='sift'] = categories$xval[categories$subpanel=='sift'] + 2
  categories$xval[categories$subpanel=='missense_cadd'] = categories$xval[categories$subpanel=='missense_cadd'] + 2
  categories$xval[categories$subpanel=='stop_cadd'] = categories$xval[categories$subpanel=='stop_cadd'] + 3
  
  # Setting text parameters
  categories$display_text = format_vep_category(categories$category)
  
  return(categories)
}

plot_category_comparison = function(categories, correction = TRUE, save_plot = F, leg = F) {
  # Plotting parameters
  panel_line = c(max(categories$xval[categories$subpanel=='main']),
                 max(categories$xval[categories$subpanel=='polyphen']),
                 max(categories$xval[categories$subpanel=='missense_cadd'])) + 1
  panel_mids = c(mean(categories$xval[categories$subpanel=='main'],na.rm=TRUE),
                 mean(categories$xval[categories$subpanel=='polyphen'],na.rm=TRUE),
                 mean(categories$xval[categories$subpanel=='missense_cadd'],na.rm=TRUE),
                 mean(categories$xval[categories$subpanel=='stop_cadd'],na.rm=TRUE))
  
  text_bottom = max(panel_line) # max(categories$xval[categories$proportion_singletons > lof_prop_singletons]) + 2
  xrange = c(min(categories$xval,na.rm=TRUE)-.5, max(categories$xval,na.rm=TRUE)+.5)
  yrange = c(min(categories$lower95, na.rm=T) - 0.02, max(categories$upper95, na.rm=T)+.02)
  
  if (save_plot) {
    if (correction) {
      pdf('figures/2d_correct_proportion_singleton.pdf', width=8, height=4)
    } else {
      pdf('figures/s_recur_uncorrected_proportion_singleton.pdf', width=8, height=4)
    }
    par(mar=c(8, 4, 1.5, 1.5))
  }
  text_size = 1
  pt_size = 2
  
  plot(NA,NA,xlim=xrange,ylim=yrange,yaxs='i',xaxs='i',axes=FALSE,xlab='',ylab=ifelse(correction, 'MAPS', 'proportion singleton'))
  # mtext(side=2, text=ifelse(correction, 'MAPS score', 'proportion singleton'), line=4, cex=text_size)
  points(categories$xval,categories$proportion_singletons,pch=20,col=categories$k, cex=pt_size)
  for (i in categories$xval) {
    points(x=rep(categories$xval[i],2),y=c(categories$lower95[i],categories$upper95[i]),type='l',lwd=3,col=categories$k[i])
  }
  axis(side=2,at=(0:10)/10,labels=(0:10)/10,lwd=0,lwd.ticks=1,las=2,cex.axis=text_size)
  library(plotrix)
  staxlab(side=1, at=categories$xval, labels=categories$display_text, las=1, cex=text_size, col=categories$k, srt=45)
  abline(v=panel_line,lwd=1)
  abline(h=0, lwd=.5, lty=2)
  mtext(side=1,at=panel_mids,text=c('all variants','PolyPhen','missense', 'nonsense'), line=6, cex=text_size)
  mtext(side=1,at=panel_mids,text=c('', '','CADD', 'CADD'), line=7, cex=text_size)
  if (leg) legend('topleft', c('LoF', 'missense', 'other'), text.col=c(color_lof, color_mis, color_syn), bty='n')
  if (save_plot) dev.off()
}