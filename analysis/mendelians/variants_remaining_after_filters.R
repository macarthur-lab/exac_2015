#### Begin preparation for Figure 4A

get_variants_remaining_after_filters = function() {
  # things to loop over in generating summary stats
  # allele frequency thresholds
  af_limits = c(1e-6, 10^-5.5, 1e-5, 10^-4.5, 1e-4, 10^-3.5, 1e-3, 10^-2.5, 1e-2, 5e-2)
  filter_xvals = c(1,2,-12,-11,3,4)
  # populations
  ancestries = names(pop_names)[c(1,2,3,5,7)]
  # filters 
  filters = c("esp_global_nonsing_filter", "esp_popmax_nonsing_filter", "kg_global_nonsing_filter", "kg_popmax_nonsing_filter", "exac_global_nonsing_filter","exac_popmax_nonsing_filter")
  filter_levels = c(rep(c("global","popmax"),3))
  filter_datasets = c(rep(c("ESP","1000 Genomes", "ExAC"),each=2))
  
  # if the "leave-out 500" data are available (which is true locally) AND the user wishes to
  # recalculate the stats for Figure 4A (which takes about 20 minutes on my MacBook Air)
  # then calculate the stats used for Figure 4A and also write a tab-delimited text file of the stats
  # to disk. If these data are not available (public release) then just read the pre-calculated stats
  # in from disk.
  recalculate_stats = FALSE
  if ('lv500.table.gz' %in% list.files('./misc_data/') & recalculate_stats) {
    lv = read.table('misc_data/lv500.table.gz',sep='\t',header=TRUE,comment.char='',quote='')
    colnames(lv) = tolower(colnames(lv))
    lv = lv[,-which(colnames(lv)=='ac_hemi')]
    lv = lv[,-which(colnames(lv)=='ac_hom')]
    
    # join relevant cols to ExAC
    poses = lv$pos[sample(1:dim(lv)[1],replace=FALSE,size=100)]
    lv$pos_id = paste(lv$chrom, formatC(lv$pos,width=9,flag='0'), lv$ref, lv$alt, sep='_')
    exac$pos_id = paste(exac$chrom, formatC(exac$pos,width=9,flag='0'), exac$ref, exac$alt, sep='_')
    
    match_indices = match(lv$pos_id, exac$pos_id)
    
    for (colname in colnames(lv)[7:25]) {
      # compute ACs in ExAC minus the leave-out 500
      lv[,paste("exac",colname,sep="_")] = exac[match_indices,colname] - lv[,colname]
    }
    
    # make allele frequency filters
    lv$kg_popmax_nonsing_filter = exac$kg_popmax_nonsing_filter[match_indices]
    lv$kg_global_nonsing_filter = exac$kg_global_nonsing_filter[match_indices]
    lv$esp_popmax_nonsing_filter = exac$esp_popmax_nonsing_filter[match_indices]
    lv$esp_global_nonsing_filter = exac$esp_global_nonsing_filter[match_indices]
    
    # re-calculate filters where necessary
    lv$exac_af_global = lv$exac_ac_adj / lv$exac_an_adj
    lv$exac_af_global[is.na(lv$exac_af_global)] = 0.0
    lv$exac_af_popmax = pmax(lv$exac_ac_afr / lv$exac_an_afr,
                             lv$exac_ac_amr / lv$exac_an_amr,
                             lv$exac_ac_eas / lv$exac_an_eas,
                             lv$exac_ac_nfe / lv$exac_an_nfe,
                             lv$exac_ac_sas / lv$exac_an_sas, na.rm=TRUE)
    lv$exac_af_popmax[is.na(lv$exac_af_popmax)] = 0.0
    lv$exac_popmax_nonsing_filter = lv$exac_af_popmax
    lv$exac_popmax_nonsing_filter[lv$exac_ac_adj==1] = 0.0
    lv$exac_global_nonsing_filter = lv$exac_af_global
    lv$exac_global_nonsing_filter[lv$exac_ac_adj==1] = 0.0
    
    lv$use = exac$use[match_indices] & lv$ac_adj > 0
    lv$category = exac$category[match_indices]
    
    # some constants and stats
    n_indiv = {} # count the number of ExAC individuals in each continental ancestry
    for (ancestry in ancestries) {
      n_indiv[[ancestry]] = max(lv[,paste("an",ancestry,sep="_")])/2
    }
    sum(n_indiv) # will be 500 
    
    # calculate summary stats
    stats = expand.grid(ancestries,filters,af_limits)
    colnames(stats) = c('ancestry','filter','af_limit')
    stats$ancestry = as.character(stats$ancestry)
    stats$total_vars = 0
    stats$filtered_vars = 0
    stats$proportion_removed = 0.0
    for (af_limit in af_limits) {
      for (filter in filters) {
        include_in_total = lv$use & !is.na(lv$category) & lv$category %in% c(mis_like,lof_like)
        include_in_filtered = include_in_total & (lv[,filter] < af_limit)
        for (ancestry in ancestries) {
          # find target row of stats table
          row = stats$ancestry == ancestry & stats$filter == filter & stats$af_limit == af_limit
          exac_colname = paste("ac",ancestry,sep="_") # column of interest in exac
          # below, these n_indiv terms cancel in the final calculation, but are included as a reminder
          # that conceptually what we are computing is the average number of variants _in each person_
          # of a given ancestry that are removed by an allele frequency filter
          average_variants_without_af_filter = sum(lv[include_in_total,exac_colname]) / n_indiv[ancestry]
          average_variants_with_af_filter = sum(lv[include_in_filtered,exac_colname]) / n_indiv[ancestry] 
          proportion_removed = 1 - average_variants_with_af_filter / average_variants_without_af_filter
          stats$proportion_removed[row] = proportion_removed # store it in the table
          stats$total_vars[row] = average_variants_without_af_filter
          stats$filtered_vars[row] = average_variants_with_af_filter
        }
      }
    }
    write.table(stats,'data/mendelians/filtering_stats.tsv',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
  } else {
    stats = read.table('data/mendelians/filtering_stats.tsv',sep='\t',quote='',header=TRUE)
  }
  
  # get colors and x axis values for plotting Figure 4A
  stats$acol = pop_colors[stats$ancestry]
  stats$xval = filter_xvals[match(stats$filter,filters)]
  
  return(stats)
}

plot_variants_remaining = function(stats, save_plot=F) {
  # Figure 4A
  if (save_plot) {
    pdf('figures/filtering_by_dataset.pdf',width=4,height=4)
    par(mar=c(6,6,1.5,1))
  }
  filters_to_plot = c("esp_global_nonsing_filter", "esp_popmax_nonsing_filter", "exac_global_nonsing_filter","exac_popmax_nonsing_filter")
  
  filter_xvals = c(1,2,-12,-11,3,4)
  filter_levels = c(rep(c("global","popmax"),3))
  
  plot(NA,NA,xlim=c(.5,4.5),ylim=c(0,1500),xlab='',ylab='',axes=FALSE,yaxs='i',xaxs='i')
  plot_data = subset(stats, af_limit == .001 & filter %in% filters_to_plot )
  points(plot_data$xval, plot_data$filtered_vars, pch=19, cex=2.5, col=alpha(plot_data$acol,.6))
  abline(h=0,lwd=2)
  axis(side=2, at=(0:3)*500, labels=(0:3)*500, lwd=2, lwd.ticks=1, las=1)
  mtext(side=1, at=filter_xvals, text=filter_levels, cex=1, line=.7)
  mtext(side=1, at=c(1.5, 3.5), text= c("ESP","ExAC"), cex=1.1, line=2.1, font=1)
  mtext(side=2, text='remaining predicted protein-altering variants\nper exome after 0.1% AF filter', line=3, cex=1.1, font=1)
  ancestries = names(pop_names)[c(1,2,3,5,7)]
  ancestry_names = pop_names[c(1,2,3,5,7)]
  text(x=rep(4.5,5),y=(15:11)*90,labels=ancestry_names,col=pop_colors[ancestries],font=1,cex=1.1,adj=c(1, 0))
  if (save_plot) dev.off()
}

#### End preparation for Figure 4A