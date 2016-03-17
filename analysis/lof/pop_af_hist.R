source('../../exac_constants.R')

load_R_libraries( "plotrix" )

pop_af_hist <- function(df, pops, sing=T, doub=T, trip=T, breaks=c(-Inf,0.001,0.01,0.05,Inf), pop_spec=F, count_alleles=F) {
  pop_acs <- data.frame(df[,paste("ac",pops,sep="_")])
  ac_adj <- rowSums(pop_acs)
  pop_ans <- data.frame(df[,paste("an",pops,sep="_")])
  pop_afs <- data.frame(pop_acs/pop_ans)
  colnames(pop_afs) <- paste("af",pops,sep="_")
  
  num_bins <- length(breaks)-1 + sing + doub + trip
  spectrum <- matrix(0,nrow=num_bins, ncol=length(pops))
  first_bin <- 1 + sing + doub + trip
  bins <- apply(pop_afs, 2, function(afs) {cut(afs, breaks, labels=first_bin:num_bins)})
  unique <- apply(pop_acs, 2, function(x) { x == ac_adj })
  
  # will be modified if sing/doub/trip are set to True
  is_singleton <- rep(F, length(pop_acs))
  is_doubleton <- rep(F, length(pop_acs))
  is_tripleton <- rep(F, length(pop_acs))
  
  for (i in 1:num_bins) {
    if (sing) {
      is_singleton <- pop_acs == 1
      if (pop_spec) {
        is_singleton <- is_singleton & unique
      } 
      
      if (count_alleles) {
        spectrum[i,] <- colSums(pop_acs*is_singleton, na.rm=T)
      } else {
        spectrum[i,] <- colSums(is_singleton, na.rm=T) 
      }
      sing = F
    } 
    
    else if (doub) {
      is_doubleton <- pop_acs == 2
      if (pop_spec) {
        is_doubleton <- is_doubleton & unique
      }
      
      if (count_alleles) {
        spectrum[i,] <- colSums(pop_acs*is_doubleton, na.rm=T)
      } else {
        spectrum[i,] <- colSums(is_doubleton, na.rm=T) 
      }
      doub = F
    } 
    
    else if (trip) {
      is_tripleton <- pop_acs == 3
      if (pop_spec) {
        is_tripleton <- is_tripleton & unique
      }
      
      if (count_alleles) {
        spectrum[i,] <- colSums(pop_acs*is_tripleton, na.rm=T)
      } else {
        spectrum[i,] <- colSums(is_tripleton, na.rm=T) 
      }
      trip = F
    }
    
    else {
      bin <- bins == i & pop_acs > 0
      bin <- bin & !is_singleton & !is_doubleton & !is_tripleton
      if (pop_spec) {
        bin <- bin & unique
      }
      if (count_alleles) {
        spectrum[i,] <- colSums(pop_acs*bin, na.rm=T)
      } else {
        spectrum[i,] <- colSums(bin, na.rm=T) 
      }
    }
  }
  colnames(spectrum) = pops
  return(spectrum)
}

plot_pop_af_hist_db <- function(m, save_plot=T, normalize=F, log=F) {
  m = m[,pops]
  if (normalize) {
    m = m/m[,'nfe']
  } else if (log) {
    m = log10(m)
  }
  
  xticks <- c("singletons", "< 0.1%", "0.1-1%", "1-5%", ">5%")
  if (save_plot) {
    pdf('figures/s_lof_pop_af_hist.pdf', width=6, height=4)
    par(mar=c(3,4,1,1))
  }
  barplot(t(m), col=pop_colors[pops], beside=T, names.arg=xticks, ylab='PTV alleles per individual')
  legend('topleft', pop_names[pops], fill=pop_colors[pops], text.col=pop_colors[pops], cex=0.8, ncol=1, bty='n')
  if (save_plot) dev.off()
}

plot_pop_af_hist <- function(data, log=F, save_plot=T) {
  xticks <- c("singletons", "< 0.1%", "0.1-1%", "1-5%", ">5%")
  
  plot_data = rbind(t(data), NA)
  if (save_plot) {
    pdf('figures/5b_lofs_per_individual_freq.pdf', width=6, height=4)
    par(mar=c(3,4,1,1))
  }
  par(bty='n', las=1)
  if (log) {
    barplot(t(data), col=pop_colors[pops], beside=T, log='y', names.arg=xticks, ylab='PTVs')
  } else {
    centers = gap.barplot(plot_data, gap=c(21, 78), ytics=c(0,10,20,80,90,100), xaxt='n',
                col=rep(c(pop_colors[pops], 'white'), 5), xlab='', ylab='', yaxs='i')
    trans_centers = matrix(centers, nrow=nrow(plot_data))
    xticks_pos = colMeans(trans_centers[2:nrow(trans_centers),]) - 1
    staxlab(1, at=xticks_pos, lab=xticks, srt=45)
    par(las=0)
    mtext('PTV alleles per individual', side=2, line=3, padj=0)
  }
  legend('topleft', pop_names[pops], fill=pop_colors[pops], text.col=pop_colors[pops], ncol=1, bty='n', cex=0.8)
  if (save_plot) dev.off()
}



