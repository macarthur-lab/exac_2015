# """
# Functions for plotting cpg-ti-tv recurrence curves.
# """

# source('exac_constants.R')
# This function assembles the sampling curves from a directory containing RData objects for each snapshots.
# It needs to be run once per sampling run. Each run yields a different set of curves for a particular sampling run.
load_mutation_class_samples <- function() {
  snapshots = seq(from=500, to=125000, by=500)
  df = data.frame(ti=0, tv=0, cpg=0)
  for (n in snapshots) {
    fname = sprintf("an_%s.RData", n)
    f = file.path('data/recurrence/ti_tv_cpg', fname)
    row = get(load(f))
    row[row==0] = -Inf
    df = rbind(df, row)
  }
  return(df)
}

plot_curves <- function(curves, with_derivative=F) {
  snapshots = seq(from=0, to=125000, by=500)
  matplot(snapshots, curves, type='l', col=c('red', 'green', 'blue'), 
          lty=1, lwd=2, xlab='an', ylab='# variants')
  if (with_derivative) {
    source('~/Desktop/MacArthur/ExAC_analysis/src/R/derivative.R')
    d = apply(curves, 2, derivative)
    par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
    colnames(d) = paste("d", colnames(d), sep='_')
    par(new = T)
    matplot(snapshots[2:length(snapshots)], d, type = "l", 
            lty=2, lwd=2, axes = FALSE, bty = "n", xlab = "", 
            ylab = "", col=c('red', 'green', 'blue'))
    axis(side=4, at=pretty(range(0:100)))
    mtext("derivative", side=4, line=3)
    par(new = F)
  }
  legend("topleft", colnames(curves), fill=c('red', 'green', 'blue'),cex=0.75, ncol=1, inset=0.1)
}

plot_mutation_class_samples = function(curves, save_plot=T, leg=T) {
  if (save_plot) {
    pdf('figures/2a_cpg_plateau.pdf', width=6, height=4)
    par(mar=c(4,5,1,1), cex=1.25)
  }
  par(bty='n', las=1)
  snapshots = seq(from=0, to=125000, by=500)
  ticks = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 120500)
  plot(snapshots, curves$ti/1000, type='l', col=color_ti, lwd=4, xlab='sampled number of chromosomes', ylab='number of variants (thousands)', xaxt='n',
       log='xy', ylim=c(10000, 1000000)/1000)
  axis(1, at=ticks, labels=format(ticks, scientific=FALSE))
  lines(snapshots, curves$tv/1000, type='l', col=color_tv, lwd=4)
  lines(snapshots, curves$cpg/1000, type='l', col=color_cpg, lwd=4)
  if (leg) legend("topleft", c('transversion', 'non-CpG transition', 'CpG transition'), lwd=4, col=c(color_tv, color_ti, color_cpg), text.col=c(color_tv, color_ti, color_cpg), bty='n')
  if (save_plot) dev.off()
}


plot_ti_tv = function(curves, save_plot=F, leg=T) {
  if (save_plot) {
    pdf('figures/s_recur_titv.pdf', width=6, height=4)
    par(mar=c(4,5,1,1), cex=1.25)
  }
  par(bty='n', las=1)
  snapshots = seq(from=0, to=125000, by=500)
  ticks = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 120500)
  plot(snapshots, (curves$ti + curves$cpg)/curves$tv, log='x', type='l', lwd=4, xlab='sampled number of chromosomes', ylab='TiTv', xaxt='n')
  axis(1, at=ticks, labels=format(ticks, scientific=FALSE))
  if (save_plot) dev.off()
}


