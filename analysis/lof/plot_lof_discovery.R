library(plyr)

snapshots = seq(from=500, to=121000, by=500)/2
discovery_pops = c("exac", "afr", "eas", "nfe", "fin", "sas", "amr", "consanguineous", "sas_non_consang")

load_and_combine_samples <- function(path) {
  runs = list.files(path, full.names=TRUE)
  final = data.frame(matrix(0, nrow=length(snapshots), ncol=length(discovery_pops)))
  colnames(final) = discovery_pops
  for (i in 1:length(runs)) {
    popcurves = list.files(runs[i], full.names=TRUE)
    for (curve in popcurves) {
      popname = sub(".txt", "", basename(curve))
      final[,popname] = final[,popname] + read.table(curve, header=T)
    }
  }
  final = final/i
  return(final)
}

load_lof_discovery_curves <- function(gene=T, hom=T, rare=T) {
  # pops = c("afr", "amr", "eas", "fin", "nfe", "sas", "consanguineous", "sas_non_consang")
  criteria = sprintf('%s_%s_%s', ifelse(gene, 'gene', 'var'),
                 ifelse(rare, 'rare', 'all'),
                 ifelse(hom, 'hom', 'lof'))
  curves = load_and_combine_samples(sprintf("data/lof/discovery/%s", criteria))
  # curves = data.frame(sapply(pops, function(x) { get_curve(x, crit=criteria) }))
  curves[curves==0] = NA
  #curves = rbind(rep(0, ncol(curves)), curves)
  corr_factor = ifelse(!gene | !hom, 1000, 1)
  curves = curves/corr_factor
  return(curves)
}

plot_lof_discovery_curves = function(curves, save_plot=T, uniform=F, gene=T, hom=T, rare=T, exac_line=T, last_point=NA) {
  criteria = sprintf('%s_%s_%s', ifelse(gene, 'gene', 'var'),
                     ifelse(rare, 'rare', 'all'),
                     ifelse(hom, 'hom', 'lof'))
  if (save_plot) {
    pdf(paste0('figures/5_', criteria, '_discovery.pdf'), width=6, height=4)
    par(mar=c(5,6,1,1))
  }
  plot_pops = discovery_pops
  if (uniform) {
    plot_pops = c(discovery_pops, 'uniform')
  }
  if (!exac_line) {
    plot_pops = plot_pops[!grepl("exac", plot_pops)]
  }
  if (hom) {
    plot_pops = plot_pops[plot_pops != 'sas']
  } else {
    plot_pops = plot_pops[!grepl('consang', plot_pops)]
  }
  if (is.na(last_point)) {
    last_point = max(apply(curves, 2, function(x) { sum(!is.na(x))}))
  } else {
    if (!(last_point %in% snapshots)) stop(sprintf('Invalid arg: last_point = %s', last_point))
    last_point = which(snapshots == last_point)
  }
  snapshots = snapshots[1:last_point]
  curves = curves[1:last_point,plot_pops]
  xlabel = 'number of individuals'
  par(bty='n', las=1)
  plot(NA, NA, log='', xlim=c(0, max(snapshots)), ylim=c(0, max(curves, na.rm=T)*1.35), lwd=0, type='l', xlab=xlabel, ylab='', xaxs='i', yaxs='i')
  sapply(plot_pops, function(x) {
    lines(snapshots, curves[[x]], col=pop_colors[[x]], lwd=4, lty=ifelse(x == 'exac', "11", 1))
  })
  par(las=0)
  text_line = ifelse(gene & !hom, 2, 3)
  if (gene) {
    mtext('number of genes with', side=2, line=text_line + 1, padj=0)
    ylabel = paste('at least one', ifelse(hom, 'homozygous', ''), 'PTV variant')
  } else {
    ylabel = paste('number of', ifelse(hom, 'homozygous', ''), 'PTV variants')
  }
  if (!gene | !hom) ylabel = paste(ylabel, '(thousands)')
  mtext(ylabel, side=2, line=text_line, padj=0)
  legend("topleft", pop_names[plot_pops], fill=pop_colors[plot_pops], text.col=pop_colors[plot_pops],
       cex=0.8, ncol=2, bty='n', x.intersp=0.2, text.width=2250)#, inset=c(0.2, 0))
  
  results = sapply(plot_pops, function(x) { 
    lfit = lm(log(curves[[x]]) ~ log(snapshots[1:length(curves[[x]])]))$coefficients
    c(lfit[[1]], lfit[[2]])
  })
  if (save_plot) dev.off()
}


