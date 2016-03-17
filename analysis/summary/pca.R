plot_pca = function(x='2', y='3', save_plot=F) {
  pca <- read.csv("data/summary/pca_data.csv", header=T)
  if (save_plot) {
    pdf("figures/1_pca.pdf")
    par()
  }  
  x = paste0('PC', x)
  y = paste0('PC', y)
  par(las=1)
  plot(pca[[x]], pca[[y]], col=alpha(pca[['pop_color']], 0.85), cex=0.2, xlab=x, ylab='')
  par(las=0)
  mtext(y, side=2, line=3)
  if (save_plot) dev.off()
}