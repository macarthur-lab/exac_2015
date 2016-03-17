# Anne's curation figure (4D) - data for this are just hard-coded here:

 curation_figure = function(save_plot = F) {
  LE = c(67,30,20,0)
  DE = c(25,2,0,0)
  BT = c(0,0,0,9)
  SDA = c(0,0,0,7)
  Variant_Review = (as.matrix(cbind(LE,DE,BT,SDA)))
  bt_genes = c('ABCA1','AGXT2','CD320','FAH','FUT2','FUT6','LIPG','LPA','SERPINA7')
  rm_genes = c('BTD','FLG','G6PD','GNE','PRSS1','SERPINA1')
  
  if (save_plot) {
    pdf('figures/curated_variants.pdf', width = 6, height = 6)
    par(mar=c(6,6,1.5,1))
  }
  ylim = c(0,125)
  colors = c("#0C97EC", "#0F8D1E", "#EC880C", "#000000")
  x_axis_labels = c("insufficient evidence", "classification error", "benign trait", "disease associated")
  classification = c("benign", "likely benign", "uncertain significance", "not reclassified")
  
  barplot(Variant_Review, names.arg = gsub(' ','\n',x_axis_labels), 
          border=NA, beside=FALSE, col = colors, cex.names = 0.8, main="", axes=FALSE,
          xlab = NA, ylab = NA, yaxs='i', xaxs='i', ylim=ylim)
  abline(h=0,lwd=3)
  axis(side=2, at = (0:5)*25, las = 2, lwd.ticks=1, lwd=2, cex=0.9) 
  mtext(side=2, "number of variants", line=3, font=1, cex=1.1)
  
  text(x=rep(3.1,5), y=(9:1)*7+3, pos=3, labels=bt_genes, font=3, cex=0.75) # third bar is from x = 2.6 to 3.6, midpoint is 3.1
  text(x=rep(4.3,2), y=(6:1)*7+3, pos=3, labels=rm_genes, font=3, cex=0.75) # fourth bar is 3.8 to 4.8, midpoint is 4.3
  
  legend('topright', legend=classification, col=colors, text.col=colors, text.font=1, cex=1, ncol=1, pch=15, bty='n')
  if (save_plot) dev.off()
}