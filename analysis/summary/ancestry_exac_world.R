# Plot for main figure depicts size difference between ExAC and other sequencing databases
# Separate plot for supplement depicts representativeness of sequencing resources in terms of ancestry (All scaled to same height)
 
# Population numbers in millions from Wikipedia, according to the simplifying assumptions I made.

# Populations and percent ancestries  # World in millions, all others actual.
ancestries <- c("East Asian", "South Asian", "European", "Middle Eastern", "African", "Latino", "Other")
world.pop <- c(1932, 2085, 1145, 410, 1022, 529, 38)
kgenomes.pop <- c(523, 494, 514, 0, 691, 355, 0)
evs.pop <- c(0, 0, 4298, 0, 2217, 0, 0)
exac.pop <- c(max(exac$an_eas), max(exac$an_sas), max(exac$an_nfe) + max(exac$an_fin), 0, max(exac$an_afr), max(exac$an_amr), max(exac$an_oth))/2

exac.colors <- c(color_eas, color_sas, color_eur, color_mde, color_afr, color_amr, color_oth)

world.scaled.to.exac <- world.pop*sum(exac.pop)/sum(world.pop)
kgenomes.scaled.to.exac <- kgenomes.pop*sum(exac.pop)/sum(kgenomes.pop)
evs.scaled.to.exac <- evs.pop*sum(exac.pop)/sum(evs.pop)

plot_panel_size = function(save_plot) {
  if (save_plot) {
    pdf('figures/Fig_1A_ExAC_1KG_ESP.pdf',height=6,width=6)
    par(cex=1.1)
  }
  barplot(as.matrix(cbind(kgenomes.pop, evs.pop, exac.pop)), col=exac.colors,
          names.arg=c("1000 Genomes", "ESP", "ExAC"), border=NA, beside=FALSE, 
          main="", axes=FALSE)
  # abline(h=0,lwd=2)
  axis(side=2, at=(0:6)*1e4, labels=formatC((0:6)*1e4,big.mark=',',format='fg'), las=2, cex.axis=1, lwd.ticks=1, lwd=0)
  
  par(las=0)
  mtext("individuals", side=2, line=4)
  legend('topleft', bty="n", ncol=1, cex=0.9,
         legend=rev(ancestries), fill=rev(exac.colors), 
         border=rev(exac.colors), text.col=rev(exac.colors))
  
  if (save_plot) dev.off()
}

plot_panel_breakdown = function(save_plot=F) {
  if (save_plot) {
    pdf('figures/Supp_Ancestral_representativeness',height=4,width=6)
    par(mar=c(3,1,3,1))
  }
  
  barplot(as.matrix(cbind(kgenomes.scaled.to.exac, evs.scaled.to.exac, world.scaled.to.exac)), col=alpha(exac.colors, 1), cex.names=.8, yaxt="n", space=c(2,2,2,2,0,0)/10,
          border=NA, beside=FALSE, cex.main=.9, main="Ancestry Composition Comparison\n(all scaled to common bar height)", names.arg=c("1KG", "ESP","ExAC",  "World","", ""))
  abline(h=0,lwd=2)
  legend(x=5.05, y=58000, inset=c(2,20), bty="n", ncol=1, 
         legend=rev(ancestries), cex=.85, fill= rev(exac.colors), 
         border= rev(exac.colors), text.col=rev(exac.colors))
         
  if (save_plot) dev.off()
}


