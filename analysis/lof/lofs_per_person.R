source('../../exac_constants.R')

load_R_libraries( "Hmisc ")
load_R_libraries( "plotrix" )

load_lofs_per_person = function() {
  options(stringsAsFactors = F)
  lofs_per_person = read.table('data/lof/lofs_pop.tsv.gz', header=F)
  names(lofs_per_person) = c('pop', 'sex', 'Hets', 'Homs', 'Rare Hets', 'Rare Homs')
  summary(lofs_per_person)
  return(lofs_per_person)
}

plot_lofs_per_person = function(save_plot = T) {
  lofs_per_person = load_lofs_per_person()
  lof_means = colMeans(subset(lofs_per_person, select=c('Hets', 'Homs', 'Rare Hets', 'Rare Homs')))
  lof_sds = apply(subset(lofs_per_person, select=c('Hets', 'Homs', 'Rare Hets', 'Rare Homs')), 2, sd)
  
  spacing = 0.1
  errcol = 'black'
  if (save_plot) {
    pdf('figures/5a_lofs_per_person.pdf', height=4, width=6)
    par(mar=c(5,4,1,1))
  }
  par(bty='n', las=1)
  barplot(lof_means, space=spacing, col=color_lof, ylim=c(-1, 110), log='', ylab='number of PTV genotypes per individual', xaxt='n')
  errbar(1:4*(1+spacing) - 0.5, lof_means, lof_means-lof_sds, lof_means+lof_sds, add=T, col=errcol, errbar.col=errcol, lwd=2, pch=NA)
  staxlab(side=1, at=1:4*(1+spacing) - 0.5, labels=c('heterozygous', 'homozygous', 'rare hets', 'rare homs'), srt=45)
  text(1:4*(1+spacing) - 0.5, max(lof_means + lof_sds + 5), labels=signif(lof_means, digits=3))
  if (save_plot) dev.off()
}

plot_lofs_per_person_split_by_pop = function() {
  lofs_per_person = load_lofs_per_person()
  text_size = 1.5
  p = ggplot(lofs_per_person) + theme(axis.title.x=element_blank(), axis.title=element_text(size=rel(text_size)), axis.text=element_text(size = rel(text_size)), panel.background = element_blank())
  p + geom_boxplot(aes(factor(pop), `Hets`, fill=tolower(pop))) + scale_fill_manual(values=pop_colors) + guides(fill=F)
  
  pop_lofs = ddply(lofs_per_person, 'pop', function(x) { 
    sds = apply(subset(x, select=-c(pop, sex)), 2, sd)*1.96/sqrt(60706) 
    names(sds) = paste0(names(sds), '_sd');
    c(colMeans(subset(x, select=-c(pop, sex))), sds) })
  p = ggplot(pop_lofs) + 
    theme(axis.title.x=element_blank(), axis.title=element_text(size=rel(text_size)), axis.text=element_text(size = rel(text_size)), panel.background = element_blank()) +
    scale_fill_manual(values=pop_colors) + 
    guides(fill=F)
  p1 = p + geom_bar(aes(x=factor(pop), y = Hets, fill=tolower(pop)), stat='identity', width=1, col='black') + 
    geom_errorbar(aes(x = factor(pop), y = Hets, ymax = Hets + Hets_sd, ymin = Hets - Hets_sd), width=0.25, lwd=2)
  p2 = p + geom_bar(aes(x=factor(pop), y = Homs, fill=tolower(pop)), stat='identity', width=1, col='black') + 
    geom_errorbar(aes(x = factor(pop), y = Homs, ymax = Homs + Homs_sd, ymin = Homs - Homs_sd), width=0.25, lwd=2)
  p3 = p + geom_bar(aes(x=factor(pop), y = `Rare Hets`, fill=tolower(pop)), stat='identity', width=1, col='black') + 
    geom_errorbar(aes(x = factor(pop), y = `Rare Hets`, ymax = `Rare Hets` + `Rare Hets_sd`, ymin = `Rare Hets` - `Rare Hets_sd`), width=0.25, lwd=2)
  p4 = p + geom_bar(aes(x=factor(pop), y = `Rare Homs`, fill=tolower(pop)), stat='identity', width=1, col='black') + 
    geom_errorbar(aes(x = factor(pop), y = `Rare Homs`, ymax = `Rare Homs` + `Rare Homs_sd`, ymin = `Rare Homs` - `Rare Homs_sd`), width=0.25, lwd=2)
  
  p = ggplot(lofs_per_person) + 
    theme(axis.title.x=element_blank(), axis.title=element_text(size=rel(text_size)), axis.text=element_text(size = rel(text_size)), panel.background = element_blank()) +
    scale_fill_manual(values=pop_colors) + scale_color_manual(values=pop_colors) + 
    guides(fill=F, color=F)
  p1 = p + geom_violin(aes(x=factor(pop), y=Hets, fill=tolower(pop), color=tolower(pop)), width=1)
  p2 = p + geom_violin(aes(x=factor(pop), y=Homs, fill=tolower(pop), color=tolower(pop)), width=1)
  p3 = p + geom_violin(aes(x=factor(pop), y=`Rare Hets`, fill=tolower(pop), color=tolower(pop)), width=1)
  p4 = p + geom_violin(aes(x=factor(pop), y=`Rare Homs`, fill=tolower(pop), color=tolower(pop)), width=1)

  multiplot(p1, p2, p3, p4, cols=2)
  #boxplot(as.numeric(lofs_per_person$Hets) ~ lofs_per_person$pop, col='lightgray')
}

# Borrowed from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  load_R_libraries( "grid" )
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}