library(binom)
library(gdata)
library(reshape)
library(Hmisc)
source('exac_constants.R')
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}

read_denovo_data = function() {
  denovos = read.delim('data/recurrence/denovos.vep.table.gz')
  names(denovos) = tolower(names(denovos))
  denovos$pos_id = paste(denovos$chrom, formatC(denovos$pos,width=9,flag='0'), denovos$ref, denovos$alt, sep='_')
  denovos$in_exac = denovos$pos_id %in% exac$pos_id
  denovos$indel = nchar(denovos$ref) != nchar(denovos$alt)
  denovos$transition = (denovos$ref == 'A' & denovos$alt == 'G') | (denovos$ref == 'G' & denovos$alt == 'A') | (denovos$ref == 'C' & denovos$alt == 'T') | (denovos$ref == 'T' & denovos$alt == 'C')
  denovos$transition[denovos$indel] = NA
  denovos$cpg = (denovos$ref == 'C' & denovos$alt == 'T' & substr(denovos$context, 2, 3) == 'CG') | (denovos$ref == 'G' & denovos$alt == 'A' & substr(denovos$context, 1, 2) == 'CG')
  denovos$type = ifelse(denovos$transition, 'non-CpG transition', 'transversion')
  denovos$type[denovos$cpg] = 'CpG transition'
  denovos$category = format_vep_category(denovos$consequence)

  # subset(denovos, cpg & consequence == 'synonymous_variant' & !in_exac) # Check coverage for these
  # table(denovos$in_exac, denovos$cpg)
  return(denovos)
}

plot_denovos = function(denovos, save_plot=T, leg=F) {
  denovos_filt = subset(denovos, !indel & consequence %in% c('synonymous_variant', 'missense_variant', 'stop_gained'))

  data = ddply(denovos_filt, c('type', 'category'), function(x) {
    data.frame(perc_in_exac=sum(x$in_exac)/nrow(x), subset(binom.confint(x=sum(x$in_exac), n=nrow(x), method='wilson'), select=c(lower, upper)))
  })
  
  variant_types = c('transversion', 'non-CpG transition', 'CpG transition')
  data$category <- reorder.factor(data$category, new.order=c('nonsense', 'missense', 'synonymous'))
  data$type <- reorder.factor(data$type, new.order=variant_types)
  data = arrange(data, category, type)
  plot_data = t(cast(data, category ~ type, value='perc_in_exac'))
  uppers = t(cast(data, category ~ type, value='upper'))
  lowers = t(cast(data, category ~ type, value='lower'))
  font_size = 1
  
  errcol = 'black'
  spacing = c(0, 1)
  if (save_plot) pdf('figures/2a_denovos_in_exac.pdf', width=6, height=4)
  a = barplot(plot_data, beside=T, col=c(color_tv, color_ti, color_cpg), ylab='proportion found in ExAC', space=spacing, ylim=c(0,1), yaxt='n')
  axis(side=2, at=(0:5)/5, labels=percent((0:5)/5), lwd=0, lwd.ticks=1, las=1)
  errbar(a, plot_data, lowers, uppers, add=T, col=errcol, errbar.col=errcol, lwd=2, pch=NA)
  if (leg) legend('topleft', variant_types, col=c(color_tv, color_ti, color_cpg), text.col=c(color_tv, color_ti, color_cpg), pch=15, bty='n')
  if (save_plot) dev.off()
#   p = ggplot(data) + geom_bar(aes(x=category, fill=Type, y=perc_in_exac), stat='identity', position='dodge', col='black') +
#     ylab('Proportion found in ExAC') + theme(axis.text=element_text(size = rel(font_size)), axis.title.x=element_blank(), axis.title=element_text(size = rel(font_size)), legend.text=element_text(size = rel(font_size)), legend.title=element_text(size = rel(font_size)), panel.background = element_blank(), legend.position='bottom')
#   if (save_plot) {
#     pdf('figures/2a_denovos_in_exac.pdf', width=6, height=4)
#     print(p)
#     dev.off()
#   } else {
#     return(p)
#   }
}
