options(stringsAsFactors=FALSE)
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}

if (!("seen" %in% colnames(exac))) {
  dbsnp = read.table('data/summary/dbsnp.table.sorted.gz',sep='\t',header=FALSE)
  colnames(dbsnp) = c('chrom','pos','ref','alt')
  dbsnp$pos_id = paste(dbsnp$chrom, formatC(dbsnp$pos,width=9,flag='0'), dbsnp$ref, dbsnp$alt, sep='_')
  exac$in_dbsnp = exac$pos_id %in% dbsnp$pos_id
  exac$in_kg_esp = exac$kg_af_popmax > 0 | exac$esp_af_popmax > 0
  exac$dbsnp_other = exac$in_dbsnp & !exac$in_kg_esp
  exac$seen = exac$in_dbsnp | exac$in_kg_esp
  
  # compute AF bins. 
  exac$log_af_bin = floor(log10(exac$af_global))
  exac$log_af_bin[exac$ac_adj == 1] = -6
  exac$log_af_bin[exac$ac_adj > 1 & exac$ac_adj < 11] = -5
  exac$log_af_bin[exac$ac_adj %in% c(11,12)] = -4
  exac$log_af_bin[exac$log_af_bin == 0] = -1 # put the 100% things in the >= 10% category
}

# Move to statistics_in_text.R?
# overall proportion novel (70.8%)
# sum(!exac$seen, na.rm=T) / nrow(exac)
# within < .01% global AF subset (78.4%)
# sum(!exac$seen & exac$af_global < .0001) / sum(exac$af_global < .0001)
# proportion < .01% global AF
# sum(exac$af_global < .0001) / dim(exac)[1]
# table(exac[,c("in_dbsnp","in_kg_esp")])
# exac[head(which(exac$in_kg_esp & !exac$in_dbsnp)),]
# table(exac$indel[exac$in_kg_esp & !exac$in_dbsnp]) # they're about 83% indels
# mean(exac$kg_ac[exac$in_kg_esp & !exac$in_dbsnp] == 0) # they're about 83% from ESP, absent from 1kg
# table(exac$indel[exac$in_kg_esp & !exac$in_dbsnp & exac$esp_ac == 0])
# table(exac$indel[exac$in_kg_esp & !exac$in_dbsnp & exac$esp_ac > 0])
# sum(exac$in_dbsnp[exac$esp_ac > 0 & exac$indel])
# sum(!exac$in_dbsnp[exac$esp_ac > 0 & !exac$indel])
# sum(exac$in_dbsnp[exac$esp_ac > 0 & !exac$indel])
# 
# sum(!exac$in_dbsnp[exac$kg_ac > 0 & !exac$indel])
# sum(exac$in_dbsnp[exac$kg_ac > 0 & !exac$indel])
# 
# sum(exac$in_kg_esp)
# sum(exac$in_dbsnp)
# sum(exac$dbsnp_other)
# test = head(dbsnp)
# paste(test$chrom, formatC(test$pos,width=9,flag='0'), test$ref, test$alt, sep='_')
# head(exac$pos_id)

calculate_af_spectrum = function(exac) {
  # aggregate and sqldf are both inexplicably and intolerably slow for this task, so I am just doing a loop
  af_spectrum = as.data.frame(table(exac$log_af_bin[exac$use]))
  colnames(af_spectrum) = c('log_af_bin','n')
  af_spectrum$log_af_bin = -6:-1 # otherwise table() makes it a factor
  af_spectrum$seen = 0
  af_spectrum$in_kg_esp = 0
  for (i in 1:dim(af_spectrum)[1]) {
    rows = exac$use & exac$log_af_bin == af_spectrum$log_af_bin[i]
    af_spectrum$seen[i] = sum(exac$seen[rows], na.rm=T)
    af_spectrum$in_kg_esp[i] = sum(exac$in_kg_esp[rows], na.rm=T)
  }
  af_spectrum$proportion = af_spectrum$n / sum(af_spectrum$n)
  af_spectrum$proportion_seenbar = af_spectrum$proportion * af_spectrum$seen / af_spectrum$n
  af_spectrum$proportion_kgespbar = af_spectrum$proportion * af_spectrum$in_kg_esp / af_spectrum$n
  return(af_spectrum)
}

plot_af_spectrum = function(af_spectrum, save_plot=F) {
  novel_color = '#003EFF'
  seen_color = '#666666'
  kg_esp_color = '#000000'
  bar_lwd = 25
  if (save_plot) {
    cairo_pdf('figures/af_spectrum_overall.pdf',width=4,height=4)
    par(mar=c(5,5,1,1))
  }
  plot(NA,NA,xlim=c(-6.5,-.5),ylim=c(0,.6),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
  points(af_spectrum$log_af_bin, af_spectrum$proportion, col=novel_color, type='h', lwd=bar_lwd, lend=1)
  points(af_spectrum$log_af_bin, af_spectrum$proportion_seenbar, col=seen_color, type='h', lwd=bar_lwd, lend=1)
  points(af_spectrum$log_af_bin, af_spectrum$proportion_kgespbar, col=kg_esp_color, type='h', lwd=bar_lwd, lend=1)
  abline(h=0)
  #axis(side=1, at=-6:-1, labels=c('singleton','2-10 alleles','0.01%','0.1%','1%','\u226510%'),srt=45,line=-.5,lwd=0,lwd.ticks=0,cex.axis=.5)
  text(x=-6:-1, y=par("usr")[3], xpd=T, labels=c('singleton','2-10 alleles','0.01%','0.1%','1%','10%+'),srt=45,adj=c(1.1,1.1),font=2,cex=.8)
  axis(side=2, at=(0:6)/10, labels=percent((0:6)/10), lwd=0, lwd.ticks=1, las=1)
  mtext(side=1, line=3.2, text='ExAC global allele frequency')
  mtext(side=2, line=4, text='proportion of ExAC alleles')
  legend('topright',c('novel','other dbSNP','in 1kG or ESP'),col=c(novel_color,seen_color,kg_esp_color),pch=15, bty='n')
  if (save_plot) dev.off()
}

plot_af_spectrum2 = function(af_spectrum, save_plot=F) {
  # 1st panel: histogram bars for # of variants in each global AF bin
  # 2nd panel: curves for syn, mis, lof as proportion of self, by AF bin
  # 3rd panel: proportion absent from ESP + 1kG, by AF bin
  
  # consider whether 1st and 3rd could be combined by stacked histogram
  lof_dist = table(exac$log_af_bin[exac$use & exac$lof_use])/sum(exac$use & exac$lof_use)
  mis_dist = table(exac$log_af_bin[exac$use & !is.na(exac$category) & exac$category=='missense_variant'])/sum(exac$use & !is.na(exac$category) & exac$category=='missense_variant')
  syn_dist = table(exac$log_af_bin[exac$use & !is.na(exac$category) & exac$category=='synonymous_variant'])/sum(exac$use & !is.na(exac$category) & exac$category=='synonymous_variant')
  
  fclass_spectrum = data.frame(bin=-6:-1, lof=as.numeric(lof_dist), mis=as.numeric(mis_dist), syn=as.numeric(syn_dist))
  if (save_plot) {
    cairo_pdf('figures/af_spectrum_byclass.pdf',width=8,height=4)
    par(mar=c(4,6,2,2))
  }
  plot(NA,NA,xlim=c(-6.5,-.5),ylim=c(0,.7),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='')
  for (fclass in c('lof','mis','syn')) {
    points(fclass_spectrum$bin, fclass_spectrum[,fclass], col=get(paste('k',fclass,sep='_')), type='l', lwd=5)
  }
  abline(h=0)
  axis(side=1, at=-6:-1, labels=c('singleton','2-10 alleles','0.01%','0.1%','1%','\u226510%'),lwd=0,lwd.ticks=1,cex.axis=.9)
  axis(side=2, at=(0:7)/10, labels=percent((0:7)/10), lwd=0, lwd.ticks=.5, las=1)
  mtext(side=1, line=2, text='ExAC global allele frequency', font=2)
  mtext(side=2, line=3.5, text='Proportion of alleles in functional class', font=2)
  if (save_plot) dev.off()
}
# by population
# q = sqldf("
# select   popmax,
#          avg(seen) proportion_seen,
#          avg(in_kg_esp) proportion_kgesp
# from     exac
# where    ac_adj = 1
# and      use
# group by 1
# order by 1
# ;")

# # while sqldf is far more elegant in code, it is extremely slow, as is any use of aggregate
# # that invokes the subset function, because subset() creates a subsetted copy of the ENTIRE data frame
# # (all cols) in memory. the fastest way i've found to do this in R is sadly the ugliest in code:
# sing_bypop = aggregate(exac$seen[exac$use & exac$ac_adj==1] ~ exac$popmax[exac$use & exac$ac_adj==1], FUN="mean")
# sing_bypop2 = aggregate(exac$in_kg_esp[exac$use & exac$ac_adj==1] ~ exac$popmax[exac$use & exac$ac_adj==1], FUN="mean")
# colnames(sing_bypop) = c('pop','proportion_seen')
# colnames(sing_bypop2) = c('pop','proportion_kgesp')
# sing_bypop$proportion_kgesp = sing_bypop2$proportion_kgesp[match(sing_bypop$pop, sing_bypop2$pop)]
# sing_bypop$xval = 1:6
# sing_bypop$acol = get(paste('k',tolower(sing_bypop$pop),sep='_'))
# 
# nonsing_bypop = aggregate(exac$seen[exac$use & exac$ac_adj > 1] ~ exac$popmax[exac$use & exac$ac_adj > 1], FUN="mean")
# nonsing_bypop2 = aggregate(exac$in_kg_esp[exac$use & exac$ac_adj > 1] ~ exac$popmax[exac$use & exac$ac_adj > 1], FUN="mean")
# colnames(nonsing_bypop) = c('pop','proportion_seen')
# colnames(nonsing_bypop2) = c('pop','proportion_kgesp')
# nonsing_bypop$proportion_kgesp = nonsing_bypop2$proportion_kgesp[match(nonsing_bypop$pop, nonsing_bypop2$pop)]
# 
# 
# 
# # just double checking this AFR thing b/c it was surprising
# sum(exac$ac_adj==1 & exac$use & exac$popmax=='AFR',na.rm=TRUE)
# sum(exac$ac_adj==1 & exac$use & exac$popmax=='AFR' & exac$seen,na.rm=TRUE)


