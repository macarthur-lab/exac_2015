
#### Begin preparation for Figure 4B

load_esp_af_data = function() {
  # If 'misc_data/exac_esp.table.gz' is present (locally), read it in along with the 
  # ESP public release data, and re-compute the allele frequency distributions
  # that are used for Figure 4B, and then write them to disk. If the file is not present (public release)
  # then just read them in from disk
  if ('exac_esp.table.gz' %in% list.files('./misc_data/')) {
    # subset of ExAC that came from ESP
    exac_esp = read.table('misc_data/exac_esp.table.gz',header=TRUE,sep='\t',comment.char='',quote='')
    colnames(exac_esp) = tolower(colnames(exac_esp))
    exac_esp$pos_id = paste(exac_esp$chrom, formatC(exac_esp$pos,width=9,flag='0'), exac_esp$ref, exac_esp$alt, sep='_')
    match_indices = match(exac$pos_id, exac_esp$pos_id)
    exac$exac_esp_ac_nfe = exac_esp$ac_nfe[match_indices]
    exac$exac_esp_an_nfe = exac_esp$an_nfe[match_indices]
    exac$af_nfe_no_esp = (exac$ac_nfe - exac$exac_esp_ac_nfe) / (exac$an_nfe - exac$exac_esp_an_nfe)
    exac$af_nfe_no_esp[(exac$an_nfe - exac$exac_esp_an_nfe)==0] = 0.0
    exac$af_nfe_no_esp[is.na(exac$af_nfe_no_esp)] = 0.0
    
    # ESP public release
    esp = read.table('data/mendelians/esp_af.tsv.gz',header=TRUE,sep='\t',comment.char='',quote='')
    esp$pos_id = paste(esp$chrom, formatC(esp$pos,width=9,flag='0'), esp$ref, esp$alt, sep='_')
    esp$not_in_exac = !(esp$pos_id %in% exac$pos_id)
    esp$indel = nchar(esp$ref) != nchar(esp$alt)
    esp$af_nfe = esp$ac_eur / esp$an_eur
    esp$af_nfe[esp$an_eur==0] = 0.0
    match_indices = match(exac$pos_id, esp$pos_id)
    exac$esp_ac_nfe = esp$ac_eur[match_indices]
    exac$esp_an_nfe = esp$an_eur[match_indices]
    exac$esp_af_nfe = exac$esp_ac_nfe / exac$esp_an_nfe
    exac$esp_af_nfe[exac$esp_an_nfe==0] = 0.0
    exac$esp_af_nfe[is.na(exac$esp_af_nfe)] = 0.0
    
    exac$af_nfe = exac$ac_nfe / exac$an_nfe
    exac$af_nfe[is.na(exac$af_nfe)] = 0.0
    
    esp_afs = c(.01,.001,mean(esp$af_nfe[esp$ac_eur==1]))
    
    af1 = sort(log10(c(exac$af_nfe_no_esp[!exac$indel & exac$esp_af_nfe > .009 & exac$esp_af_nfe < .011], 
                       rep(0,sum(esp$af_nfe > .009 & esp$af_nfe < .011 & esp$not_in_exac & !esp$indel)))))
    af1[af1==-Inf] = -5.55
    af01 = sort(log10(c(exac$af_nfe_no_esp[!exac$indel & exac$esp_af_nfe > .0009 & exac$esp_af_nfe < .0011],
                        rep(0,sum(esp$af_nfe > .0009 & esp$af_nfe < .0011 & esp$not_in_exac & !esp$indel)))))
    af01[af01==-Inf] = -5.60
    af001 = sort(log10(c(exac$af_nfe_no_esp[!exac$indel & exac$esp_af_nfe > 0 & exac$esp_ac == 1],
                         rep(0,sum(esp$af_nfe > 0 & esp$ac_eur==1 & esp$not_in_exac & !esp$indel)))))
    af001[af001==-Inf] = -5.65
    write.table(data.frame(af1),'data/mendelians/fig4b_af1.tsv',row.names=F,col.names=T,quote=F)
    write.table(data.frame(af01),'data/mendelians/fig4b_af01.tsv',row.names=F,col.names=T,quote=F)
    write.table(data.frame(af001),'data/mendelians/fig4b_af001.tsv',row.names=F,col.names=T,quote=F)
  } else {
    af1 = read.table('data/mendelians/fig4b_af1.tsv',header=T)$af1
    af01 = read.table('data/mendelians/fig4b_af01.tsv',header=T)$af01
    af001 = read.table('data/mendelians/fig4b_af001.tsv',header=T)$af001
    esp_afs = c(.01,.001,1/8600)
  }
  return(list(af1, af01, af001, esp_afs))
}

#### End preparation for Figure 4B

plot_exac_esp_af = function(esp_af_data, save_plot=F) {
  # color constants for Figure 4B
  color_1 = '#111111' # color for AF 1% histogram
  color_01 = '#156615' # color for AF 0.1% histogram
  color_001 =  '#2566FF' # color for singleton histogram
  esp_colors = c(color_1,color_01,color_001)
  
  af1 = esp_af_data[[1]]
  af01 = esp_af_data[[2]]
  af001 = esp_af_data[[3]]
  esp_afs = esp_af_data[[4]]
  # Figure 5B - AF_NFE distribution in ExAC of variants at AF_NFE of 1%, .1% and AC=1 in ESP
  
  # Decide where to break this histogram and what y axis limit to set
  h_last = hist(af001,breaks=100,plot=FALSE)
  broken_hist_value = h_last$density[1]
  ymax = 5 #+ broken_hist_value - floor(broken_hist_value)
  
  # Labels on the x axis (log10 allele frequency)
  x_axis_labels = c(-5:-1, log10(3*10^(-5:-1))) 
  
  if (save_plot) {
    pdf('figures/nfe_af_by_esp_eur_af.pdf',width=4,height=4)
    par(mar=c(6,6,1.5,1))
  }
  hist(af1,xlim=c(-5.8,-1),breaks=100,
       prob=TRUE,axes=FALSE,ylim=c(0,ymax),ylab='',xlab='',col=alpha(color_1,.7),border=NA,main='',yaxs='i')
  par(new=TRUE)
  hist(af01,xlim=c(-5.8,-1),breaks=100,
       prob=TRUE,axes=FALSE,ylim=c(0,ymax),xlab='',ylab='',col=alpha(color_01,.7),border=NA,main='',yaxs='i')
  par(new=TRUE)
  h_last = hist(af001,xlim=c(-5.8,-1),breaks=100,
                prob=TRUE,axes=FALSE,ylim=c(0,ymax),xlab='',ylab='',col=alpha(color_001,.7),border=NA,main='',yaxs='i',plot=TRUE)
  broken_hist_value = h_last$density[1]
  # Because the breaks are in units of .05, it is the case that for all i, h_last$density[i] * .05 == h_last$counts[i] / sum(h_last$counts)
  # So to plot the y axis in terms of percentages instead of probability density, multiply by .05
  axis(side=2,at=0:5,labels=percent(c(0:4,broken_hist_value)*.05),lwd=2,lwd.ticks=1,las=1)
  axis.break(axis=2,breakpos=4.5,style='slash')
  abline(h=0,lwd=2)
  mtext(side=2,text='percentage of variants',line=2.9,font=1,cex=1.1)
  mtext(side=1,text='ExAC European AF excluding ESP',line=3.7,font=1,cex=1.1)
  segments(x0=log10(esp_afs),y0=rep(0,4),y1=rep(3.9,4),col=esp_colors,lty=3,lwd=2)
  text(x=log10(esp_afs),y=rep(4,4),pos=3,labels=c('1%','0.1%','AC=1'),col=esp_colors,font=1,cex=1)
  text(x=mean(log10(esp_afs)),y=4.5,pos=3,labels='ESP European AF',font=1,cex=1.1)
  segments(x0=min(log10(esp_afs))-.2,x1=max(log10(esp_afs))+.1,y0=4.5,lwd=1,col='#000000')
  axis(side=1,at=c(-5.6,x_axis_labels),labels=NA,lwd=0,lwd.ticks=1,cex.axis=.7,las=2,tck=-.01,line=0,hadj=1)
  axis(side=1,at=c(-5.6,x_axis_labels),labels=c('0%',percent(10^x_axis_labels)),tick=F,cex.axis=1,las=2,line=-.5,hadj=1)
  axis.break(1,-5.4,style='slash')
  if (save_plot) dev.off()
}