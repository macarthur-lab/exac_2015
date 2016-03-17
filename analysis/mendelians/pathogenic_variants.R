#### Begin preparation for Figure 4C

load_pathogenic_variants = function() {
  # If HGMD dataset is available (locally) then re-calculate the stats for Figure 4C
  # and write them to disk. If not (public release) then just re-loads the stats from disk
  if ('hgmd_201402_normalized.tsv' %in% list.files('./misc_data/')) {
    
    # Load ClinVar dataset and join to ExAC
    if (!("clinvar" %in% ls(globalenv()))) {
      # We are using the 2015-07 freeze of ClinVar.
      # This is a local copy of this file: https://github.com/macarthur-lab/clinvar/blob/b218de6b26c5008ab0a18a372ba165c01059eb69/output/clinvar.tsv
      # And info on the generation of this file can be found in: https://github.com/macarthur-lab/clinvar
      clinvar = read.table('data/mendelians/clinvar_2015_07.tsv',sep='\t',header=TRUE,quote='',comment.char='')
    }
    
    clinvar$pos_id = paste(clinvar$chrom, formatC(clinvar$pos,width=9,flag='0'), clinvar$ref, clinvar$alt, sep='_')
    # consider only ClinVar variants that are asserted as pathogenic, are non-conflicted, and refer
    # to the ALT allele
    exac$clinvar_use = exac$pos_id %in% clinvar$pos_id[clinvar$pathogenic==1 & clinvar$conflicted==0 & clinvar$mut=='ALT']
    
    # load HGMD dataset or, if not available (because we don't have permission to publicly release it), 
    # then just load summary stats
    if (!("hgmd" %in% ls(globalenv()))) {
      hgmd = read.table('misc_data/hgmd_201402_normalized.tsv',sep='\t',header=TRUE,quote='',comment.char='')
      hgmd$pos_id = paste(hgmd$chrom, formatC(hgmd$pos,width=9,flag='0'), hgmd$ref, hgmd$alt, sep='_')
    }
    require(sqldf)
    mendelian = sqldf("
                      select   pos_id, max(symbol) symbol from (
                      select   pos_id, symbol
                      from     clinvar
                      where    pathogenic=1 and conflicted=0
                      union
                      select   pos_id, symbol
                      from     hgmd
                      ) group by 1 order by 1;")
    mendelian$clinvar = mendelian$pos_id %in% clinvar$pos_id[clinvar$pathogenic==1 & clinvar$conflicted==0]
    mendelian$hgmd = mendelian$pos_id %in% hgmd$pos_id
    mendelian$log_af_bin = exac$log_af_bin[match(mendelian$pos_id, exac$pos_id)]
    mendelian$log_af_bin[is.na(mendelian$log_af_bin)] = -7 # default is zeroton
    mendelian$inheritance = as.character(NA)
    mendelian$inheritance[mendelian$symbol %in% dominant_list] = 'autosomal dominant'
    mendelian$inheritance[mendelian$symbol %in% recessive_list] = 'autosomal recessive'
    mendelian$use = as.integer((mendelian$pos_id %in% exac$pos_id[exac$use]) | (!(mendelian$pos_id %in% exac$pos_id)))
    mendelian$het_only = as.integer(mendelian$pos_id %in% exac$pos_id[exac$use & exac$ac_hom == 0])
    
    fig4c_data = sqldf("
                       select   log_af_bin,
                       inheritance,
                       case when inheritance = 'autosomal dominant' then -1 else 1 end x_offset,
                       sum(use) total,
                       sum(het_only) het_only
                       from     mendelian
                       where    use
                       and      inheritance is not null
                       group by 1, 2
                       order by 1, 2
                       ;")
    
    write.table(fig4c_data,'data/mendelians/fig4c_data.tsv',sep='\t',row.names=F,col.names=T,quote=F)
  } else {
    fig4c_data = read.table('data/mendelians/fig4c_data.tsv',sep='\t',header=T)
  }
  return(fig4c_data)
}

#### End preparation for Figure 4C

plot_pathogenic_variants = function(fig4c_data, save_plot = F) {
  # color constants for Figure 4C
  dom_color = '#D89516'
  rec_color = '#169AAE'
  block_color = '#FFFFFF80' 
  
  # Anne's figure 4C
  hist_width = .8
  fig4c_data$xleft = fig4c_data$log_af_bin - (hist_width/2) * (fig4c_data$inheritance=='autosomal dominant')
  fig4c_data$xright = fig4c_data$log_af_bin + (hist_width/2) * (fig4c_data$inheritance=='autosomal recessive')
  fig4c_data$color = ''
  fig4c_data$color[fig4c_data$inheritance=='autosomal dominant'] = dom_color
  fig4c_data$color[fig4c_data$inheritance=='autosomal recessive'] = rec_color
  
  if (save_plot) {
    cairo_pdf('figures/4c.pdf',width=6,height=6)
    par(mar=c(6,6,1.5,1))
  }
  plot(NA, NA, xlim=c(-7.5, -0.5), ylim=c(0, 7000), xaxs='i', yaxs='i', xlab='', ylab='', main='', axes=FALSE)
  rect(xleft=fig4c_data$xleft[3:14], xright=fig4c_data$xright[3:14], ybottom=0, ytop=fig4c_data$total[3:14], col=alpha(fig4c_data$color[3:14],.7), border=NA)
  abline(h=0,lwd=2)
  mtext(side=1, at=-7:-1, text=c('absent','singleton','AC=2-10','0.01%', '0.1%', '1%', '10%+'), cex=.65)
  axis(side=2, at=(0:6)*1000, labels=formatC(c((0:4)*1000, fig4c_data$total[2], fig4c_data$total[1]),big.mark=',',format='fg'), las=2, cex.axis=1, lwd=2)
  # plot the zerotons manually
  rect(xleft=-7 - hist_width/2, xright=-7, ybottom=0, ytop=6000, col=alpha(dom_color,.7), border=NA)
  rect(xleft=-7, xright=-7 + hist_width/2, ybottom=0, ytop=5000, col=alpha(rec_color,.7), border=NA)
  axis.break(2, 1000*4.5, style='gap')
  axis.break(2, 1000*5.5, style='gap')
  legend('topright',c('autosomal dominant','autosomal recessive'),col=alpha(c(dom_color, rec_color),.7),text.col=alpha(c(dom_color, rec_color),.7),pch=15,bty='n',cex=1,text.font=1)
  mtext(side=1, line=2, text='ExAC global allele frequency', cex=1.1)
  mtext(side=2, line=3.8, text='ClinVar and HGMD pathogenic variants', font=1, cex=1.1)
  if (save_plot) dev.off()
}