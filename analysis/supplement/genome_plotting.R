# source("http://bioconductor.org/biocLite.R")
# biocLite("ggbio")

library(ggbio)
data(hg19Ideogram, package = "biovizBase")
library(GenomicRanges)

plot_snp_density = function(exac_data, resolution=1000) {
  intervals=(0:(250000000/resolution))*resolution
  # Generating bins of genomic positions
  exac_data$pos_bin = cut(exac_data$pos, intervals)
  exac_data$bin_start = as.numeric( sub("\\((.+),.*", "\\1", exac_data$pos_bin))
  exac_data$bin_end = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", exac_data$pos_bin))
  data = subset(exac_data, select=c(chrom, bin_start, bin_end))
  
  # Counting variants in these bins
  counts = plyr::count(data, vars = c('chrom', 'bin_start', 'bin_end'))
  
  # Creating GRanges objects from your data
  density <- with(counts, GRanges(paste0("chr", chrom), IRanges(bin_start, bin_end), strand="+", freq))
  # Consistent chromosome lengths/which chromosomes to use
  seqlengths(density) <- seqlengths(hg19Ideogram)[names(seqlengths(density))]
  density <- keepSeqlevels(trim(density), paste0("chr", c(1:22, "X", "Y")))
  #autoplot(density, layout = "karyogram")
  
  autoplot(density, layout = "karyogram", aes(color=freq, fill=freq, size=freq)) + # cheating a bit since it increases size of box based on frequency
    scale_colour_gradient(low = "#132B43", high = "#56B1F7", space = "Lab", na.value = "#031B33", guide = "colourbar")
}

get_multiallelic_density = function(exac_data, min_allelic_level=3, resolution=1000) {
  intervals=(0:(250000000/resolution))*resolution
  # Only get multi-allelic sites
  allelic_state = plyr::count(subset(exac_data, select=c(chrom, pos)))
  multi_allelics = subset(allelic_state, freq > min_allelic_level, select=c(chrom, pos))
  # Create bins and count
  multi_allelics$pos_bin = cut(multi_allelics$pos, intervals)
  multi_allelics$bin_start = as.numeric( sub("\\((.+),.*", "\\1", multi_allelics$pos_bin))
  multi_allelics$bin_end = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", multi_allelics$pos_bin))
  multiallelic_counts = plyr::count(multi_allelics, vars = c('chrom', 'bin_start', 'bin_end'))
  #head(multiallelic_counts[order(multiallelic_counts$freq, decreasing = T),], 10)
  return(multiallelic_counts)
}
  
plot_multiallelic_density = function(multiallelic_counts) {
  # GRanges and Plot
  multiallelic_density <- with(multiallelic_counts, GRanges(paste0("chr", chrom), IRanges(bin_start, bin_end), strand="+", freq))
  seqlengths(multiallelic_density) <- seqlengths(hg19Ideogram)[names(seqlengths(multiallelic_density))]
  multiallelic_density <- keepSeqlevels(trim(multiallelic_density), paste0("chr", c(1:22, "X", "Y")))
  
  # Overlay chromosomes with rectangles
  autoplot(multiallelic_density, layout = "karyogram", aes(color=freq, fill=freq, size=freq)) + # cheating a bit since it increases size of box based on frequency
    scale_colour_gradient(low = "#132B43", high = "#56B1F7", space = "Lab", na.value = "#031B33", guide = "colourbar")
  
  # Plot on top of chromosomes with cytobands
#   hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
#   ggplot(hg19) + layout_karyogram(cytoband=TRUE, ylim=c(0,10)) + layout_karyogram(multiallelic_density, aes(color=freq, x=start, y=freq), geom = "line", ylim = c(11, 21))
}
