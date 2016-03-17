source('exac_constants.R')
require(plotrix)

# Load ExAC data
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data(reload=T)
}

# There are NA values in instances where the variant was absent from ESP, i.e. ESP AF = 0,
# so assign AF for these as zero
exac$esp_af_global[is.na(exac$esp_af_global)] = 0.0
exac$esp_af_popmax[is.na(exac$esp_af_popmax)] = 0.0
exac$kg_af_global[is.na(exac$kg_af_global)] = 0.0
exac$kg_af_popmax[is.na(exac$kg_af_popmax)] = 0.0
exac$esp_ac[is.na(exac$esp_ac)] = 0.0
exac$kg_ac[is.na(exac$kg_ac)] = 0.0

# make new filters based on 1kg and ESP
# nonsing_filter is observed allele frequency, but never filtering on a singleton
exac$kg_popmax_nonsing_filter = exac$kg_af_popmax
exac$kg_popmax_nonsing_filter[exac$kg_ac==1] = 0.0
exac$kg_global_nonsing_filter = exac$kg_af_global
exac$kg_global_nonsing_filter[exac$kg_ac==1] = 0.0
exac$esp_popmax_nonsing_filter = exac$esp_af_popmax
exac$esp_popmax_nonsing_filter[exac$esp_ac==1] = 0.0
exac$esp_global_nonsing_filter = exac$esp_af_global
exac$esp_global_nonsing_filter[exac$esp_ac==1] = 0.0
exac$exac_popmax_nonsing_filter = exac$af_popmax
exac$exac_popmax_nonsing_filter[exac$ac_adj==1] = 0.0
exac$exac_global_nonsing_filter = exac$af_global
exac$exac_global_nonsing_filter[exac$ac_adj==1] = 0.0

# Handle the binning of ExAC allele frequencies. This is also done in analysis/properties/af_spectrum.R
# but just in case that hasn't been done yet, recompute it here.
if (!('log_af_bin' %in% colnames(exac))) {
  exac$log_af_bin = floor(log10(exac$af_global))
  exac$log_af_bin[exac$ac_adj == 1] = -6
  exac$log_af_bin[exac$ac_adj > 1 & exac$ac_adj < 11] = -5
  exac$log_af_bin[exac$ac_adj %in% c(11,12)] = -4
  exac$log_af_bin[exac$log_af_bin == 0] = -1 # put the 100% things in the >= 10% category
}

