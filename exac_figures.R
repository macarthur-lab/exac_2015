options(stringsAsFactors=FALSE)
#install.packages(c('binom', 'plyr', 'ggplot2', 'data.table', 'reshape', 'plotrix', 'dplyr', 'Hmisc', 'gdata', 'magrittr', 'vioplot'))
library(plyr)
library(dplyr)
library(plotrix)

source('exac_constants.R')
exac = load_exac_data()
use_data = subset(exac, use)
constraint = load_constraint_data()

# Figure 1 - Summary statistics
pdf('figures/figure01.pdf', height=10, width=8.5)
layout(rbind(c(1, 2), c(3, 4), c(5, 5), c(6, 6)), heights=c(1, 1, 0.33, 0.33))
par(cex=1, mar=c(5, 5.5, 1, 2), oma=c(1, 0, 1, 0), pch=20, bty='n')
# 1A - Population breakdown
source('analysis/summary/ancestry_exac_world.R')
plot_panel_size(save_plot=F)
mtext('a', side=3, cex=2, adj = -0.2, line = 0.3)
# 1B - PCA
source('analysis/summary/pca.R')
plot_pca(save_plot=F)
mtext('b', side=3, cex=2, adj = -0.2, line = 0.3)
# 1C - SFS
source('analysis/summary/af_spectrum.R')
af_spectrum = calculate_af_spectrum(exac)
plot_af_spectrum(af_spectrum, save_plot=F)
mtext('c', side=3, cex=2, adj = -0.2, line = 0.3)
# 1D - SFS by category
# Generate counts of observed variants for Figure 1D
if (!file.exists('./data/summary/observed_cutoffs.txt.gz')) {
  source('analysis/summary/prepare_proportion_observed.R')
  observed_cutoffs = generate_observed_df(exac)
  outfh = gzfile('data/summary/observed_cutoffs.txt.gz', 'w')
  write.table(observed_cutoffs, quote=F, sep='\t', row.names=F, file=outfh)
  close(outfh)
  observed_cutoffs = generate_observed_df(use_data)
  outfh = gzfile('data/summary/observed_use_cutoffs.txt.gz', 'w')
  write.table(observed_cutoffs, quote=F, sep='\t', row.names=F, file=outfh)
  close(outfh)
  observed_cutoffs = generate_observed_df(subset(exac, filter == 'PASS'))
  outfh = gzfile('data/summary/observed_pass_cutoffs.txt.gz', 'w')
  write.table(observed_cutoffs, quote=F, sep='\t', row.names=F, file=outfh)
  close(outfh)
}
source('analysis/summary/percent_possible_observed.R')
possible_categories = prepare_all_observed('use')
plot_proportion_observed(possible_categories, leg=T, save_plot=F)
mtext('d', side=3, cex=2, adj = -0.2, line = 0.3)
# 1E - Indels
source('analysis/summary/indel_stats.R')
indels = load_indel_data(exac)
par(mar=c(2, 5.5, 0, 2))
plot_indel_hist(indels)
mtext('e', side=3, cex=2, adj = -0.08, line = 0.3)
par(mar=c(0, 5.5, 0, 2))
plot_indel_af(indels)
mtext('f', side=3, cex=2, adj = -0.08, line = 0)
dev.off()

# Figure 2 - Recurrence
source('analysis/recurrence/plot_curves.R')
source('analysis/recurrence/denovo_recurrence.R')
source('analysis/recurrence/recurrence.R')
source('analysis/summary/category_stats.R')
pdf('figures/figure02.pdf', height=10, width=8.5)
layout(rbind(c(1, 2), c(3, 4), c(5, 5)))
par(cex=1, mar=c(5, 4.1, 2, 2), oma=c(0, 0, 0, 0))
# 2A - Recurrence in de novos
denovos = read_denovo_data()
plot_denovos(denovos, save_plot=F, leg=T)
mtext('a', side=3, cex=2, adj = -0.2, line = 0.3)
# 2B - Saturation curves
mutation_class_samples = load_mutation_class_samples()
plot_mutation_class_samples(mutation_class_samples, save_plot=F, leg=F)
mtext('b', side=3, cex=2, adj = -0.2, line = 0.3)
# 2C - SFS by type
plot_recurrence_sfs(save_plot=F)
mtext('c', side=3, cex=2, adj = -0.2, line = 0.3)
# 2D - Mu vs cross-population between doubletons
doubleton_data = load_doubleton_data()
doubleton_cross_pop(doubleton_data, save_plot=F)
mtext('d', side=3, cex=2, adj = -0.2, line = 0.3)
# 2E - Recurrence-corrected stats on categories
categories = category_comparison(exac)
filtered_categories = filter_categories(categories)
par(mar=c(8.2, 4.1, 2, 2))
plot_category_comparison(filtered_categories, save_plot=F)
mtext('e', side=3, cex=2, adj = -0.08, line = 0.3)
dev.off()

# Figure 3 - Constraint
pdf('figures/figure03.pdf', height=13, width=8.5)
layout(rbind(c(1, 2), c(3, 3), c(4, 5), c(6, 6)))
par(cex=1, oma=c(0, 0, 0, 0))
par(mar=c(5, 2, 2, 2))
source('analysis/constraint/overall_constraint_figure_script.R')
# 3A - Constraint
constraint_density_plot(save_plot=F)
mtext('a', side=3, cex=2, adj = -0.05, line = 0.3)
# 3B - Gene lists
par(mar=c(2, 4, 2, 2))
constraint_gene_lists(save_plot=F)
mtext('b', side=3, cex=2, adj = -0.2, line = 0.3)
# 3C - Expression
source('analysis/constraint/expression_constraint.R')
gtex_data = prepare_gtex_data()
par(mar=c(2, 4, 2, 0))
plot_expression_constraint(gtex_data[[1]], gtex_data[[2]], save_plot=F)
mtext('c', side=3, cex=2, adj = -0.08, line = 0.3)
# 3D - eQTLs
source('analysis/constraint/eqtl_constraint.R')
par(mar=c(5, 4, 2, 2))
plot_eqtl_tissue(save_plot=F)
mtext('d', side=3, cex=2, adj = -0.2, line = 0.3)
# 3E - GWAS
source('analysis/constraint/GWAS_set_constraint.R')
gwas_data = prepare_gwas_data()
plot_gwas_figure(gwas_data, save_plot=F)
mtext('e', side=3, cex=2, adj = -0.2, line = 0.3)
# 3F - Categories with constraint - warning: this takes a long time
source('analysis/summary/category_stats.R')
source('analysis/constraint/category_stats_with_constraint.R')
exac_mutation = get_exac_with_mutation_data()
categories = category_comparison(exac_mutation)
filtered_categories = filter_categories(categories, min_n_category = 10000)
categories_detail = prepare_category_constraint_data(exac_mutation)
par(mar=c(10, 4, 1, 2))
plot_category_stats_constraint(filtered_categories, categories_detail, save_plot=F)
mtext('f', side=3, cex=2, adj = -0.08, line = 0.3)
dev.off()

# Figure 4 - Mendelians
source('analysis/mendelians/mendelian_figures.R')
pdf('figures/figure04.pdf',width=10,height=8.5)
par(mfrow = c(2, 2))
par(cex=1, mar=c(5, 5, 2, 1), oma=c(0, 0, 0, 0))
# 4A - ESP, ExAC AF filter comparison
source('analysis/mendelians/variants_remaining_after_filters.R')
variants_remaining = get_variants_remaining_after_filters()
plot_variants_remaining(variants_remaining, save_plot = F)
mtext('a', side=3, cex=2, adj = -0.2, line = 0.3)
# 4B - AF Estimate comparison
source('analysis/mendelians/exac_esp_afs.R')
esp_af_data = load_esp_af_data()
plot_exac_esp_af(esp_af_data, save_plot = F)
mtext('b', side=3, cex=2, adj = -0.2, line = 0.3)
# 4C - Pathogenic variants
source('analysis/mendelians/pathogenic_variants.R')
pathogenic_variants = load_pathogenic_variants()
plot_pathogenic_variants(pathogenic_variants)
mtext('c', side=3, cex=2, adj = -0.2, line = 0.3)
# 4D - Curation
source('analysis/mendelians/curated_associations.R')
curation_figure()
mtext('d', side=3, cex=2, adj = -0.2, line = 0.3)
dev.off()

# Figure 5 - LoFs
source('analysis/lof/load_downsampled_to_hist.R')
source('analysis/lof/pop_af_hist.R')
source('analysis/lof/lofs_per_person.R')
source('analysis/lof/plot_lof_discovery.R')
pdf('figures/figure05.pdf', height=8, width=8.5)
par(mfrow = c(2, 2))
par(cex=1, mar=c(5, 4, 2, 2), oma=c(0, 0, 0, 0))
# 5A - Number of LoFs per person
plot_lofs_per_person(save_plot=F)
mtext('a', side=3, cex=2, adj = -0.12, line = 0.3)
# 5B - LoFs across populations
downsampled_data = load_downsampled_to_hist(pops, count_alleles=T)/3000
plot_pop_af_hist(downsampled_data, save_plot=F)
mtext('b', side=3, cex=2, adj = -0.12, line = 0.3)
# 5C - Number of LoF genes per population
lof_discovery_data = load_lof_discovery_curves(hom=F)
plot_lof_discovery_curves(lof_discovery_data, hom=F, save_plot=F, last_point=10000)
mtext('c', side=3, cex=2, adj = -0.12, line = 0.3)
# 5D - Number of hom LoF genes per population
lof_discovery_data = load_lof_discovery_curves(hom=T)
plot_lof_discovery_curves(lof_discovery_data, hom=T, save_plot=F, last_point=10000)
mtext('d', side=3, cex=2, adj = -0.12, line = 0.3)
dev.off()

# Supplemental Figures/Data

# QC Metrics
## Table of QC measurement
## "Use" criteria
# source('analysis/supplement/use_analysis.R')
## Multi-allelics
# source('analysis/supplement/genome_plotting.R')

# Summary stats
## PCA 1 vs 2
## Number of variants per person by class and ancestry

# Supplementary Table X - Proportion variants observed by category
source('analysis/summary/percent_possible_observed.R')
covered_possible_observed = subset(load_all_observed('use'), cutoff == 0.8 & cov == 10)
library(dplyr)
library(magrittr)
cpgs_possible = covered_possible_observed %>% 
  subset(type == 'CpG transition' & freq.poss > 200) %>% 
  arrange(desc(freq.obs/freq.poss)) %>%
  select(consequence, `CpGs possible` = freq.poss, `CpGs observed` = freq.obs) %>%
  mutate(`% CpGs` = percent(`CpGs observed`/`CpGs possible`, 3))
transitions_possible = covered_possible_observed %>% 
  subset(type == 'non-CpG transition') %>% 
  select(consequence, `non-CpG transitions possible` = freq.poss, `non-CpG transitions observed` = freq.obs) %>%
  mutate(`% non-CpG transition` = percent(`non-CpG transitions observed`/`non-CpG transitions possible`, 3))
transversions_possible = covered_possible_observed %>% 
  subset(type == 'transversion') %>% 
  select(consequence, `transversions possible` = freq.poss, `transversions observed` = freq.obs) %>%
  mutate(`% transversion` = percent(`transversions observed`/`transversions possible`, 3))
cpgs_possible %>% 
  left_join(transitions_possible) %>% 
  left_join(transversions_possible) %>%
  mutate(consequence = format_vep_category(consequence)) %>%
  write.table(quote=F, row.names=F, sep='\t')

# Recurrence section
# SFigure R1
pdf('figures/s_figure_recurrence.pdf', height=10, width=8.5)
layout(rbind(c(1, 2), c(3, 3), c(4, 5)))
par(cex=1, mar=c(5, 4, 2, 2), oma=c(0, 0, 0, 0))
source('analysis/recurrence/plot_curves.R')
mutation_class_samples = load_mutation_class_samples()
plot_ti_tv(mutation_class_samples, save_plot=F)
mtext('a', side=3, cex=2, adj = -0.2, line = 0.3)
# SFigure R2
source('analysis/recurrence/recurrence.R')
# doubleton_data = load_doubleton_data()
doubleton_cross_pop(doubleton_data, distance=T, save_plot=F)
mtext('b', side=3, cex=2, adj = -0.2, line = 0.3)
# SFigure R3
par(mar=c(8, 4, 2, 2))
source('analysis/supplement/recurrence_correction.R')
source('analysis/summary/category_stats.R')
categories_uncorr = category_comparison(exac, correction = F)
filtered_categories_uncorr = filter_categories(categories_uncorr)
plot_category_comparison(filtered_categories_uncorr, correction = F, save_plot=F)
mtext('c', side=3, cex=2, adj = -0.08, line = 0.3)
# SFigure R4
par(mar=c(7, 4, 2, 2))
plot_mu_vs_singleton(save_plot=F)
mtext('d', side=3, cex=2, adj = -0.2, line = 0.3)
# SFigure R5
titv_categories = get_titv_singleton_by_category()
plot_titv_singleton_by_category(titv_categories, save_plot=F)
mtext('e', side=3, cex=2, adj = -0.2, line = 0.3)
dev.off()

# LoF section

# Frequency spectrum of LoFs in high pLI genes, broken down by population
source('analysis/lof/load_downsampled_to_hist.R')
source('analysis/lof/pop_af_hist.R')
afs_constr = load_downsampled_to_hist(pops, count_alleles=T, constrained_only=T)/3000
pdf('figures/s_figure_ptv.pdf', width=6, height=4)
par(mar=c(3,4,1,1))
plot_pop_af_hist_db(afs_constr, save_plot=F)
dev.off()
