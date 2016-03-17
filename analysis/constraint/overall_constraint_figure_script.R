# Script to recreate figures 3a and 3b from Lek et al in review

library(binom)
# Load data
source('exac_constants.R')
if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}

# Figure 3A (Z score dist)
constraint_density_plot = function(save_plot=T) {
  if (save_plot) {
    pdf('figures/constraint_z_score_dists.pdf', width=10, height=6)
    par(mar=c(5.1,4.1,4.1,2.1))
  }
  plot(density(constraint$syn_z), xlim=c(-5,10), main='', xlab='Z Score', ylab='', lwd=3, col=alpha(k_syn,0.5), bty='n', yaxt='n', cex.axis=1.2, cex.lab=1.2)
  polygon(density(constraint$syn_z), col=alpha(k_syn,0.5), border='NA')
  lines(density(constraint$mis_z), col=alpha(k_mis,0.5), lwd=3)
  polygon(density(constraint$mis_z), col=alpha(k_mis,0.5), border='NA')
  lines(density(constraint$lof_z), col=alpha(k_lof,0.5), lwd=3)
  polygon(density(constraint$lof_z), col=alpha(k_lof,0.5), border='NA')
  if (save_plot) dev.off()
}

# Figure 3B (pLI + gene lists)
# Importing gene lists of interest
aodl_severe_hi <- read.table('data/constraint/gene_lists/aodl_severe_hi_2015_07_20.txt')
aodl_moderate_hi <- read.table('data/constraint/gene_lists/aodl_moderate_hi_2015_07_20.txt')
aodl_mild_hi <- read.table('data/constraint/gene_lists/aodl_mild_hi_2015_08_04.txt')
genes_essential <- read.table('data/constraint/gene_lists/core_essentials_hart.tsv')
genes_ad <- read.table('data/constraint/gene_lists/all_ad_gencode.tsv')
genes_ar <- read.table('data/constraint/gene_lists/all_ar_gencode.tsv')
genes_olfactory <- read.table('data/constraint/gene_lists/olfactory_receptors_gencode.tsv')

# Finding fraction of gene lists with pLI >= 0.9
all_genes_0.9 <- c(
  nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_severe_hi$V1))/nrow(subset(constraint, gene %in% aodl_severe_hi$V1)),
  nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_moderate_hi$V1))/nrow(subset(constraint, gene %in% aodl_moderate_hi$V1)),
  nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_mild_hi$V1))/nrow(subset(constraint, gene %in% aodl_mild_hi$V1)),
  nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_essential$V1))/nrow(subset(constraint, gene %in% genes_essential$V1)),
  nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_ad$V1))/nrow(subset(constraint, gene %in% genes_ad$V1)),
  nrow(subset(constraint, pLI >= 0.9))/nrow(constraint),
  nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_ar$V1))/nrow(subset(constraint, gene %in% genes_ar$V1)),
  nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_olfactory$V1))/nrow(subset(constraint, gene %in% genes_olfactory$V1))
)

all_gene_names <- c("olfactory", "recessive", "all", "dominant", "essential", "mild HI", "moderate HI", "severe HI")
all_genes_0.9 <- sort(all_genes_0.9)

# Finding upper and lower 95% CI
lower_95_ci_genes <- c(
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_olfactory$V1)),n=nrow(subset(constraint, gene %in% genes_olfactory$V1)),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_ar$V1)),n=nrow(subset(constraint, gene %in% genes_ar$V1)),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9)),n=nrow(constraint),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_ad$V1)),n=nrow(subset(constraint, gene %in% genes_ad$V1)),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_essential$V1)),n=nrow(subset(constraint, gene %in% genes_essential$V1)),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_mild_hi$V1)),n=nrow(subset(constraint, gene %in% aodl_mild_hi$V1)),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_moderate_hi$V1)),n=nrow(subset(constraint, gene %in% aodl_moderate_hi$V1)),method='wilson')$lower,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_severe_hi$V1)),n=nrow(subset(constraint, gene %in% aodl_severe_hi$V1)),method='wilson')$lower
)

upper_95_ci_genes <- c(
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_olfactory$V1)),n=nrow(subset(constraint, gene %in% genes_olfactory$V1)),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_ar$V1)),n=nrow(subset(constraint, gene %in% genes_ar$V1)),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9)),n=nrow(constraint),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_ad$V1)),n=nrow(subset(constraint, gene %in% genes_ad$V1)),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% genes_essential$V1)),n=nrow(subset(constraint, gene %in% genes_essential$V1)),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_mild_hi$V1)),n=nrow(subset(constraint, gene %in% aodl_mild_hi$V1)),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_moderate_hi$V1)),n=nrow(subset(constraint, gene %in% aodl_moderate_hi$V1)),method='wilson')$upper,
  binom.confint(x=nrow(subset(constraint, pLI >= 0.9 & gene %in% aodl_severe_hi$V1)),n=nrow(subset(constraint, gene %in% aodl_severe_hi$V1)),method='wilson')$upper
)

# Plot
constraint_gene_lists = function(save_plot=T) {
  if (save_plot) {
    pdf('figures/constraint_pLI_gene_list_enrichment.pdf', width=16, height=14)
    par(mar=c(6.1,10.1,3.1,2.1))
  }
  text_size = ifelse(save_plot, 2, 1)
  barplot(all_genes_0.9, beside=T, horiz=T, col=c(rep('#558AC2',2), '#969696',rep('#558AC2',5)), border=c(rep('#558AC2',2),'#969696',rep('#558AC2',5)), xaxt='n', xlim=c(0,1))
  col_mids = 1:length(all_genes_0.9)*(1.2)-0.5
  segments(lower_95_ci_genes, col_mids, upper_95_ci_genes, col_mids, lwd=6)
  epsilon = 0.2
  segments(lower_95_ci_genes, col_mids-epsilon, lower_95_ci_genes, col_mids+epsilon, lwd=6)
  segments(upper_95_ci_genes, col_mids-epsilon, upper_95_ci_genes, col_mids+epsilon, lwd=6)
  axis(side=1, at=seq(0,1,0.1), labels=F, lwd=2, lwd.ticks=2, tck=-0.02)
  axis(side=1, at=seq(0,1,0.2), labels=seq(0,1, 0.2), tick=F,line=ifelse(save_plot, 1, 0), cex.axis=text_size)
  abline(v=0, lwd=5)
  axis(side=2, at=c(0.6,1.9,3.1,4.3,5.5,6.7,7.9,9.1), labels=all_gene_names, lwd=0, tick=F, las=1, pos=0.01, cex.axis=text_size)
  if (save_plot) {
  axis(side=2, at=c(0.3,1.6,2.8,4.0,5.2,6.4,7.6,8.8), labels=c(
    paste0("n = ",nrow(subset(constraint, gene %in% genes_olfactory$V1))),
    paste0("n = ",nrow(subset(constraint, gene %in% genes_ar$V1))),
    paste0("n = ",nrow(constraint)),
    paste0("n = ",nrow(subset(constraint, gene %in% genes_ad$V1))),
    paste0("n = ",nrow(subset(constraint, gene %in% genes_essential$V1))),
    paste0("n = ",nrow(subset(constraint, gene %in% aodl_mild_hi$V1))),
    paste0("n = ",nrow(subset(constraint, gene %in% aodl_moderate_hi$V1))),
    paste0("n = ",nrow(subset(constraint, gene %in% aodl_severe_hi$V1)))),
    lwd=0, tick=F, las=1, pos=0.01, cex.axis=text_size)
  }
  mtext('fraction of genes with pLI >= 0.9', side=1, line=ifelse(save_plot, 4, 2.5), cex=text_size*1.2)
  if (save_plot) dev.off()
}
