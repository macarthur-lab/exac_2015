source('pop_af_hist.R')

load_downsampled_to_hist <- function(pops, breaks=c(-Inf,0.001,0.01,0.05,Inf), s=T, d=F, t=F, constrained_only=F, canonical_only=F, return_seperate=F, count_alleles=F, complete_only=F, filter=NA) {
  num_bins <- length(breaks) - 1 + s + d + t
  tot_hist <- matrix(0, nrow=num_bins, ncol=length(pops))

  down_df <- get(load('data/lof/lofs_downsampled_an6000.RData'))
  down_df[,paste("an",pops,sep="_")] <- down_df$ans
  
  # if no filter is provided, but a filter is desired, get that filter
  if (constrained_only | complete_only | canonical_only) {
    if (is.na(filter)) {
      if (constrained_only) {
        if (!("constraint" %in% ls(globalenv()))) {
          constraint = load_constraint_data()
        }
        constrained_genes = subset(constraint, lof_phi_cut == '(0.9,1]')$gene
        filter = down_df$symbol %in% constrained_genes
      } else if (complete_only) {
        source('analysis/lof/partial.R')
        filter = sapply(down_df$feature, is_not_partial)
      } else if (canonical_only) {
        filter = down_df$canonical == 'YES'
      } 
    }
  }

  if (!is.na(filter)) {
    down_df = down_df[filter,] 
  }
  
  # get histogram for current downsampled dataframe
  h = pop_af_hist(down_df, pops, sing=s, doub=d, trip=t, breaks=breaks, count_alleles=count_alleles)
  colnames(h) = pops
  if (return_seperate) {
    tot_hist = lapply(pops, function(pop) { h[,pop] })
    names(tot_hist) = pops
  } else {
    tot_hist <- tot_hist + h 
  }
  
  if (!return_seperate) {
    tot_hist <- data.frame(tot_hist)
    colnames(tot_hist) <- pops
  }
  return(tot_hist)
}
