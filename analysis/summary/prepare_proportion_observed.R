library(data.table)
library(plyr)

source('exac_constants.R')
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}
if (!("use_data" %in% ls(globalenv()))) {
  use_data = subset(exac, use)
}

load_all_possible_variants = function() {
  if (file.exists('data/synthetic.vep.cov.table')) {
    all_possible = fread('data/synthetic.vep.cov.table', data.table=F)
    names(all_possible) = tolower(colnames(all_possible))
    names(all_possible)[1] = 'chrom'
    all_possible$cpg = (all_possible$ref == 'C' & all_possible$alt == 'T' & substr(all_possible$context, 2, 3) == 'CG') | (all_possible$ref == 'G' & all_possible$alt == 'A' & substr(all_possible$context, 1, 2) == 'CG')
    all_possible$transition = (all_possible$ref == 'A' & all_possible$alt == 'G') | (all_possible$ref == 'G' & all_possible$alt == 'A') | (all_possible$ref == 'C' & all_possible$alt == 'T') | (all_possible$ref == 'T' & all_possible$alt == 'C')
    return(all_possible)
  } else {
    print('Error - you do not have synthetic.header.vep.cov.table in your data/ directory.')
    print('You can download it from ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/manuscript_data/all_possible_variants/synthetic.vep.cov.table.gz.')
    print('However, be advised that this is a 2.5G download, and needs to be un-gzipped (to 36G) to run this command.')
    print('This code is largely for reproducibility/documentation purposes, and you can find the much smaller summary statistics in data/summary/all_possible_cutoffs.txt.gz which can be loaded with load_all_possible_variants().')
  }
}

generate_possible_df = function(all_possible) {
  # Calculate possible variants at all coverages
  cutoffs = (0:20)/20
  table_data_detail_30 = ldply(cutoffs, function(x) { 
    print(paste('30', x))
    a = plyr::count(subset(all_possible, coverage_30 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  table_data_detail_20 = ldply(cutoffs, function(x) { 
    print(paste('20', x))
    a = plyr::count(subset(all_possible, coverage_20 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  table_data_detail_15 = ldply(cutoffs, function(x) { 
    print(paste('15', x))
    a = plyr::count(subset(all_possible, coverage_15 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  table_data_detail_10 = ldply(cutoffs, function(x) { 
    print(paste('10', x))
    a = plyr::count(subset(all_possible, coverage_10 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  table_data_detail_10$cov = 10
  table_data_detail_15$cov = 15
  table_data_detail_20$cov = 20
  table_data_detail_30$cov = 30
  table_data_detail = rbind(table_data_detail_10, table_data_detail_15, table_data_detail_20, table_data_detail_30)
  outfh = gzfile('data/all_possible_cutoffs.txt.gz', 'w')
  write.table(table_data_detail, quote=F, sep='\t', row.names=F, file=outfh)
  close(outfh)
}

# Observed
generate_observed_df = function(exac) {
  cutoffs = (0:20)/20
  observed_30 = ldply(cutoffs, function(x) { 
    print(paste('30', x))
    a = plyr::count(subset(exac, coverage_30 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  observed_20 = ldply(cutoffs, function(x) { 
    print(paste('20', x))
    a = plyr::count(subset(exac, coverage_20 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  observed_15 = ldply(cutoffs, function(x) { 
    print(paste('15', x))
    a = plyr::count(subset(exac, coverage_15 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  observed_10 = ldply(cutoffs, function(x) { 
    print(paste('10', x))
    a = plyr::count(subset(exac, coverage_10 >= x, select=c(consequence, cpg, transition)))
    a$cutoff = x
    return(a)
  })
  observed_10$cov = 10
  observed_15$cov = 15
  observed_20$cov = 20
  observed_30$cov = 30
  observed_cutoffs = rbind(observed_10, observed_15, observed_20, observed_30)
  return(observed_cutoffs)
}

