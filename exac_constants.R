# ExAC constants file
# Use by source('exac_constants.R')

load_bioc_libraries = function( libName )
{
	if ( !require( libName, character.only = TRUE ) )
	{
		source( "http://www.bioconductor.org/biocLite.R" )
		biocLite( libName, suppressUpdates=TRUE )
		
		if ( !require( libName, character.only = TRUE ) )
		{
			stop( "Couldn't install", libName, "from Bioconductor. Please install manually." )
		}
    }
}

load_R_libraries = function( libName )
{
	if ( !require( libName, character.only=TRUE ) )
	{
		install.packages( libName )

		if ( !require( libName, character.only=TRUE ) )
		{
			stop( "Couldn't install", libName, "from CRAN packages. Please install manually." )
		}

	}
}

devtools_load =  function( package, github )
{
	# install package
	if (!require( package, character.only=TRUE ) )
	{
		# devtools
		load_R_libraries( 'devtools' )
		
		devtools::install_github( github )
		if ( !require( package, character.only=TRUE ) )
		{
			stop( "Couldn't install", package, "using devtools package. Please install manually." )
		}
	}
}

options(stringsAsFactors=FALSE)

load_R_libraries( 'plyr' )
load_R_libraries( 'dplyr' )
load_R_libraries( 'magrittr' )

calling_interval_length = 59684884
num_samples = 60706

# The colors!
color_amr = k_amr = '#ED1E24'
color_eur = k_eur = '#6AA5CD'
color_afr = k_afr = '#941494'
color_sas = k_sas = '#FF9912'
color_eas = k_eas = '#108C44'
color_oth = k_oth = '#ABB9B9'
color_mde = k_mde = '#000080'

color_nfe = k_nfe = color_eur
color_fin = k_fin = '#002F6C'

color_syn = k_syn = '#AAAAAA'
color_mis = k_mis = '#FF6103'
color_lof = k_lof = '#9D1309'

color_cpg = '#2E9FFE'
color_ti = '#458B00'
color_tv = '#EA4444'

lof_like = c('frameshift_variant','essential_splice','stop_gained')
mis_like = c('missense_variant','inframe_indel','stop_lost',
             'mature_miRNA_variant','start_lost')
syn_like = c('synonymous_variant','3_prime_UTR_variant','5_prime_UTR_variant',
             'extended_splice','stop_retained_variant','non_coding_transcript_exon_variant',
             'intron_variant','intergenic_variant','regulatory_region_variant')

format_vep_category = function(category_list) {
  return(category_list %>%
    gsub("_"," ", .) %>%
    gsub('stop gained', 'nonsense', .) %>%
    gsub("inframe indel", "in-frame indel", .) %>%
    gsub("initiator codon", "start lost", .) %>%
    gsub(" variant", "", .) %>%
    gsub("transcript exon", "transcript", .) %>%
    gsub(" prime ","'", .) %>%
    gsub("probably damaging", "prob damaging", .) %>%
    gsub("possibly damaging", "poss damaging", .))
}
variant_types = c('transversion', 'non-CpG transition', 'CpG transition')
variant_type_colors = c(color_tv, color_ti, color_cpg)

# example usage: alpha(k_lof,.5) gives you a 50% transparent LoF maroon
alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

pops <- c('afr', 'amr', 'eas', 'fin', 'nfe', 'sas')
pop_colors = c('afr' = color_afr,
               'amr' = color_amr,
               'eas' = color_eas,
               'fin' = color_fin,
               'nfe' = color_nfe, 
               'oth' = color_oth,
               'sas' = color_sas,
               'mde' = color_mde,
               'uniform' = 'pink',
               'consanguineous' = 'pink',
               'sas_non_consang' = 'orange')
pop_names = c('afr' = 'African',
             'amr' = 'Latino',
             'eas' = 'East Asian',
             'fin' = 'Finnish',
             'nfe' = 'European',
             'oth' = 'Other',
             'sas' = 'South Asian',
             'mde' = 'Middle Eastern',
             'uniform' = 'Uniform',
             'sas_non_consang' = 'South Asian (F < 0.05)',
             'consanguineous' = 'South Asian (F > 0.05)')


## Population sizes
# Populations and percent ancestries
# World in millions, all others actual
# From Laramie
populationAncestries <- data.frame(
  row.names=c("East Asian", "South Asian", "European", "Middle Eastern", "African", "Latino", "Oceanic", "DiverseOther", "AfricanEuropeanAdmixed"),
  world=c(1932, 2085, 1145, 410,  1022, 529, 38, NA, NA),
  kgenomes=c(523, 494, 514, 0, 691, 355, 0, NA, NA),
  esp=c(0, 0, 4298, 0, 2217,0,0, NA, NA),
  exac=c(4327, 8256, 36677, 0, 5203, 5789, 0,454, NA),  # Numbers from Monkol on 3/4/15
  pcgcGwas=c(5219, 0, 116766, 0, 0, 0, 0, NA, NA),   # PGC Published 4 (SCZ, BIP, MDD, ADHD)
  ptsd=c(106, 0, 8393, 0, NA, 829, 0, 1295, 9845)
)

# Loading data
data_url = 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/manuscript_data/'
open_file_exac = function(fname) {
  con = gzcon(url(paste0(data_url, fname)))
  dat = readLines(con)
  return(read.delim(textConnection(dat), header=T))
}

load_exac_data = function(type='', private=FALSE, canonical=FALSE, reload=FALSE) { # This takes ~10 minutes (download, then ~5 mins of processing)
  time_start <- proc.time()
  transcripts = ifelse(type != '', paste0('.', type), '')
  canon = ifelse(canonical, '.canonical', '')
  data_format = ifelse(private, '.full', '')
  fname_unzipped = paste0('ExAC.r0.3.1.sites.vep', data_format, transcripts, canon, '.table')
  fname = paste0(fname_unzipped, '.gz')
#   if (reload | (!require('data.table') & !file.exists(fname)) | (require('data.table') & !file.exists(fname_unzipped))) {
  if (reload | !file.exists(fname)) {
    print(paste('Step (0/6): Downloading', fname, '...'))
    download.file(paste0(data_url, fname), fname)
  } else {
    print('Using locally loaded exac file. Run exac = load_exac_data(reload=T) to update to newest file.')
  }
  print('Step (1/6): Loading data...')
#   if (require('data.table')) {
#     if (!file.exists(fname_unzipped)) {
#       print('Unfortunately, fread must operate on unzipped files :-(')
#       print('Unzipping (will save for future loads)')
#       print(paste('If space is tight, you may want to delete', fname_unzipped))
#       system(paste('gzip -d -c', fname, '>', fname_unzipped))
#     }
#     exac = suppressWarnings(fread(fname_unzipped, data.table=F))
#   } else {
#     print('data.table not found. Using normal read.delim... install.packages("data.table")!')
    exac = read.delim(fname, header=T)
#   }
  print('Now processing data (2/6): Determining sequence contexts...')
  colnames(exac) = tolower(colnames(exac))
  exac$pos_id = paste(exac$chrom, formatC(exac$pos,width=9,flag='0'), exac$ref, exac$alt, sep='_')
  exac$indel = nchar(exac$ref) != nchar(exac$alt)
  exac$bases_inserted = nchar(exac$alt) - nchar(exac$ref)
  exac$transition = (exac$ref == 'A' & exac$alt == 'G') | (exac$ref == 'G' & exac$alt == 'A') | (exac$ref == 'C' & exac$alt == 'T') | (exac$ref == 'T' & exac$alt == 'C')
  exac$transition[exac$indel] = NA
  exac$cpg = (exac$ref == 'C' & exac$alt == 'T' & substr(exac$context, 2, 3) == 'CG') | (exac$ref == 'G' & exac$alt == 'A' & substr(exac$context, 1, 2) == 'CG')
  exac$cpg[exac$indel] = NA
  exac$alt_is_ancestral = exac$alt == exac$ancestral
  exac$alt_is_ancestral[exac$ref != exac$ancestral & exac$alt != exac$ancestral] = NA
  
  print('(3/6): Computing allele frequencies...')
  exac$af_popmax = exac$ac_popmax / exac$an_popmax
  exac$af_popmax[exac$an_popmax == 0 | is.na(exac$an_popmax) | is.na(exac$ac_popmax)] = 0.0
  exac$af_global = exac$ac_adj / exac$an_adj
  exac$af_global[exac$an_adj == 0 | is.na(exac$an_adj) | is.na(exac$ac_adj)] = 0.0
  exac$singleton = exac$ac_adj == 1
  exac$maf_global = pmin(exac$af_global, 1-exac$af_global) # minor allele frequency
  exac$daf_global = ifelse(exac$alt_is_ancestral, 1-exac$af_global, exac$af_global)
  exac$daf_popmax = ifelse(exac$alt_is_ancestral, 1-exac$af_popmax, exac$af_popmax)
  
  print('(4/6): Calculating bad regions of the genome...')
  resolution = 1000
  number_of_regions = 10
  intervals=(0:(250000000/resolution))*resolution
  exac$pos_bin = cut(exac$pos, intervals)
  exac$bin_start = as.numeric( sub("\\((.+),.*", "\\1", exac$pos_bin))
  allelic_state = plyr::count(subset(exac, select=c(chrom, pos)))
  multiallelics = subset(allelic_state, freq > 3, select=c(chrom, pos))
  multiallelics$pos_bin = cut(multiallelics$pos, intervals)
  multiallelics$bin_start = as.numeric( sub("\\((.+),.*", "\\1", multiallelics$pos_bin))
  multiallelic_counts = plyr::count(multiallelics, vars = c('chrom', 'bin_start'))
  bad_sectors = subset(head(multiallelic_counts[order(multiallelic_counts$freq, decreasing = T),], number_of_regions), select=c(chrom, bin_start))
  bad_sectors$sector = paste(bad_sectors$chrom, bad_sectors$bin_start, sep='_')
  
  print('(5/6): Determining what to _use_...')
  exac$sector = paste(exac$chrom, exac$bin_start, sep='_')
  exac$bad = exac$sector %in% bad_sectors$sector
  exac$use = exac$an_adj > .80*max(exac$an_adj, na.rm=TRUE) & exac$ac_adj > 0 & exac$filter=='PASS' & !exac$bad
  exac$lof_use = !is.na(exac$lof) & exac$lof == 'HC' & is.na(exac$lof_flags)
  
  print('(6/6): Parsing SIFT and PolyPhen scores, and separating functional categories...')
  # for the idea behind these named group regexes, see http://stackoverflow.com/a/2969666/3806692
  sift_regex = "([a-z]*)\\(([0-9\\.]*)\\)"
  exac$sift_word = sub(sift_regex, "\\1", exac$sift)
  exac$sift_word[exac$sift_word == ''] = NA
  exac$sift_score = as.numeric(sub(sift_regex, "\\2", exac$sift)) 
  polyphen_regex = "([a-z_]*)\\(([0-9\\.]*)\\)"
  exac$polyphen_word = sub(polyphen_regex, "\\1", exac$polyphen) # parse words
  exac$polyphen_word[exac$polyphen_word == ''] = NA
  exac$polyphen_score = as.numeric(sub(polyphen_regex, "\\2", exac$polyphen)) # parse scores
  exac$category = exac$consequence
  exac$category[exac$consequence=='inframe_insertion'] = 'inframe_indel'
  exac$category[exac$consequence=='inframe_deletion'] = 'inframe_indel'
  exac$category[exac$consequence=='splice_donor_variant'] = 'essential_splice'
  exac$category[exac$consequence=='splice_acceptor_variant'] = 'essential_splice'
  exac$category[exac$consequence=='splice_region_variant'] = 'extended_splice'
  
  print('Done! Here you go. Took:')
  print(proc.time() - time_start)
  return(exac)
}

load_all_possible_variants = function() {
  return(read.table('data/summary/all_possible_cutoffs.txt.gz', header=T))
}

load_constraint_data = function(reload=FALSE) {
  fname = 'forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz'
  if (reload | !file.exists(fname)) {
    print(paste('Downloading', fname, '...'))
    download.file(paste0(data_url, fname), fname)
  } else {
    print('Using locally loaded constraint file. Run constraint = load_constraint_data(reload=T) to update to newest file.')
  }
  constraint = read.delim(fname, header=T)
  constraint$feature = sapply(constraint$transcript, function(x) { strsplit(x, '.', fixed=TRUE)[[1]][[1]] })
  
  boundaries = c(0, 0.25, 0.75, 1)
  # boundaries = c(0, 0.25, 0.5, 0.75, 1)
  constraint$lof_cut = cut(constraint$lof_z, quantile(constraint$lof_z, probs=boundaries))
  constraint$mis_cut = cut(constraint$mis_z, quantile(constraint$mis_z, probs=boundaries))
  constraint$syn_cut = cut(constraint$syn_z, quantile(constraint$syn_z, probs=boundaries))
  
  constraint$lof_phi_cut = cut(constraint$pLI, breaks=c(0, 0.1, 0.9, 1))
  
  # Confirmed to Kaitlin's criteria (3230 transcripts)
  constraint$lof_constrained = constraint$pLI > 0.9
  constraint$mis_constrained = constraint$mis_z > 3.09 & constraint$syn_z < 3.09
  return(constraint)
}

load_cadd_exac = function(reload=F) {
  if (reload | !file.exists('data/ExAC.r0.3.tsv.gz')) {
    print('Downloading CADD data...')
    ret = system('bgzip')
    gzip_type = ifelse(ret <= 1, 'b', '')
    system(paste0("curl http://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3.tsv.gz | zcat | sed '1d' | ", gzip_type, "gzip -c > data/ExAC.r0.3.tsv.gz"))
  } else {
    print('Using locally loaded CADD file. Run exac = load_cadd_exac(reload=T) to update to newest file.')
  }
  print('Loading CADD data...')
  cadd_data = read.table('data/ExAC.r0.3.tsv.gz', header=T, comment.char="")
  colnames(cadd_data) = tolower(colnames(cadd_data))
  colnames(cadd_data)[1] = 'chrom'
  cadd_data$pos_id = paste(cadd_data$chrom, formatC(cadd_data$pos,width=9,flag='0'), cadd_data$ref, cadd_data$alt, sep='_')
  cadd_data = subset(cadd_data, select=c(pos_id, rawscore, phred))
  print('Done!')
  return(cadd_data)
}

load_gene_lists = function() {
#   list_metadata = read.table(textConnection("
# filename|display|xval
# universe|All genes|-2
# all_ar|Autosomal recessive|-3
# all_ad|Autosomal dominant|-4
# core_essentials_hart|Essential in culture|-5
# haploinsufficient|Haploinsufficient|-6
# mgi_essential|Essential in mice|-7
# fda_approved_drug_targets|Drug targets|-8
# grep_or|Olfactory receptors|-9
# gwascatalog|Near GWAS hits|-10
# homozygous_lof_tolerant_twohit|Homozygous LoF tolerant|-11
# "),sep='|',header=TRUE)
  
  if (!file.exists('gene_lists')) {
    system('git clone git@github.com:macarthur-lab/gene_lists.git')
  } else {
    print('Using locally available gene_lists... Run git pull to update.')
  }
  gene_lists = ldply(Sys.glob('gene_lists/lists/*'), function(x) {
#     gene_list = read.table(paste('gene_lists/lists/', x, '.tsv',sep=''))
    gene_list = read.table(x)
    names(gene_list) = 'symbol'
    gene_list$list = strsplit(strsplit(x, '/', fixed=T)[[1]][3], '.', fixed=T)[[1]][1]
    gene_list
  })
  
  all_genes = dcast(gene_lists, symbol ~ list, fun.aggregate = function(x) { length(x) >= 1 })
}


load_mutation_data = function(exac) {
  mutations = load_raw_mutation_data()
  print('Calculating metrics for synonymous SNPs...')
  synonymous_snps = subset(exac, !indel & use & nchar(context) == 3 & category == 'synonymous_variant')
  bymutation = ddply(synonymous_snps, c('context', 'ref', 'alt'), function(x) {
    data.frame(n=nrow(x), 
               singletons=sum(x$singleton), 
               doubletons=sum(x$ac_adj == 2),
               tripletons=sum(x$ac_adj == 3),
               quad=sum(x$ac_adj == 4),
               quint=sum(x$ac_adj == 5),
               ac_gt_five=sum(x$ac_adj >= 5))
  })
  mutation_data = merge(bymutation, mutations)
  mutation_data$color = ifelse(mutation_data$transition, ifelse(mutation_data$cpg, color_cpg, color_ti), color_tv)
  
  print('Done loading mutation data!')
  return(mutation_data)
}

get_k_run <- function(run_col) {
  v = unlist(sapply(run_col, strsplit, ":"))
  v = v[seq(2, length(v), 2)]
  v = as.numeric(v)
  return(v)
}

# constraint = load_constraint_data()
# constrained_transcripts = subset(constraint, lof_phi_cut == '(0.9,1]')$feature
is_phi_constr <- function(feature) {
  feats <- strsplit(feature, split=",")[[1]]
  which_feats = feats %in% constrained_transcripts
  return(any(which_feats))
}

load_raw_mutation_data = function() {
  mutations = read.table('data/recurrence/fordist_1KG_mutation_rate_table.txt', header=T)
  
  mutations$ref = substr(mutations$from, 2, 2)
  mutations$alt = substr(mutations$to, 2, 2)
  mutations$context = mutations$from
  mutations$cpg = (mutations$ref == 'C' & mutations$alt == 'T' & substr(mutations$context, 2, 3) == 'CG') | (mutations$ref == 'G' & mutations$alt == 'A' & substr(mutations$context, 1, 2) == 'CG')
  mutations$transition = (mutations$ref == 'A' & mutations$alt == 'G') | (mutations$ref == 'G' & mutations$alt == 'A') | (mutations$ref == 'C' & mutations$alt == 'T') | (mutations$ref == 'T' & mutations$alt == 'C')
  return(mutations)
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}