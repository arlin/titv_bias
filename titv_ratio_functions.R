# compute mean effect of ti or tv for paths defined by source & dest amino acid 
#

source("titv_functions.R")
source("aa_wt_matrix_functions.R")
source("codon_usage_functions.R")

# some defaults for testing 
cu_species <- human
aa_matrix <- exs

# this is intended for a symmetric aa replacement matrix

# this is the main function for computing expected effects on titv ratio from 
#  * the genetic code (via compute_codon_paths)diag.present=\1,\r\tis.log.odds=FALSE\r\tis.distance=FALSE
#  * codon usage (via metadata object from codon_usage_metadata.R)
#  * aa replacement cost (via metadata object from aa_matrices_metadata.R)
# 
# If you want to compute effects of codon use alone, use aa_matrix = uniform_aa. 
# Likewise, for exchangeability effects alone, use cu_species = uniform_cu. 
# 
# To use a custom codon usage or aa matrix, you must add metadata according to the 
# pattern in codon_usage_metadata.R and aa_matrices_metadata.R. 
# 
# below is a function to iterate use of this 

# the function for using a symmetric aa_matrix (unordered paths) is just a thin shell 
# around the function for getting possibly ordered paths

weight_codon_paths <- function( cu_species, aa_matrix, paths ) 
{
	# compute paths if the user hasn't supplied the output of this command: 
	if ( missing(paths) ) { 
		paths <- compute_codon_paths()
	}
	
	# get paths, remove synonymous & nonsense, add the 2-char code for each replacement 	
	paths <- paths[paths[,'from_Aa']!=paths[,'to_Aa'] & paths[,'to_Aa'] != '*' & paths[,'from_Aa'] != '*',]
	if ( aa_matrix$is.square ) { 
		aa_paths <- apply(paths,1,function(x,y){paste(x["from_Aa"],x["to_Aa"],sep='')})
	}
	else {
		aa_paths <- apply(paths,1,function(x,y){a<-sort(c(x["from_Aa"],x["to_Aa"]));paste(a[1],a[2],sep='')})
	}
	paths <- cbind(paths[,c(1,7,5)],aa_paths)
	colnames(paths) <- c("from_Codon", "to_Codon", "mut_Type", "Repl_code")
	
	# simplify to only one type of tv
	paths[paths[,'mut_Type']!='ti','mut_Type'] <- 'tv'

	# get the aa weight matrix in 2 cols (2-char replacement code, then value)
	aa_repl_wts <- flattened_aa_matrix( aa_matrix )
	
	# in case this is a distance matrix rather than a similarity matrix, we'll 
	# flip the number line around so that max is min, and min is max
	if ( aa_matrix$is.distance ) { 
		aa_repl_wts <- rescale_distance_to_similarity( aa_repl_wts )
	}

	# in case these are log-of-odds values, we'll exponentiate so that ratio of 
	# sums isn't misleading due to negative values. 
	if ( aa_matrix$is.log.odds ) { 
		aa_repl_wts <- exp( aa_repl_wts )
	}

	# now merge, 
	aa_weighted <- merge(paths, aa_repl_wts, by.x=4,by.y="row.names")
	
	# now get codon usage 
	cu <- read_codon_usage( cu_species ) 
	
	# now merge again 
	both <- merge(aa_weighted,cu,by.x="from_Codon",by.y="row.names")
	combined <- both[,5] * both[,6]
	both <- cbind(both,combined)
	return( both ) 
} 

titv_ratio_by_codon_paths <- function( cu_species, aa_matrix, paths ) 
{
	out <- weight_codon_paths( cu_species, aa_matrix, paths ) 
	ti <- sum(out[out[,'mut_Type']=="ti","combined"])
	tv <- sum(out[out[,'mut_Type']=="tv","combined"])
	return( ti/tv ) 
}

# cus - a list of codon usage objects as defined in codon_usage_metadata.R 
# aa_matrix - an aa matrix object as defined in aa_matrices_metadata.R 
# paths - table of codon paths generated by compute_codon_paths()

iterate_matrix_with_cus <- function( aa_matrix, cus, paths ) {
	# compute paths if the user hasn't supplied the output of this command: 
	if ( missing(paths) ) { 
		paths <- compute_codon_paths()
	}
	result <- matrix(,nrow=1,ncol=3)
	colnames(result) <- c("AA_matrix", "CodonUsage", "TiTvRatio") 
	for (i in 1:length( cus )) { 
		cu_species <- cus[[i]]
		out <- titv_ratio_by_codon_paths( cu_species, aa_matrix, paths)
		result <- rbind(result,c(aa_matrix$name, cu_species$name, out))
	}
	return(result[-1,])
}	


# this computes the ratio by counting each aa replacement only once, as either a 
# replacement that takes place only by transitions, or only by transversions.  
# Replacements that can take place by both transitions and transversions are ignored. 
# 
titv_ratio_by_aa_paths <- function( aa_matrix, distance=FALSE) 
{

	# compute all singlet codon paths with ti, tv or ambig
	paths <- compute_aa_paths()
	unord_paths <- paths[!duplicated(paths[,4]),]
	
	# get the matrix in 2 cols, unordered pair and value
	tmp <- flattened_symmetric_aa_matrix( aa_matrix )

	# in case this is a distance matrix rather than a similarity matrix, we'll 
	# flip the number line around so that max is min, and min is max
	if ( distance ) { 
		tmp[,2] <- rescale_distance_to_similarity( as.numeric(as.character(tmp[,2])) )
	}
	# do the merge and the calculations
	tmp <- merge(tmp,unordered_paths,by.x=1,by.y=4)[,c(1,2,5)]
	dat <- data.frame("weight"=as.numeric(as.character(tmp[,2])),"group"=factor(tmp[,3]))
#	aggregate(weight ~ group, data=dat, FUN = function(x) c(M=mean(x), SD=sd(x)))
	sum(dat[dat[,2]=="ti",1])/sum(dat[dat[,2]=="tv",1])
}

# end of file 
