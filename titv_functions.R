# titv_functions.R

# we'll be using biostrings but also raw characters 
library('Biostrings')

# some definitions
pyr <- c('T', 'C')
pur <- c('A', 'G')
nts <- c(pyr,pur)
weak <- c('A', 'T')
strong <- c('C', 'G')

#########
#
# compute_codon_paths 
#
########

compute_codon_paths <- function() {
	# start with a grid that has all 64*9 codon-to-different-codon paths
	paths <- expand.grid(nts,nts,nts,1:3,c('ti','tv_Cons','tv_WS'))
	colnames(paths)[4:5] <- c("codon_Pos","mut_Type")
	from_Codon <- paste(paths[,1],paths[,2],paths[,3],sep='')
	paths <- cbind(paths, from_Codon)
	
	# now identify the nt to be changed, and what it will be changed to
	# (see assign_mutation for explanation of the 3 types of mutations)
	from_Nt <- apply(paths,1,function(x,y) { substring(x[6],x[4],x[4]) })
	paths <- cbind(paths,from_Nt)
	to_Nt <- apply(paths,1,function(x,y) { assign_mutation(x[7],x[5])})
	
	# we need to do some type-casting here 
	paths <- cbind(paths,to_Nt)
	i <- sapply(paths, is.factor)
	paths[i] <- lapply(paths[i], as.character)
	
	# then compute the new codon resulting from the mutation 
	to_Codon <- apply(paths,1,function(x,y,z) { mutate_codon(x[6],x[4],x[8])})
	paths <- cbind(paths,to_Codon)
	
	# change order for convenience, then get the translations  
	paths <- paths[,c(6,4,7,5,8,9)]
	from_Aa <- apply(paths,1,function(x) { as.character(translate(DNAString(x[1]))) })
	to_Aa <- apply(paths,1,function(x) { as.character(translate(DNAString(x[6]))) })
	paths <- cbind(paths[,1],from_Aa,paths[,2:6],to_Aa)

	# need to convert the translations
	i <- sapply(paths, is.factor)
	paths[i] <- lapply(paths[i], as.character)
	return( paths )
}

#########
#
# compute_aa_paths 
#
# get the set of { from_Aa, to_Aa, mut_type } singlet paths where mut_type is 
#  * "ti" if from_Aa can convert to to_Aa only by a transition mutation 
#  * "tv" if from_Aa can convert to to_Aa only by a transversion mutation 
#  * "ambig" if from_Aa can convert to to_Aa by BOTH ti AND tv paths 
#
# there are 150 ordered singlet paths (75 unordered) 
#
########

compute_aa_paths <- function() {
# get codon paths without unneeded cols; prune out synonymous & nonsense rows
	paths <- compute_codon_paths()[,c('from_Aa','to_Aa','mut_Type')]
	pruned <- paths[paths[,'from_Aa']!=paths[,'to_Aa'] & paths[,'to_Aa'] != '*' & paths[,'from_Aa'] != '*',]

# simplify mut_Type to 2 levels (ti,tv), remove redundancies ==> 158 paths 
	pruned[pruned[,'mut_Type']!='ti','mut_Type'] <- 'tv'
	pruned <- pruned[!duplicated(pruned[,c('from_Aa','to_Aa','mut_Type')]),]

# Now we need to recode the paths that have both ti and tv.  In the cases to remove, there will
# be a path from_Aa --> to_AA path with a ti, and another with a tv.  The "duplicated" function 
# will identify one of these.   These ambig paths will come in forward-reverse pairs, i.e., if 
# G-->R occurs by both ti and tv, the R-->G also occurs by both.  We start by giving every 
# path a string for the unordered (alphabetized) pair of from_Aa, to_Aa, i.e., "AB" is the 
# string for A-->B and B-->A. 

	unordered_Pair <- apply(pruned,1,function(x,y){a<-sort(c(x["from_Aa"],x["to_Aa"]));paste(a[1],a[2],sep='')})
	pruned <- cbind(pruned,unordered_Pair)
	ambig_paths <- pruned[duplicated(pruned[,c("from_Aa","to_Aa")]),"unordered_Pair"]
	ambig_paths <- ambig_paths[!duplicated(ambig_paths)]
	for ( pair in ambig_paths ) {
		pruned[pruned[,"unordered_Pair"]==pair,"mut_Type"] <- "ambig"
	}
	pruned <- pruned[!duplicated(pruned),]
	return( pruned )
}

#########
#
# assign_mutation 
#
#########

# note that each nt has 3 possible mutations, 1 transition and 2 transversions. 
# One of the tv's (tv_Cons) conserves along the WS or weak-strong dimension 
# (i.e., AT vs. GC), and the other (tv_WS) changes along this dimension.  
# For instance T-->C (ti), T-->A (tv_Cons), T-->G (tv_WS). 
	
assign_mutation <- function(x,y) { 
	if ( y=='tv_Cons' ) { 
		if ( x %in% weak ) { return( weak[weak!=x] ) }
		if ( x %in% strong ) { return( strong[strong!=x] )}
	}
	if ( y=='ti' ) { 
		if ( x %in% pur ) { return( pur[pur!=x] ) }
		if ( x %in% pyr ) { return( pyr[pyr!=x] ) }
	}
	if ( y=='tv_WS') { 
		if ( x %in% weak ) { 
			if ( x %in% pyr ) { return( strong[strong!=pyr[pyr!=x]] ) }
			if ( x %in% pur ) { return( strong[strong!=pur[pur!=x]] ) }
		}		
		if ( x %in% strong ) { 
			if ( x %in% pyr ) { return( weak[weak!=pyr[pyr!=x]] ) }
			if ( x %in% pur ) { return( weak[weak!=pur[pur!=x]] ) }
		}		
	}
}


########
#
# mutate_codon
#
#######

mutate_codon <- function( codon, pos, new_nt ) { 
	codon_nts <- strsplit( as.character(codon), NULL)[[1]]
	codon_nts[as.integer(pos)] <- as.character(new_nt)
	return( paste( codon_nts[1],codon_nts[2],codon_nts[3], sep=''))
}


# end of file 