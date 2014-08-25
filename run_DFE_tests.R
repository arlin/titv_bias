# DFE_U_tests.R -- repeat Mann-Whitney U test and report stats from DFE studies
# 
# this generates the data for Table 2.  Just do this:
#
# > source( "run_DFE_tests.R" )
# 
# this generates output files "DFE_U_tests.tab", which you can easily turn
# into a formatted table using the spreadsheet "DFE_results_formatted_table.xslx".  
# Just paste the results into the "raw" sheet and they will appear formatted in 
# the other sheet.  Fisher's method is easier to implement in the spreadsheet.
#
# This also runs the meta-statistical tests and puts these in a separate table.   

# get the full version of the Mann-Whitney from conditional inference package
library(coin)
source("DFE_tests/DFE_functions.R")

# get the data from all the studies (in separate file)
source( "DFE_tests/DFE_data.R" )
all_dfes <- list( firnberg, roscoe, hietpas, lind, domingo_calap, carrasco, peris, sanjuan, maclean, rokyta )

# utility function for std error of mean
sem <- function(x) sqrt(var(x)/length(x))

# set up results table (one row of NAs that we'll discard later)
result <- matrix(,nrow=1,ncol=8)
colnames(result) <- c("Study", "MeanTi", "StdTi", "CountTi", "MeanTv", "StdTv", "CountTv", "U test" ) 

# compute results 
for (i in 1:length(all_dfes)) { 
	data <- all_dfes[[i]]
	titv <- as.factor(c(rep("ti",length(data$ti)),rep("tv",length(data$tv))))
	# "greater" indicates one-sided test with alt ::= x shifted to the right of y
	# Apparently x refers to ti because it is the factor level with the lower numeric encoding  
	out <- wilcox_test(c(data$ti,data$tv) ~ titv, alternative="greater")
	result <- rbind( result, c(data$name, mean(data$ti), sem(data$ti), length(data$ti), mean(data$tv), sem(data$tv), length(data$tv), pvalue(out) ) ) 
}

# first row in result is NA, so leave it out
result <- result[-1,]
write.table( result,"DFE_tests/DFE_U_tests.tab", sep="\t", row.names=FALSE )

# now let's do the other tests 
# define another results table 
meta_result <- matrix(,nrow=1,ncol=6)
colnames(meta_result) <- c("Group", "n", "Fisher_Xsq", "Fisher_P", "Stouffer_Z", "Stouffer_P" )

# define a test function 
#
test_meta <- function( group, result, meta_result ) {
	# set up 
	names <- lapply( group, function(x) { x$name })
	test_name <- do.call("paste",c(names,sep="_"))
	p_values <- unlist(lapply( names, function(x) { as.numeric(result[ result[,1]==x, 8 ]) }))
	sizes <- unlist(lapply( names, function(x) { sum(as.numeric(result[ result[,1]==x, c(4,7) ])) }) )

	# run tests 
	fisher <- Fisher.test(p_values)
	stouffer <- Stouffer.test(p_values, sizes)

	# output an updated meta table 
	rbind( meta_result, c(test_name, length(p_values), fisher["Xsq"], fisher["p.value"], stouffer["Z"], stouffer["p.value"] ) )
}

# run tests on groups (sorry didn't streamline this)
# studies that use competition and study the whole DFE
comp_subgroup <- list( lind, carrasco, domingo_calap, peris, sanjuan )
meta_result <- test_meta( comp_subgroup, result, meta_result )

# high-throughput studies that use a surrogate for fitness
ht_subgroup <- list( firnberg, roscoe, hietpas )
meta_result <- test_meta( ht_subgroup, result, meta_result )

# studies that use competition and are only looking at beneficial mutations
adaptive_subgroup <- list( maclean, rokyta )
meta_result <- test_meta( adaptive_subgroup, result, meta_result )

# all studies 
meta_result <- test_meta( all_dfes, result, meta_result )

write.table( meta_result[-1,],"DFE_tests/DFE_meta_tests.tab", sep="\t", row.names=FALSE )

# end 