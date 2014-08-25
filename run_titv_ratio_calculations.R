# compute titv ratios weighted by codon use & aa exchangeability 

# to compute the results in Table 1, 
#
# > source("run_titv_ratio_calculations.R")
#
# NOTE that the point of this is to see if a weighting increases the ti/tv ratio. 
# To ensure that the ti/tv ratio will always increase when ti's are favored, 
# distance matrices (e.g., Grantham) are converted to similarities, and log-odds 
# match scores (e.g., PAM, BLOSUM) are exponentiated.  

source("aa_wt_matrix_functions.R")
source("titv_ratio_functions.R")
source("codon_usage_functions.R")

mats_to_test <- list( uniform_aa, exs, tangs_u, blosum, pam, lg )

# set up our results tables
result <- matrix(,nrow=1,ncol=3)
colnames(result) <- c("AA_matrix", "CodonUsage", "TiTvRatio") 
result2 <- matrix(,nrow=1,ncol=5)
colnames(result2) <- c("AA_matrix", "Uniform CU", "CU mean", "CU low", "CU high") 

# this is slow step so we'll do it once then re-use
paths <- compute_codon_paths()

for (i in 1:length(mats_to_test)) { 
	out <- iterate_matrix_with_cus( mats_to_test[[i]], all_cus, paths )
	result <- rbind(result, out)
	uniform <- out[out[,2]=="Uniform",3]
	mean <- mean(as.numeric(out[out[,2]!="Uniform",3]))
	range <- range(as.numeric(out[out[,2]!="Uniform",3]))
	result2 <- rbind(result2, c(mats_to_test[[i]]$name, uniform, mean, range) )
}

write.table( result[-1,],"weighted_titv_ratios_raw.tab", sep="\t", row.names=FALSE )
write.table( result2[-1,],"weighted_titv_ratios.tab", sep="\t", row.names=FALSE )

# end 