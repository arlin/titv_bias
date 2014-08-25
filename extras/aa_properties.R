# aa_properties_data
# 
# we're just going to create a matrix of 0s and 1s to indicate conservation-- 
#   a replacement is charge-conserved if is_charged(AA1) = is_charged(AA2)
#   a replacement is hyrophobicity-conserved if is_hydrophobic(AA1) = is_hydrophobic(AA2)

charged <- strsplit( "RKDE", NULL )[[1]]
polar <- strsplit( "QNHSTYCMW", NULL)[[1]] 
hydrophobic <- strsplit( "AILFVPG", NULL)[[1]]
all_aas <- strsplit("ARNDCQEGHILKMFPSTWYV",NULL)[[1]]

get_aa_props_matrix <- function() { 
	pairs <- expand.grid(src=all_aas, dest=all_aas)
	labels <- apply(expand.grid(src=all_aas, dest=all_aas),1, function(x) { paste(x[1],x[2], sep="") })
	
	charge_cons <- apply(pairs, 1, function(x) { if( x[1] %in% charged == x[2] %in% charged ) {return(1)} else {return(0)} })

	hydrop_cons <- apply(pairs, 1, function(x) { if( x[1] %in% hydrophobic == x[2] %in% hydrophobic ) {return(1)} else {return(0)} })

	aa_props <- matrix( c(charge_cons,hydrop_cons), ncol=2, byrow=FALSE)
	rownames(aa_props) <- labels
	colnames(aa_props) <- c("charge_cons", "hydrop_cons") 
	return( aa_props ) 
}

# end of file 


aa_repl_wts <- flattened_aa_matrix( aa_matrix )
aa_props <- get_aa_props_matrix() 

# merge this with a matrix by row 
dat <- merge(aa_props, aa_repl_wts, by.x="row.names",by.y="row.names")

# run the test 
out <- wilcox_test(dat[,"Weight"] ~ as.factor(dat[,"charge_cons"]))
