# aa_wt_matrix_functions

source("aa_matrices_metadata.R")

# reading files - this is a bit convoluted but effective.  Note that there is never 
#     any attempt to guess the format of a file.   For a file to be parsed 
#     it needs a metadata block declared in aa_matrices_metadata.R
# 
# approach: 
#  first use readLines to navigate variable number of header lines
#  (warn=FALSE stifles warnings about not finding a return at the end of file) 
#  then use read.table for its handy "fill" feature (specify num cols via length 
#   of colnames); 
#  then convert to matrix because previous command actually produces an item=textrow list

# note:
# it should be possible to re-write this without readLines, using the "skip" and "nrows" arguments to read.table

read_aa_matrix <- function( dat ) { 
	start <- dat$file.header.lines + 1
	if ( dat$is.square | dat$diag.present ) { 
		n_rows <- 20
		colnames <- dat$aa.order
		rownames <- dat$aa.order
	}
	else { 
		n_rows <- 19
		colnames <- dat$aa.order[1:19]
		rownames <- dat$aa.order[2:20]
	}
	end <- start + n_rows - 1
	lines <- readLines( dat$file, warn=FALSE )[start:end]
	tmp <- read.table(textConnection(lines), fill=T,col.names=colnames, colClasses="numeric")
	tmp <- matrix( unlist(tmp), ncol=n_rows )
	rownames(tmp) <- rownames
	colnames(tmp) <- colnames
	return( tmp )
}

# implicitly, a square matrix is asymmetric.  a square matrix with symmetric values 
# is not an error condition, but just a square matrix with symmetric values.    

flattened_aa_matrix <- function( aa_matrix ) 
{ 
	# read matrix from file, flatten it
	aa_matrix$mat <- read_aa_matrix( aa_matrix )
	tmp <- as.data.frame(as.table(aa_matrix$mat))
	tmp <- tmp[!is.na(tmp[,3]),]
	# set the ordered 2-char code for this_AA --> that_AA
	if ( aa_matrix$is.square ) { 
		labels <- apply(tmp,1,function(x){paste(x[1],x[2],sep='')})
	}
	else { 
		labels <- apply(tmp,1,function(x){a<-sort(c(x[1],x[2]));paste(a[1],a[2],sep='')})
	}
	tmp <- as.matrix(tmp[,3])
	rownames(tmp) <- labels
	colnames(tmp) <- "Weight"
	return(tmp)
}

rescale_distance_to_similarity <- function ( weight ) 
{
	upper <- sum(range(weight))
	return( upper - weight ) 
}

compute_dim <- function( num_elements, is_square, diag_present ) { 
# assume a square matrix has diagonal values, even if null or empty
	if ( is_square ) { 
		return( sqrt( num_elements ) )
	}
# the dimension of a triangle matrix can be computed from this formula
	d <- ( sqrt(8 * elements + 1) - 1 )/2
	if ( diag_present ) { 
		return( d )
	}
	else { 
		return( d + 1 )
	}
}


# end of file
