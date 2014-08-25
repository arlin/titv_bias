# codon usage functions

source("codon_usage_metadata.R") 
colnames <- c("AmAcid","Codon","Number","/1000","Fraction")

# note that we can read from URL as well 
# read.table(file=human$url,blank.lines.skip = TRUE,skip=15,col.names=colnames,nrows=64)

read_codon_usage <- function( species ) 
{
	cu <- read.table(file=species$file,blank.lines.skip = TRUE,skip=1,col.names=colnames)
	# cu[,5] <- cu[,3]/sum(cu[,3])
	freq <- as.matrix( cu[,3]/sum(cu[,3]) )
	colnames( freq ) <- "Frequency"
	rownames( freq ) <- cu[,2]
	return( freq )
}

# end 
