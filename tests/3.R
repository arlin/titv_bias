# methods to test handling codon frequency tables

source("codon_usage_metadata.R")
source("codon_usage_functions.R")

# read in the file, compare the R-processed matrix with the raw file accessed via system (grep)
#   * note that the grep test based on a codon string is not necessarily safe 
#   * and the grep test assumes the GCC style of codon table output from CUTG database
#   * a possibly better test would be to read raw data via url from CUTG (see note in codon_usage_functions)

# default
default_species <- human
default_num <- 10

# we'll allow fractional slop in equality test, because input file has low-precision freqs
# (e.g., "12.17" for freq per 1000 = 0.01217) & we recompute high-precision freqs from counts
slop <- 0.05

test.codon.freq.read <- function( species=default_species, num=default_num ) 
{
	# step 1: read file using my R method 
	cu <- read_codon_usage( species )

	for ( iter in 1:num ) { 
		# step 2: pick a random value 
		codon <- sample(rownames(cu),1)

		# step 3: compare processed matrix with raw file using unix 
		this <- cu[rownames(cu)==codon]

		command <- paste("grep", codon, species$file, "| tr -s ' ' | cut -d ' ' -f4", sep=" ")
		that <- as.numeric(system(command, intern=TRUE))/1000

		message <- paste("Count for codon", codon, "from processed matrix matches grep-cut result" )
		out <- checkTrue(abs(this - that)/that < slop, message)
		print( paste( codon, this, "=", that, out, sep=" " ) )
	}
}



