# tests - amino acid weight matrices

# intro to RUnit test module: 
#  http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/

# define some test values to be sure that we have read matrices correctly 

# get singlet codon paths, reduce to unique ordered aa paths
source('titv_functions.R') 
paths <- compute_codon_paths()
paths <- paths[paths[,'from_Aa']!=paths[,'to_Aa'] & paths[,'to_Aa'] != '*' & paths[,'from_Aa'] != '*',]
singlet_all <- paths[!duplicated(paths[,c('from_Aa','to_Aa')]),c('from_Aa','to_Aa')]

# get a set to represent the 75 unordered pairs (reverse = forward[,2:1])
two.letter <- apply(singlet_all,1,function(x) { a <- sort(x); paste(a[1],a[2],sep="")})
singlet_forward <- singlet_all[!duplicated(two.letter),] 

# randomly chosen pairs
test_pairs <- matrix(strsplit("SADGLRCLYDIM", NULL)[[1]], ncol=2,nrow=6)
blo62_test_values <- c(-1, -1, -3, -1,  2, -1)

# singlet occupancy = singlet pairs with defined values (max=150)
# non-null values for every pair in singlet_all

# symmetry = (weight[a,b]==weight[b,a])
# compare values for singlet_for with those for singlet_rev

# sum = sum of all values in matrix
# check the sum of the matrix against a known value of sum

# test values = values for 6 pairs
# check known value for 6 pairs chosen at random from 150 ordered singlet pairs



######
# matrix of AA distances are read correctly 
######
test.aa.matrices <- function () 
{
# EX matrix is read correctly 

# U matrix is read correctly 

# codon usage is read correctly 

# PAM matrix is read correctly 

# BLO matrix is read correctly 

# GRA matrix is read correctly 

}
