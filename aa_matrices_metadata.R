# aa_matrices_metadata - metadata file for aa matrices 
# 
# this is the only file you need to change if you want to include a new matrix. 
# 
# To add a new matrix, 
#   1. store the matrix file locally 
#   2. fill out a metadata list (include path to matrix file)
# 
# Notes: 
# 
# name = the display name to be used to denote this matrix in any output 
#
# file.header.lines = number of lines to be skipped from start of file 
# 
# is.square = the data in $file are given in a square (not triangle) matrix 
# 
# diag.present = the triangle matrix in $file has a diagonal (not relevant for 
#   square matrices, which always have a diagonal, though it may contain blank 
#   values or NAs)
# 
# is.log.odds = values are log of odds-ratio for true vs. false match.  We need 
#    to exponentiate these values in order to compute a comparable ti/tv ratio
#
# is.distance = values are distances, negatively correlated with similarity

aa_order_3lc <- strsplit("ARNDCQEGHILKMFPSTWYV",NULL)[[1]]
aa_order_chem <- strsplit("CSTPAGNDEQHRKMILVFYW",NULL)[[1]]

# uniform (all 1s for testing purposes) 
uniform_aa <- list( 
	name="Uniform",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="a matrix all of 1s, for testing and no-effect controls", 
	file="aa_matrices/uniform.dat",
	file.header.lines=0 )
# end of uniform_aa 

# test (various values for testing purposes) 
test <- list( 
	name="Test",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="a matrix with assigned values for testing purposes", 
	file="aa_matrices/test.dat",
	file.header.lines=0 )
# end of test 

test_sq <- list( 
	name="Test (square)",
	is.square=TRUE,
	aa.order=aa_order_chem,
	diag.present=TRUE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="a matrix with assigned values for testing purposes", 
	file="aa_matrices/test_sq.dat",
	file.header.lines=0 )
# end of test 

# EX
ex <- list( 
	name="EX",
	is.square=TRUE,
	aa.order=aa_order_chem,
	diag.present=TRUE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="Asymmetric EX matrix of Yampolsky & Stoltzfus, 2005", 
	file="aa_matrices/ex.dat",
	file.header.lines=0 )
# end of ex 

# EX uncertainty
ex_uncertainty <- list( 
	name="EX uncertainty",
	is.square=TRUE,
	aa.order=aa_order_chem,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="Bootstrapped uncertainty values for EX (Yampolsky & Stoltzfus, 2005)", 
	file="aa_matrices/ex_uncertainty.dat",
	file.header.lines=0 )
# end of ex_uncertainty

# EXS
exs <- list( 
	name="EXS",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="Symmetric version of EX matrix of Yampolsky & Stoltzfus, 2005", 
	file="aa_matrices/excs.dat",
	file.header.lines=0 )
# end of exs 

# Tang's U matrix  
tangs_u <- list( 
	name="U",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="U matrix of Tang, et al., 2004.  Data from Table 2.  Missing values for doublets & triplets are given as NA",
	file="aa_matrices/u.dat",
	file.header.lines=0 )
# end of tangs_u 

# LG matrix  
lg <- list( 
	name="LG",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=FALSE,
	meta="source: http://www.atgc-montpellier.fr/download/datasets/models/lg_LG.PAML.txt
		downloaded: April 29, 2014
		citation: Le SQ, Gascuel O. 2008. Mol Biol Evol. 2008 Jul;25(7):1307-20. 
		doi: 10.1093/molbev/msn067
		PubMed link: http://www.ncbi.nlm.nih.gov/pubmed/18367465?dopt=Abstract", 
	file="aa_matrices/lg.dat",
	file.header.lines=6 )
# end of lg 

# Grantham matrix 
grantham <- list( 
	name="Grantham",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=TRUE,
	meta="from the original paper by Grantham, entered by Arlin", 
	file="aa_matrices/gra-orig.dat",
	file.header.lines=0 )
# end of grantham 

# Miyata matrix 
miyata <- list( 
	name="Miyata",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=FALSE,
	is.log.odds=FALSE,
	is.distance=TRUE,
	meta="", 
	file="aa_matrices/miyata.dat",
	file.header.lines=0 )
# end of miyata

# BLOSUM62 matrix 
blosum <- list( 
	name="BLOSUM62",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=TRUE,
	is.log.odds=TRUE,
	is.distance=FALSE,
	meta="This from AAIndex.  Metadata are included in the source file.", 
	file="aa_matrices/blo62.dat",
	file.header.lines=8 )
# end of blosum

# BLOSUM62 matrix  
blosum2 <- list( 
	name="BLOSUM62 (Biostrings)",
	is.square=TRUE,
	aa.order=aa_order_3lc,
	diag.present=TRUE,
	is.log.odds=TRUE,
	is.distance=FALSE,
	meta="this is the BLOSUM62 matrix included in the biostrings package", 
	file="aa_matrices/other_blo62.dat",
	file.header.lines=0 )
# end of blosum2


# PAM250 from Benner, et al 
pam <- list( 
	name="PAM250",
	is.square=FALSE,
	aa.order=aa_order_3lc,
	diag.present=TRUE,
	is.log.odds=TRUE,
	is.distance=FALSE,
	meta="This is from AAIndex.   Metadata are included in the source file.", 
	file="aa_matrices/pam250.dat",
	file.header.lines=9 )
# end of pam
# end of file 
