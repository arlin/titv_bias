# run_tests.R

# following example in 
# http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/

# libraries
library('RUnit')  
library('Biostrings') 

# sources for this project
source('titv_functions.R')
# source('aa_wt_matrix_functions.R')
# source('codon_usage_tables.R')

test.suite <- defineTestSuite("TiTv project tests",
								dirs = file.path("tests"),
								testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)