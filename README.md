README for analysis of ti:tv bias
=================================

This file provides information on the set of files provided with Stoltzfus & Norris, 2014 (in prep).  

These files contain data and computer code in R, a data analysis environment that is widely used by statisticians (and bioinformaticians), and is available for all major platforms for free.  To find out more about R, go to R-project.org.  

The main point here is to carry out basic calculations to address various versions of a transition:transversion ratio for non-synonymous changes.  The scripts could be modified to address some other questions, but they are currently written to focus on the expected ratio of transitions to transversions among non-synonymous singlets (1-nt changes), using the canonical genetic code.   

CONTENTS
========

The contents should be obvious from the names of things.  At the top level are various R scripts that you might use.  See below. There is an "aa_matrices" subdirectory for input files for amino acid matrices, and a "codon_usage" subdirectory for codon usage tables.  

Data and code for carrying out tests on experimental studies of mutant effects are in "DFE_tests" (distribution of fitness effects).  

INSTALL 
=======

To use the files, just put them in the directory where you want to work, then launch R.  Typing "source(<file>)" will load the contents of the file into the R environment. 

There are a series of tests that were used during development of this package to be sure that files are ready correctly and calculations are performed correctly.  

TESTS 
=====

Executing the following 

> source("run_tests.R")

will run a test suite.  The most important bits for me to verify were some rather complicated counting operations on the genetic code.  These tests are in tests/1.R.  Read that file to get a sense of what is being verified.  

The other tests (2.R and 3.R) are not fully implemented.  They are intended to test data processing for amino acid matrices and codon usage tables.  I have done a number of spot tests to ensure that numbers processed out of foreign input files are as expected. 

CALCULATIING TABLES 1, 2 and 4 
==============================

The main results of the paper are analyses of studies of the distribution of fitness effects (DFE), or distribution of other mutant effects (i.e., surrogates of fitness) shown in tables 1, 2 and 4.  Executing this command 

> source( "run_DFE_tests.R" )

generates the output file "DFE_U_tests.tab". Note that "run_DFE_tests.R" is at the top level, but all the other files are in the "DFE_tests" subdirectory. The same script runs the meta-statistical tests (output in "DFE_meta_tests.tab") referenced in the text of Stoltzfus & Norris. 

The DFE_tests subdirectory also has a spreadsheet ("DFE_results_formatted_table.xslx") that will generate a formatted table. Just paste the results in "DFE_U_tests.tab" into the "raw" sheet and they will appear formatted in another sheet.  

Adding a study is simple (once you have processed the data from sources provided by the authors of the study, which can be difficult).   First, add a stanza to DFE_data.R that includes the data.  Then you need to edit "run_DFE_tests.R" in two places.  Adding your study to this line: 
 
> all_dfes <- list( firnberg, starita, roscoe, . . .  

will ensure that a Mann-Whitney U test is performed.  But if you want to group this study with other studies in a meta-analysis, you need to look toward the end of the "run_DFE_tests.R" file and also add the new study to one of the groups mentioned there. 

CALCULATIING TABLE 3 
====================

Table 3 shows the weighted ratios (weighted by codon use and amino acid exchangeabilities) of transitions to transversions calculated from a codon model.  Executing this command 

> source("run_titv_ratio_calculations.R")

will generate all the data for Table 3, which will appear in a text file named "weighted_titv_ratios.tab".  You can paste this into the provided spreadsheet ("titv_ratios_formatted_table.xlsx") and get a formatted table. 

To add your own cost matrix and have it used in these calculations, you first need to (1) put the file in the aa_matrices subdirectory (the other files will indicate which formats are allowable) and (2) add a stanza to the aa_matrices_metadata file, using the other examples to guide you.  That is, if your input file is formatted like the Grantham file, then copy the grantham stanza to use as a template.  Once you have these things done, you just need to specify in "run_titv_ratio_calculations.R" that you want to include your new matrix. This is the command to edit: 

> mats_to_test <- list( uniform_aa, exs, tangs_u, blosum, pam, lg )

To add a different codon usage is similar.  First, add the file in CUTG database format.  Next, put a stanza in the codon_usage_metadata file.  Then, at the end of the metadata file, edit this line: 

> all_cus <- list(uniform_cu, lambda, hiv1, e_coli, . . . 

so that it includes the new species.  The "run_titv_ratio_calculations.R" automatically iterates over every species included in all_cus. 

ADAPTING OR EXTENDING THIS CODE FOR OTHER CALCULATIONS
======================================================

This code is written in a somewhat general way, but it is really focused on one paper and has not been tested for anything else.  It would be great if you can find other ways to use the code, but please be cautious, and ask Arlin if things don't seem to be working.  

See notes for Table 1 above if you want to add an amino acid cost matrix or codon usage from a different species.  See the notes for Tables 2, 3 and 4 if you want to add another DFE study. 

