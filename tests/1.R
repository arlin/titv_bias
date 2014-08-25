# tests - codon paths 

# intro to RUnit test module: 
#  http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/

######
# counts of unique changes by codon-to-codon paths
######
test.codon.path.counts <- function() 
{
# get paths 
	paths <- compute_codon_paths()
	
# SYNONYMOUS TRANSITIONS: 
# The object is to count all codon changes that do not change the amino acid. There are
# 4*4*2 [TCAG][TCAG][RY] = 32 half-blocks in the code.  The triplets in each half-block 
# are interconvertible via a transition at pos=3, whereas there is no transversion within 
# a half-block.  Of the 32 possible half-blocks, 29 of them-- all but TGR, ATR, TAR-- 
# encode the same AA, which means there are 29*2 = 58 synonymous ti at pos=3.  To this,
# add 4 synonymous TTR<-->CTR Leucine mutations at pos=1. There are no synonymous ti other 
# than these 62.  

# there are 62 synonymous ti
	t <- paths[paths[,'from_Aa']==paths[,'to_Aa'] & paths[,'mut_Type']=='ti' & paths[,'from_Aa'] != '*',]
	checkTrue(dim(t)[1] == 62, 'Expect 62 synonymous ti paths.')

# synonymous transitions are 4, 0, and 58 by position 1, 2, and 3, resp
# (clumsy but effective) 
	a <- dim(t[t[,'codon_Pos']==1,])[1]
	b <- dim(t[t[,'codon_Pos']==2,])[1]
	c <- dim(t[t[,'codon_Pos']==3,])[1]
	checkTrue(a==4 & b==0 & c==58, 'Expect (4, 0, 58) synonymous ti paths by codon position.')

# SYNONYMOUS TRANSVERSIONS: 
# each of 8 full blocks (L,V,S,P,T,A,R,G) contributes 8 syn tv at pos=3. 
# Add 4 Isoleucine ATA<-->ATY at pos=3 and 4 Arginine CGR<-->AGR at pos=1

# there are 72 synonymous tv 
	t <- paths[paths[,'from_Aa']==paths[,'to_Aa'] & paths[,'mut_Type']!='ti' & paths[,'from_Aa'] != '*',]
	checkTrue(dim(t)[1] == 72, 'Expect 72 synonymous tv paths.')

# synonymous tvs are 4, 0, 68 by position 1, 2 and 3, resp.
# (clumsy but effective) 
	a <- dim(t[t[,'codon_Pos']==1,])[1]
	b <- dim(t[t[,'codon_Pos']==2,])[1]
	c <- dim(t[t[,'codon_Pos']==3,])[1]
	checkTrue(a==4 & b==0 & c==68, 'Expect (4, 0, 58) synonymous ti paths by codon position.')

# NONSENSE CHANGES: 
# It's easier to count the reverse class of nonsense-to-sense changes, and then invoke
# forward-reverse symmetry to impute the sense-to-nonsense class. Of the 3*9=27 total 
# mutations to nonsense codons, 4 change one nonsense to another, and all 4 are 
# transitions: TGA<-->TAA at pos=2, TAA<-->TAG at pos=3.  The remaining 27 - 4 = 23 
# are all nonsense-to-sense.  In terms of codon position, the { 9, 9, 9 } distribution 
# shifts to { 9, 7, 7 }.  In terms of ti:tv, the 9:18 distribution shifts to 5:18.  
# The combined distribution is { 3:6, 1:6, 1:6 }.  The sense-to-nonsense class is 
# simply the class of reversals of these 23 paths, therefore it has the same distribution. 

# there are 23 changes from sense to nonsense 
	t <- paths[paths[,'to_Aa'] == '*' & paths[,'from_Aa'] != '*',]
	checkTrue(dim(t)[1] == 23, 'Expect 23 sense-to-nonsense paths.')

# and these should be split 9, 7, 7 among 1st, 2nd, and 3rd positions
	a <- dim(t[t[,'codon_Pos']==1,])[1]
	b <- dim(t[t[,'codon_Pos']==2,])[1]
	c <- dim(t[t[,'codon_Pos']==3,])[1]
	checkTrue(a==9 & b==7 & c==7, 'Expect (9, 7, 7) nonsense paths by codon position.')

# REPLACEMENTS: 
# By subtraction, we know how many replacements to expect.  There are 61*9=549 
# total changes to sense codons, and we know (above) that 62+72=134 are syn 
# and 23 are nonsense, for a total of 157. This leaves 492 replacement changes.  
# Enumerating these is tedious and error-prone.  I have done it manually and 
# with a Perl script, with results converging on 116 ti (54, 60, 2 by phase) 
# and 276 tv (112, 116, 48 by phase).   

# there are 116 ti replacements 
	t <- paths[paths[,'from_Aa']!=paths[,'to_Aa'] & paths[,'to_Aa'] != '*' & paths[,'from_Aa'] != '*',]
	checkTrue(dim(t[t[,'mut_Type']=='ti',])[1] == 116, 'Expect 116 ti replacements.')
	
# ti replacements are 54, 60, 2 by positions 1, 2 and 3, resp.  
	a <- dim(t[t[,'codon_Pos']==1 & t[,'mut_Type']=='ti',])[1]
	b <- dim(t[t[,'codon_Pos']==2 & t[,'mut_Type']=='ti',])[1]
	c <- dim(t[t[,'codon_Pos']==3 & t[,'mut_Type']=='ti',])[1]
	checkTrue(a==54 & b==60 & c==2, 'Expect (54, 60, 2) ti replacements by codon position.')

# there are 276 tv replacements 
	checkTrue(dim(t[t[,'mut_Type']!='ti',])[1] == 276, 'Expect 276 tv replacements.')

# tv replacements are 112, 116, 48 by positions 1, 2 and 3, resp.  
	a <- dim(t[t[,'codon_Pos']==1 & t[,'mut_Type']!='ti',])[1]
	b <- dim(t[t[,'codon_Pos']==2 & t[,'mut_Type']!='ti',])[1]
	c <- dim(t[t[,'codon_Pos']==3 & t[,'mut_Type']!='ti',])[1]
	checkTrue(a==112 & b==116 & c==48, 'Expect (4, 0, 58) tv replacements by codon position.')

# FINAL CHECK: 
# Note that these numbers all add up so that we have the expected 61*3=183 
# changes at each position:  
#
# pos=1: 8 (4 ti + 4 tv) syn, 166 (54 ti + 112 tv) rep, 9 nonsense = 183
# pos=2: 0 syn, 176 (60 ti + 116 tv) rep, 7 nonsense = 183
# pos=3: 126 (58 ti + 68 tv) syn, 50 (2 ti + 45 tv) rep, 7 nonsense = 183 
# 
# Among these, there are 183 transitions and 2*183 transversion Above, we 
# didn't split the nonsense change by ti or tv, but the syn + rep 
# counts add up to 178 ti and 345 tv, thus the remaining 5 ti and 18 tv must 
# be in the sense-to-nonsense class.  This is consistent with above, when we 
# counted the reverse class of nonsense-to-sense changes: of 3*9=27 mutations, 
# we excluded 4 transitions (TGA<-->TAA, TAA<-->TAG) as occurring between 2 
# stop codons.  Thus, the nonsense-to-sense class has 5 ti and 18 tv, so also 
# the class of reversals (sense-to-nonsense) must have 5 ti and 18 tv.  
# 
# As one final test, we can look explicitly for this distribution

	t <- paths[paths[,'to_Aa'] == '*' & paths[,'from_Aa'] != '*',]
	a <- dim(t[t[,'mut_Type'] == 'ti',])[1]
	b <- dim(t[t[,'mut_Type'] != 'ti',])[1]
	checkTrue(a == 5 & b == 18, 'Expect 5 ti and 18 tv nonsense paths.')
}

######
# counts of unique changes by Aa-to-Aa paths
######
test.aa.path.counts <- function () 
{
# get paths 
	paths <- compute_codon_paths()
	pruned <- paths[paths[,'from_Aa']!=paths[,'to_Aa'] & paths[,'to_Aa'] != '*' & paths[,'from_Aa'] != '*',]

# there are 150 ordered Aa-to-Aa paths 
	pruned_more <- pruned[!duplicated(pruned[,c('from_Aa','to_Aa')]),]
	checkTrue(dim(pruned_more)[1] == 150, 'Expect 150 ordered Aa-to-Aa paths.')

# there are 158 ordered Aa-to-Aa paths when marked by ti or tv
	# first, simplify mut_Type to just 2 levels, ti and tv
	pruned[pruned[,'mut_Type']!='ti','mut_Type'] <- 'tv'
	pruned_more <- pruned[!duplicated(pruned[,c('from_Aa','to_Aa','mut_Type')]),]
	checkTrue(dim(pruned_more)[1] == 158, 'Expect 158 ordered Aa-to-Aa paths marked by ti or tv.')

# the exceptions (with both ti and tv paths) are F<-->L, G<-->R, I<-->M and R<-->W  
	exp <- c("FL", "GR", "IM", "RW" )
	# get the paths, turn them into sorted 2-letter strings, then sort-uniq those
	pruned_more <- pruned_more[duplicated(pruned_more[,c('from_Aa','to_Aa')]),c('from_Aa','to_Aa')]
	pruned_more <- sort(apply(pruned_more,1,function(x,y){a<-sort(c(x[1],x[2]));paste(a[1],a[2],sep='')}))
	pruned_more <- pruned_more[!duplicated(pruned_more)]
	checkTrue(!(FALSE %in% (pruned_more==exp)), paste('Expect 4 ambiguous Aa paths', exp[1],exp[2],exp[3],exp[4],sep=" "))
} 


######
# matrix of transitions and transversions is correct 
#####
test.titv.paths <- function() 
{
# (T<-->C, G<-->A) is a ti, all others are tv

}
# end 