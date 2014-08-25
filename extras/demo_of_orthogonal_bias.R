# a question that came up during this study was whether the DFE for a sample of 
# transitions and transversions will be biased if the sampling had a mutational
# bias.  For instance, MacLean, et al. are studying adaptive mutants that 
# arise in the lab via a process that is biased by mutation.  Transitions occur 
# at higher rates so we see more of them in such studies.  Will this bias the 
# DFE?  Arlin initially thought that, because transitions occur at higher rates, 
# they have an unfair advantage of being recovered, and that beneficial 
# transitions recovered from mutation-driven evolution would be *less* 
# beneficial (on average) than transversions.  David McCandlish was skeptical 
# that this would turn out to be correct.  
#
# The code below generates a landscape of mutational possibilities including 
# transitions and transversions, each with a fitness drawn from a similar 
# distribution, then samples from it using origin-fixation dynamics.  This 
# is the same as the "mutational landscape" concept of Gillespie and Orr, 
# though this particular version allows beneficial, neutral and deleterious 
# fixations according to a standard formula for the probability of fixation.  
#
# The results indicate that David was right to be skeptical.  There is no 
# effect of biased sampling.  We can explain it like this.  Under the null 
# hypothesis, transitions and tranversions have the same fitness distribution. 
# Biased sampling via mutation results in more transitions, but the mutation 
# bias is an orthogonal factor.  We get more transitions, but not a shift in 
# the distribution of fitness effects.    

# need this library
library(coin)

# defaults (size = number of mutational possibilities)
def_Ne <- 100000
def_ti_bias <- 3
def_size <- 10000
def_sample_size <- 1000

# probability of fixation 
p_fix <- function( s, N ) { 2*s / (1 - exp(-4*N*s))}

# run test, e.g., run_test(), run_test(ti_bias=5), run_test(ti_bias=5, size=1000)
#
run_test <- function( Ne=def_Ne, ti_bias=def_ti_bias, size=def_size, sample_size=def_sample_size ) {

	# selection coefficients for tis and tvs 
	ti_s <- rnorm(size,0,sd=0.01)
	tv_s <- rnorm(size,0,sd=0.01)

	# fixation probabilities
	ti_p <- p_fix( ti_s, Ne)
	tv_p <- p_fix( tv_s, Ne)

	# the mut landscape probabilities for all ti & tv mutations, under origin-fixation, 
	# given ti_bias and the above probabilities of fixation, normalized to sum to 1
	mut_landscape <- c( ti_p*ti_bias, tv_p)/sum(ti_p*ti_bias, tv_p)
	
	# think of these as the other columns in the mut landscape table, giving the index, 
	# the mutation type, and the selection coefficient
	indices <- 1:(2*size)
	sel_coeffs <- c( ti_s, tv_s )
	mut_types <- as.factor( c(rep("ti",size), rep("tv",size)))

	# now we'll take a sample and analyze it 
	out <- sample(indices,sample_size,replace=TRUE,prob=mut_landscape)

	# see how many tis and tvs we got
	tab <- table(mut_types[out])

	# test for equality 
	test <- wilcox_test( sel_coeffs[out] ~ mut_types[out])

	# output
	c(tab, test)
}
