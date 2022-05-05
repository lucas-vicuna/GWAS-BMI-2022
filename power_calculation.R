
# This script was used to estimate the statistical power of the cross-sectional GWAS on BMI to detect significant associations. 
# We followed the indications provided here: https://genome.sph.umich.edu/wiki/Power_Calculations:_Quantitative_Traits

# In the following example, we estimated the statistical power of associations between rs269511 and BMI across age stratum 15.5-16.5 years 
# (see S5 Table of the manuscript). 

# N: number of individuals
# M: number of independent genetic variants examined
# alpha: significance level
# H2: effect size of the genotype

N = 587
M = 656059
alpha = 0.05/M
H2 = 0.0070804356

threshold = qchisq(alpha, df = 1, lower.tail = FALSE)
power = pchisq(threshold, df = 1, lower.tail = FALSE, ncp = N * H2)
power
