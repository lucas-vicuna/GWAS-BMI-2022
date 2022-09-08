
# This script was used to estimate the statistical power of the cross-sectional GWAS on BMI to detect significant associations. 
# We followed the indications provided here: https://genome.sph.umich.edu/wiki/Power_Calculations:_Quantitative_Traits

# In the following example, we estimated the statistical power SNPs from Table 1 of the manuscript

# N: number of individuals
# M: number of independent genetic variants examined
# alpha: significance level
# H2: effect size of the genotype

library(purrr)
library(broom)
library(readr)
library(dplyr)
library(data.table)

mifuncion <- function(chr,snp,intervalo){
  
snp2=paste0('GT_',snp)

df <- read_delim('~/proyectos/bmi/xAge/TramosxEdad.csv') %>% 
  filter(INT==intervalo)

f2 <- paste0('/data/genetica/tomas_localancestry/gocs_pel_ceu_pel95_phen904/rfmix15/rfmix15.chr',chr,'.stats.phen.pel_ceu.txt')
df2 <- fread(f2, sep='\t', header=T,nrows = 1)
num_snp <- (1:ncol(df2))[names(df2)==snp2]              

df22 <- fread(f2, sep='\t', header=T,select=c(1,6,num_snp,num_snp+1))
colnames(df22) <- c('CODE','ga','geno_snp','locanc_snp')

df3 <- df %>% left_join(df22)

modelo <- lm(zBMI ~ Sex + ga + geno_snp + locanc_snp, data = df3)

beta_snp=modelo$coefficients[4]
M = 609122
alpha = 0.05/M
H2=summary(modelo)$r.squared
N=nrow(df3)
threshold = qchisq(alpha, df = 1, lower.tail = FALSE)
power = pchisq(threshold, df = 1, lower.tail = FALSE, ncp = N * H2)

A=tibble(Snp=snp,Chr=chr,Age_stratum=intervalo,Sample_size=N,Effect_size=beta_snp, Statistical_power=power)

return(A)
}

A2=bind_rows(
  mifuncion(chr=5,snp='rs269511',intervalo='(1.5,2.5]'),
  mifuncion(chr=6,snp='rs9275582',intervalo='(1.5,2.5]'),
  mifuncion(chr=6,snp='rs9275593',intervalo='(1.5,2.5]'),
  mifuncion(chr=6,snp='rs9275595',intervalo='(1.5,2.5]'),
  mifuncion(chr=12,snp='rs7134291',intervalo='(8.5,9.5]'),
  mifuncion(chr=10,snp='rs7896870',intervalo='(0.5,1.5]'),
  mifuncion(chr=20,snp='rs474169',intervalo='(10.5,11.5]'),
  mifuncion(chr=15,snp='rs1495271',intervalo='(4.5,5.5]'),
  mifuncion(chr=10,snp='rs11244839',intervalo='(0.5,1.5]'),
  mifuncion(chr=8,snp='rs13257360',intervalo='(14.5,15.5]')
)
A2 %>% View
A2 %>% 
write_csv('~/proyectos/bmi/xAge/power_snps/power_snps.csv')

