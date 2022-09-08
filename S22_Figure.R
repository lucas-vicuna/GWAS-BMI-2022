library(readr)
library(dplyr)
library(stringr)

setwd('~/proyectos/')

df <- read_tsv('bmi/Adiposity_Rebound/data/data_AR_trans.tsv')
# f1 <- 'bmi/Adiposity_Rebound/nuevos_analisis/qqplot_AgeAR.pdf'
# f2 <- 'bmi/Adiposity_Rebound/nuevos_analisis/qqplot_BMIAR.pdf'


# df1 <- df %>% filter(abs(zAge_AR)<3)
# 
# pdf(f1) 
# qqnorm(df$zAge_AR, main = '');qqline(df$zAge_AR, col = 2)
# dev.off()
# 
# 
# df2 <- df %>% filter(abs(zBMI_AR)<3)
# pdf(f2) 
# qqnorm(df$zBMI_AR, main = '');qqline(df$zBMI_AR, col = 2)
# dev.off()

df1 <- df %>% filter(abs(zAge_AR)<3)
df2 <- df %>% filter(abs(zBMI_AR)<3)
f1 <- 'bmi/Adiposity_Rebound/nuevos_analisis/qqplot_AgeAR_BMIAR.pdf'
pdf(file=f1, width=10, height=10)
par(mfrow=c(2,1), mar=c(3,3,3,3))

f1 <- 'bmi/Adiposity_Rebound/nuevos_analisis/qqplot_AgeAR_BMIAR2.pdf'
pdf(file=f1, width=5, height=10)
par(mfrow=c(2,1))
qqnorm(df$zAge_AR, main = '', xlab="Theoretical Quantile", ylab='Sample Quantile');qqline(df$zAge_AR, col = 2)
qqnorm(df$zBMI_AR, main = '', xlab="Theoretical Quantile", ylab='Sample Quantile');qqline(df$zBMI_AR, col = 2)
dev.off()
