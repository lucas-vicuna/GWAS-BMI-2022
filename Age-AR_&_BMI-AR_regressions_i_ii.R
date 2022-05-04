# This script was used to perform the GWAS on Age-AR and BMI-AR, including 2 combinations of covariates: 1. With local ancestry, without global ancestry and without the first 10 genetic PCs; and 2. Without local ancestry, with global ancestry and without the first 10 genetic PCs. 

library(data.table)
library(dplyr)
library(broom)
library(readr)

args <- commandArgs(TRUE)
chr <- args[1]

link <- '~/proyectos/'
setwd(link)

file <- paste0('path/data_AR_trans.tsv')
file2 <- paste0('path/chr',chr,'_SNP_selected.tsv')

f1 <- 'path/educacion_madre2.csv'

file3 <- paste0('/path/rfmix15.chr',chr,'.stats.phen.txt')
file4 <- paste0('path/AgeAR_i_chr',chr,'.txt')
file5 <- paste0('path/BMIAR_i_chr',chr,'.txt')
file6 <- paste0('path/AgeAR_ii_chr',chr,'.txt')
file7 <- paste0('path/BMIAR_ii_chr',chr,'.txt')
file8 <- paste0('path/chr',chr,'_AR_failed.txt')

df <- fread(file, sep='\t', header=T)
df <- df %>% select(CODE, Sex, Age_AR, BMI_AR, zAge_AR, zBMI_AR) %>% mutate(Sex = factor(Sex))

df1 <- read_csv(f1)
df <- df %>% inner_join(df1) %>% mutate(edu=factor(edu))

df2 <- fread(file2,header = T,sep='\t')
k <- df2 %>% nrow


df3 <- fread(file3,
             sep='\t', 
             select =c(1,6,df2$POSITION, df2$POSITION+1),
             header=T)
colnames(df3)[1:2] <- c('CODE','ga')

df4 <- df %>% inner_join(df3)

for (i in 1:k){
  
  temp <- df4[,.SD,.SDcols=c(1:8,i+8,i+8+k)]
  names(temp)[9:10] <- c('geno_snp', 'locanc_snp')
  
  tryCatch({
  
     #i With local ancestry, without global ancestry, without first 10 PCs, with maternal education level
    
    m_Age_i <- lm(zAge_AR ~ poly(zBMI_AR,3):Sex + Sex + edu + geno_snp:Sex + locanc_snp, data = temp %>% filter(abs(zAge_AR)<3))
    m_BMI_i <- lm(zBMI_AR ~ poly(zAge_AR,3):Sex + Sex + edu + geno_snp:Sex + locanc_snp, data = temp %>% filter(abs(zBMI_AR)<3))
    
     #ii Without local ancestry, with global ancestry, without first 10 PCs, with maternal education level
 
    m_Age_ii <- lm(zAge_AR ~ poly(zBMI_AR,3):Sex + Sex + edu + geno_snp:Sex + ga, data = temp %>% filter(abs(zAge_AR)<3))
    m_BMI_ii <- lm(zBMI_AR ~ poly(zAge_AR,3):Sex + Sex + edu + geno_snp:Sex + ga, data = temp %>% filter(abs(zBMI_AR)<3))
    
    m1 <- bind_cols(SNP = df2$SNP[i], m_Age_i %>% tidy)
    m2 <- bind_cols(SNP = df2$SNP[i], m_BMI_i %>% tidy)
    m3 <- bind_cols(SNP = df2$SNP[i], m_Age_ii %>% tidy)
    m4 <- bind_cols(SNP = df2$SNP[i], m_BMI_ii %>% tidy)

    write.table(m1,file=file4,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')
    write.table(m2,file=file5,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')
    write.table(m3,file=file6,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')
    write.table(m4,file=file7,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')

    },
    error=function(e){ write.table(df2$POSITION[i],file=file12,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
  
}
