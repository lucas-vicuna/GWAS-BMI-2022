# This script was used to perform the GWAS on Age-AR and BMI-AR, including 2 combinations of covariates: 3. With local ancestry, without global ancestry and with the first 10 genetic PCs; and 4. Without local ancestry, without global ancestry and with the first 10 genetic PCs. 

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

f12 <- 'path/educacion_madre2.csv'
f13 <- 'path/GOCS_904PCA.txt'

file3 <- paste0('path/rfmix15.chr',chr,'.stats.phen.txt')
file4 <- paste0('path/AgeAR_iii_chr',chr,'.txt')
file5 <- paste0('path/BMIAR_iii_chr',chr,'.txt')
file6 <- paste0('path/AgeAR_iv_chr',chr,'.txt')
file7 <- paste0('path/BMIAR_iv_chr',chr,'.txt')
file8 <- paste0('path/chr',chr,'_AR_failed.txt')

df <- fread(file, sep='\t', header=T)
df <- df %>% select(CODE, Sex, Age_AR, BMI_AR, zAge_AR, zBMI_AR) %>% mutate(Sex = factor(Sex))

df12 <- read_csv(f12)
df13 <- read_delim(f13, delim='\t',col_names = F)
names(df13) <- c('CODE',paste0('X',1:10))

df <- df %>% inner_join(df12) %>% inner_join(df13) %>% mutate(edu=factor(edu))

df2 <- fread(file2,header = T,sep='\t')
k <- df2 %>% nrow

df3 <- fread(file3,
             sep='\t', 
             select =c(1,6,df2$POSITION, df2$POSITION+1),
             header=T)
colnames(df3)[1:2] <- c('CODE','ga')

df4 <- df %>% inner_join(df3)

for (i in 1:k){
  
  temp <- df4[,.SD,.SDcols=c(1:18,i+18,i+18+k)]
  names(temp)[19:20] <- c('geno_snp', 'locanc_snp')
  
  tryCatch({
    
    #iii With local ancestry, without global ancestry, with first 10 PCs, with maternal education level
    
    m_Age_iii <- lm(zAge_AR ~ poly(zBMI_AR,3):Sex + Sex + edu + geno_snp:Sex + locanc_snp+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = temp %>% filter(abs(zAge_AR)<3))
    m_BMI_iii <- lm(zBMI_AR ~ poly(zAge_AR,3):Sex + Sex + edu + geno_snp:Sex + locanc_snp+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = temp %>% filter(abs(zBMI_AR)<3))
    
    #iv Without local ancestry, without global ancestry, with first 10 PCs, with maternal education level
    
    m_Age_iv <- lm(zAge_AR ~ poly(zBMI_AR,3):Sex + Sex + edu + geno_snp:Sex + X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = temp %>% filter(abs(zAge_AR)<3))
    m_BMI_iv <- lm(zBMI_AR ~ poly(zAge_AR,3):Sex + Sex + edu + geno_snp:Sex + X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = temp %>% filter(abs(zBMI_AR)<3))
    
    m1 <- bind_cols(SNP = df2$SNP[i], m_Age_iii %>% tidy)
    m2 <- bind_cols(SNP = df2$SNP[i], m_BMI_iii %>% tidy)
    m3 <- bind_cols(SNP = df2$SNP[i], m_Age_iv %>% tidy)
    m4 <- bind_cols(SNP = df2$SNP[i], m_BMI_iv %>% tidy)

    write.table(m1,file=file4,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')
    write.table(m2,file=file5,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')
    write.table(m3,file=file6,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')
    write.table(m4,file=file7,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE, sep='\t')

    },
    error=function(e){ write.table(df2$POSITION[i],file=file12,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
  
}
