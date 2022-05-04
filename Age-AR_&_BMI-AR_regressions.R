# This script was used to perform the main GWAS on Age-AR and BMI-AR, including local ancestry and global ancestry, but excluding the first 10 PCs. 

library(data.table)
library(dplyr)
library(broom)

args <- commandArgs(TRUE)
chr <- args[1]

link <- '~/proyectos/bmi/'
setwd(link)

file <- paste0('path/data_AR_trans.tsv')
file2 <- paste0('path/chr',chr,'_SNP_selected.tsv')
file3 <- paste0('/path/rfmix15.chr',chr,'.stats.phen.txt')
file6 <- paste0('path/output Age_AR2/chr',chr,'_Age_AR2.txt')
file7 <- paste0('path/output BMI_AR2/chr',chr,'_BMI_AR2.txt')
file8 <- paste0('path/output Age_AR/chr',chr,'_AR_II.txt')

df <- fread(file, sep='\t', header=T)
df <- df %>% select(CODE, Sex, Age_AR, BMI_AR, zAge_AR, zBMI_AR) %>% mutate(Sex = factor(Sex))

df2 <- fread(file2,header = T,sep='\t')
k <- df2 %>% nrow

df3 <- fread(file3,
             sep='\t', 
             select =c(1,6,df2$POSITION, df2$POSITION+1),
             header=T)
colnames(df3)[1:2] <- c('CODE','ga')

df4 <- df %>% inner_join(df3)

for (i in 1:k){
  
  temp <- df4[,.SD,.SDcols=c(1:7,i+7,i+7+k)]
  names(temp)[8:9] <- c('geno_snp', 'locanc_snp')
  
  tryCatch({
    
    m_Age2 <- lm(zAge_AR ~ BMI_AR + Sex + ga + Sex:geno_snp + locanc_snp, data = temp %>% filter(abs(zAge_AR)<3)) %>% tidy
    m_BMI2 <- lm(zBMI_AR ~ Age_AR + Sex + ga + Sex:geno_snp + locanc_snp, data = temp %>% filter(abs(zBMI_AR)<3)) %>% tidy
    
    m3 <- bind_cols(SNP = df2$SNP[i], tipo = 'Age2', m_Age2)
    m4 <- bind_cols(SNP = df2$SNP[i], tipo = 'BMI2', m_BMI2)
    
    write.table(m3,file=file6,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE)
    write.table(m4,file=file7,append=TRUE,col.names = ifelse(i==1, TRUE, FALSE),row.names=F,quote=FALSE)
    },
    error=function(e){ write.table(df2$POSITION[i],file=file8,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
  
}
