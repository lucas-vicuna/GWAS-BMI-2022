# This script was used to perform the cross-sectional GWAS on BMI, including 4 combinations of covariates: 1. With local ancestry, without global ancestry and without the first 10 genetic PCs; 2. Without local ancestry, with global ancestry and without the first 10 genetic PCs; 3. With local ancestry, without global ancestry and with the first 10 genetic PCs; 4. Without local ancestry, without global ancestry and with the first 10 genetic PCs. 

library(data.table)
library(dplyr)
library(nlme)
library(lme4)
library(lmerTest)
library(parallel)
library(readr)

args <- commandArgs(TRUE)
chr <- args[1]

setwd('~/proyectos/')

# f1 <- 'Genetica/data/adolescents_5.5_15.5_years_bmi.txt'
# f2 <- 'Genetica/data/detail_v2_chr22.tsv'
# f3 <- 'Genetica/data/rfmix15.chr22.stats.phen.pel_ceu.txt'
# f12 <- 'Genetica/data/educacion_madre2.csv'
# f13 <- 'Genetica/data/GOCS_904PCA.txt'

f1 <- '/path/adolescents_5.5_15.5_years_bmi.txt'
f2 <- paste0('/path/detail_v2_chr',chr,'.tsv')
f3 <- paste0("/path/rfmix15.chr",chr,".stats.phen.txt")

f12 <- 'data/educacion_madre2.csv'
f13 <- 'data/GOCS_904PCA.txt'

f4i <- paste0('/path/outputi/chr',chr,'.txt')
f4ii <- paste0('/path/outputii/chr',chr,'.txt')
f4iii <- paste0('/path/outputiii/chr',chr,'.txt')
f4iv <- paste0('/path/outputiv/chr',chr,'.txt')

f5i <- paste0('/path/outputi/chr',chr,'.txt')
f5ii <- paste0('/path/outputii/chr',chr,'.txt')
f5iii <- paste0('/path/outputiii/chr',chr,'.txt')
f5iv <- paste0('/path/outputiv/chr',chr,'.txt')


data1 <- fread(f1, sep=' ',select=c(1,2,5,6),header = T)

data1 <- data1 %>% 
  add_count(CODE) %>% 
  filter(n>2) %>% 
  select(-n)

df12 <- read_csv(f12)
df13 <- read_delim(f13, delim='\t',col_names = F)
names(df13) <- c('CODE',paste0('X',1:10))

data1 <- data1 %>% inner_join(df12) %>% inner_join(df13) %>% mutate(Sex=factor(Sex), edu=factor(edu))

data2 <- fread(f2, sep='\t', header=T)
data2 <- data2 %>% 
  filter(SIG==1) %>% 
  select(CHR,SNP,POSITION)

data3 <- fread(f3, 
               sep='\t', 
               select =c(1,6,data2$POSITION, data2$POSITION+1),
               header=T)

colnames(data3)[1:2] <- c('CODE','ga')

data4 <- data1 %>% inner_join(data3)

ctrl=lmeControl(opt='optim',maxIter = 50, msMaxIter = 50, tolerance = 1e-6, niterEM = 25,msMaxEval = 200,msTol = 1e-7)


k <- nrow(data2)

Myoutput <- function(x, i){
  
  a = broom.mixed::tidy(x) %>% filter(effect=='fixed') %>% select(-c(effect,group))
  a = bind_cols(CHR = chr, SNP= data2$SNP[i], a)
  
  return(a)
}

f=function(x){
  
  if(x==1){
    
    #i With local ancestry, without global ancestry, without first 10 PCs, with maternal education level
    
    for(i in 1:k){
      tempi <- data4[,.SD,.SDcols=c(1:16,i+16,i+16+k)]
      names(tempi)[17] <- 'geno_snp'
      names(tempi)[18] <- 'locanc_snp'
      
      tryCatch({
        mod_i <- lme(log(BMI) ~ (Age) + I(Age^2) + I(Age^3) + Sex + Sex:Age+ geno_snp+locanc_snp+edu,
                     random = ~ Age+I(Age^2)|CODE, data=tempi, method = "ML",control=ctrl)
        
        output_i = Myoutput(mod_i, i)
        write.table(output_i,file=f4i,append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
        error=function(e){ write.table(data2$POSITION[i],file=f5i,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
    }
  }
  

  if(x==2){
    #ii Without local ancestry, with global ancestry, without first 10 PCs, with maternal education level
    
    for(i in 1:k){
      tempii <- data4[,.SD,.SDcols=c(1:16,i+16,i+16+k)]
      names(tempii)[17] <- 'geno_snp'
      names(tempii)[18] <- 'locanc_snp'    
      tryCatch({
        mod_ii <- lme(log(BMI) ~ (Age) + I(Age^2) + I(Age^3) + Sex + Sex:Age+ga+geno_snp+edu,
                      random = ~ Age+I(Age^2)|CODE, data=tempii, method = "ML",control=ctrl)
        
        output_ii = Myoutput(mod_ii, i)
        write.table(output_ii,file=f4ii,append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
        error=function(e){ write.table(data2$POSITION[i],file=f5ii,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
    }
  }

  if(x==3){
    #iii With local ancestry, without global ancestry, with first 10 PCs, with maternal education level
    
    for(i in 1:k){
      tempiii <- data4[,.SD,.SDcols=c(1:16,i+16,i+16+k)]
      names(tempiii)[17] <- 'geno_snp'
      names(tempiii)[18] <- 'locanc_snp'
      tryCatch({
        mod_iii <- lme(log(BMI) ~ (Age) + I(Age^2) + I(Age^3) + Sex + Sex:Age+geno_snp+locanc_snp+edu+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                       random = ~ Age+I(Age^2)|CODE, data=tempiii, method = "ML",control=ctrl)
        
        output_iii = Myoutput(mod_iii, i)
        write.table(output_iii,file=f4iii,append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
        error=function(e){ write.table(data2$POSITION[i],file=f5iii,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
    }
  }

  if(x==4){
    #iv Without local ancestry, without global ancestry, with first 10 PCs, with maternal education level
    
    for(i in 1:k){
      tempiv <- data4[,.SD,.SDcols=c(1:16,i+16,i+16+k)]
      names(tempiv)[17] <- 'geno_snp'
      names(tempiv)[18] <- 'locanc_snp'
      
      tryCatch({
        mod_iv <- lme(log(BMI) ~ (Age) + I(Age^2) + I(Age^3) + Sex + Sex:Age+geno_snp+edu+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                      random = ~ Age+I(Age^2)|CODE, data=tempiv, method = "ML",control=ctrl)
        
        output_iv = Myoutput(mod_iv, i)
        write.table(output_iv,file=f4iv,append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
        error=function(e){ write.table(data2$POSITION[i],file=f5iv,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
      
    }
  }
}

#t <- proc.time()
results <- mclapply(1:4, f, mc.cores = 4)
#proc.time()-t

