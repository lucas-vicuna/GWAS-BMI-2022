# This R script was used to perform the main cross-sectional GWAS on BMI, including local ancestry and global ancestry, but excluding the first 10 PCs. 

library(data.table)
library(dplyr)
library(parallel)
library(purrr)  
library(broom)

args <- commandArgs(TRUE)
chr <- args[1]

link <- '~/proyectos/bmi/'
setwd(link)

file <- paste0('path/TramosxEdad.csv')
file1 <- paste0('path/detail_chr',chr,'.tsv')
file2 <- paste0('path/chr',chr,'_categorias_SNP.tsv')
file3 <- paste0('/path/rfmix15.chr',chr,'.stats.phen.txt')
file4 <- paste0('path/chr',chr,'_BMIxAge.txt')
file5 <- paste0('path/chr',chr,'_BMIxAgeII.txt')

df <- fread(file, sep='\t', header=T)
df <- df %>% select(CODE,INT,orden,zBMI,Sex,Age)

df1 <- fread(file1,header = T,sep='\t')
df2 <- fread(file2,header = T,sep='\t')

df1 <- df1 %>% left_join(df2) %>% select(-OK,-FAILED)

df2 <- df1 %>% 
  filter(VAR0==0) %>%  
  rowwise() %>% 
  mutate(min = min(c(freq0,freq1,freq2)),
         median = median(c(freq0,freq1,freq2))) %>% 
  ungroup() %>% 
  filter((min==0 & median >= 5/904) | (min>=5/904)) %>% 
  select(-min,-median)

output <- function(results,i){
  A <- data.frame()
  for (j in results %>% names())
    {
    m <- (df %>% count(orden,INT) %>% filter(orden == j))$INT
    B <- bind_cols(data.frame(INT = m, orden = j, SNP = df2$SNP[i]),
                   results[j])
    A <- bind_rows(A,B)
  }
  return(A)
}

df3 <- fread(file3,
             sep='\t', 
             select =c(1,6,df2$POSITION, df2$POSITION+1),
             header=T)

colnames(df3)[1:2] <- c('CODE','ga')

df4 <- df %>% inner_join(df3)

k <- df2 %>% nrow()

for (i in 1:k){
  
  temp <- df4[,.SD,.SDcols=c(1:7,i+7,i+7+k)]
  names(temp)[8] <- 'geno_snp'
  names(temp)[9] <- 'locanc_snp'
  
  tryCatch({
    
    results <- 
      temp %>% 
      split(.$orden) %>% 
      map( ~lm(zBMI ~ Sex + ga + geno_snp + locanc_snp, data = .)) %>% 
      map(tidy)  
    
    write.table(output(results,i),file=file4,append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
    error=function(e){ write.table(df2$POSITION[i],file=file5,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
  
}
