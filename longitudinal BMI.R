# This script was used to perform the longitudinal GWAS on BMI, including local ancestry and global ancestry. 

library(data.table)
library(dplyr)
library(nlme)
library(lme4)
library(lmerTest)
library(parallel)

args <- commandArgs(TRUE)
chr <- args[1]

# t <- proc.time()
link = '/path/bmi/'
link2 = '/path/'

file1 <- 'adolescents_5.5_15.5_years_bmi.txt'
file2 <- paste0('detail/detail_v2_chr',chr,'.tsv')
file3 <- paste0("rfmix15.chr",chr,".stats.phen.txt")

data1 <- fread(paste0(link,file1), sep=' ',select=c(1,2,5,6),header = T)

data1 <- data1 %>% 
  add_count(CODE) %>% 
  filter(n>2) %>% 
  select(-n)

data2 <- fread(paste0(link,file2), sep='\t', header=T)

data2 <- data2 %>% 
  filter(SIG==1) %>% 
  select(CHR,SNP,POSITION)

data3 <- fread(paste0(link2,file3), 
               sep='\t', 
               select =c(1,6,data2$POSITION, data2$POSITION+1),
               header=T)

colnames(data3)[1:2] <- c('CODE','ga')

data4 <- data1 %>% inner_join(data3)

ctrl=lmeControl(opt='optim',maxIter = 50, msMaxIter = 50, tolerance = 1e-6, niterEM = 25,msMaxEval = 200,msTol = 1e-7)

# k <- nrow(data2)
# 
# t <- proc.time()
# 
# k <- 6
# for(i in 1:k){
#   
#   temp <- data4[,.SD,.SDcols=c(1:5,i+5,i+5+k)]
#   names(temp)[6] <- 'geno_snp'
#   names(temp)[7] <- 'locanc_snp'
#   
#   tryCatch({
#     lme.def <- lme(log(BMI) ~ (Age) + I(Age^2) + I(Age^3) + as.factor(Sex) + as.factor(Sex)*Age+ga+geno_snp+locanc_snp, 
#                    random = ~ Age+I(Age^2)|CODE, data=temp, method = "ML",control=ctrl)
#     myRes=summary(lme.def)$tTable
#     rownamesmyRes=row.names(myRes)
#     row.names(myRes)=c(rownamesmyRes[1:6],paste0('snp.',data2$SNP[i]),'locanc_snp',rownamesmyRes[9])
#     write.table(myRes,file=paste0("chr",chr,"_1SNP_locanc.txt",sep=""),append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
#     error=function(e){ write.table(data2$POSITION[i],file=paste("chr",chr,"_1SNP_locancfailed.txt",sep=""),
#                                    append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })
#   
# }
# 
# proc.time()-t

numCores <- detectCores()
k <- nrow(data2)
i <-  1:k

f <- function(i){
  temp <- data4[,.SD,.SDcols=c(1:5,i+5,i+5+k)]
  names(temp)[6] <- 'geno_snp'
  names(temp)[7] <- 'locanc_snp'

  tryCatch({
    lme.def <- lme(log(BMI) ~ (Age) + I(Age^2) + I(Age^3) + as.factor(Sex) + as.factor(Sex)*Age+ga+geno_snp+locanc_snp,
                   random = ~ Age+I(Age^2)|CODE, data=temp, method = "ML",control=ctrl)
    myRes=summary(lme.def)$tTable
    rownamesmyRes=row.names(myRes)
    row.names(myRes)=c(rownamesmyRes[1:6],paste0('snp.',data2$SNP[i]),'locanc_snp',rownamesmyRes[9])
    write.table(myRes,file=paste0("chr",chr,"_1SNP_locanc.txt",sep=""),append=TRUE,col.names=FALSE,row.names=TRUE,quote=FALSE)},
    error=function(e){ write.table(data2$POSITION[i],file=paste("chr",chr,"_1SNP_locancfailed.txt",sep=""),
                                   append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) })

}

t <- proc.time()
results <- mclapply(i, f, mc.cores = numCores)
proc.time()-t
