library(readr)
library(dplyr)
library(stringr)

setwd('~/proyectos/')
# df <- read_tsv('bmi/BMIxAge/data/TramosxEdad.csv')
df <- read_tsv('bmi/BMIxAge/data/TramosxEdad.csv')


f1 <- paste0('bmi/BMIxAge/nuevo_analisis/qqplot/tramo 9-16_v2.pdf')


# k=1
pdf(file=f1, width=5, height=10)
# par(mfrow=c(4,2), mar=c(3,3,3,3))
par(mfrow=c(4,2),oma=c(2,2,1,1), mar=c(3,3,3,3))
for(i in 9:16){
  xlab=''
  ylab=''
  df1 <- NULL
  titulo <- NULL
  df1 <- df %>% filter(orden==i)
  intervalo <- df1$INT %>% unique
  intervalo <- intervalo %>% str_replace_all('\\(|\\]','') %>% str_replace(',','-')
  titulo <- paste(intervalo, 'years old')
  # f1 <- paste0('bmi/BMIxAge/nuevo_analisis/qqplot/tramo',i,'.pdf')
  # pdf(f1) 
    
  # if(k %in% c(1,3,5,7)){
  #   ylab='Sample Quantile'
  # }
  # if(k %in% c(7, 8)){
  #   xlab='Theoretical Quantile'
  # }
  qqnorm(df1$zBMI, main = titulo, xlab=xlab, ylab=ylab);qqline(df1$zBMI, col = 2)
  # k=k+1
  # dev.off()
}
mtext('Theoretical Quantile',side=1,line=0,outer=TRUE,cex=1.3)
mtext('Sample Quantile',side=2,line=0,outer=TRUE,cex=1.3,las=0)
dev.off()

