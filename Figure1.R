library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(nlme)
library(lme4)
library(lmerTest)

# setwd('C:/Users/ebarr/Downloads/')

# setwd('C:/Users/ebarr/Onedrive/Documentos/Genetica')
setwd('C:/Users/ebarr/Onedrive/Descargas')


file1 <- paste0("rfmix15.chr",22,".stats.phen.pel_ceu.txt")
df <- fread(file1, sep='\t', header=T, select = 1, col.names = 'CODE')
# df <- fread(file1, sep='\t', header=T)


df1 <- fread('BaseMarzo2020.csv') 
# df1 <- readr::read_csv('BaseMarzo2020.csv') 

colnames(df1) <- c("CODE","date_birth","Sex","vis_date","Weight","Height","BMI","BAZ")

df1 <- df1 %>% inner_join(df)

df2 <- df1 %>% 
  mutate(vis_date = mdy(vis_date),
         date_birth = mdy(date_birth),
         Age = as.numeric(difftime(vis_date, date_birth, units = c("days"))/365)
  ) %>% 
  select(CODE, vis_date, Sex, Weight, Height, BMI, Age) %>% 
  filter(!is.na(BMI)) %>% 
  group_by(CODE) %>% 
  arrange(vis_date) %>% 
  mutate(time = row_number()) %>% 
  ungroup()


df2 %>% 
  ggplot() +
  # aes(x = Age, y = BMI) + 
  aes(x = Age, y = BMI,color=Sex) +
  geom_smooth(se = F) +
  coord_cartesian(xlim=c(3,18))



df2 %>%
  ggplot() +
  aes(x = Age, y = BMI) +
  # geom_smooth(se = F) +
  geom_point(alpha = 1/30) +
  geom_smooth(se=F,lty='dotted', size=1.2,color='red') +
  # geom_path(data = df2 %>% filter(CODE %in% c( 5,46,168,1006,955)) %>% mutate(CODE=factor(CODE)),aes(color=CODE),size = 0.8) +
  # geom_path(data = df2 %>% filter(CODE %in% c( 647,701,905,750)) %>% mutate(CODE=factor(CODE)),aes(color=CODE),size = 0.8) +
  # geom_path(data = df2 %>% filter(CODE %in% c( 19, 140,640, 1103)) %>% mutate(CODE=factor(CODE)),aes(color=CODE),size = 0.8) +
  geom_path(data = df2 %>% filter(CODE %in% c( 216,1175,916,313)) %>% mutate(CODE=factor(CODE)),aes(color=CODE),size = 0.8) +
  
  # geom_smooth(data = df2 %>% filter(CODE==701 & between(Age,2,10)),color='red',size = 1, method='lm',formula = y ~ poly(x, 2)) +
  # geom_path(aes(color=ind, group=CODE))
  coord_cartesian(xlim=c(0,15), ylim=c(15,40))

#46, 272,608,701

## AR <0: 5,46,168,1006,955
## AR >15: 647,701,905,750
## AR between(AR,5.5,7.5): 19, 140,640, 1103



df3 <- df2 %>% mutate(logBMI = log(BMI),
                      Sex1=ifelse(Sex=='Masculino',1,0))

ctrl=lmeControl(opt='optim',maxIter = 50, msMaxIter = 50, tolerance = 1e-6, niterEM = 25,msMaxEval = 200,msTol = 1e-7)

# mod1 <- lme(logBMI ~ Age + I(Age^2), random= ~ 1+ Age +I(Age^2) | CODE, data=df3%>% filter(between(Age, 2,10))   , method = "ML",control=ctrl)
# mod2 <- lme(logBMI ~ Sex + Age + I(Age^2), random= ~ 1+ Age +I(Age^2) | CODE, data=df3%>% filter(between(Age, 2,10))   , method = "ML",control=ctrl)
# mod3 <- lme(logBMI ~ Age + I(Age^2) + I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,3,18)) , method = "ML",control=ctrl)
# mod3 <- lmer(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3)+ ( 1+ Age +I(Age^2) +I(Age^3) | CODE ), data=df4, REML = FALSE)
mod1 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,1,18)), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )
mod2 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,2,18)), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )
mod3 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3)+ Sex:I(Age^4), random= ~ 1+ Age +I(Age^2) +I(Age^3)+I(Age^4) | CODE, data=df3 %>% filter(between(Age,1,18)), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )
mod4 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,3,18)) , method = "ML",control=ctrl)
mod5 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,1,10) & CODE %in% lista), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )
mod6 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,1,10) & CODE %in% n1), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )
mod7 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,1,10) & CODE %in% n2), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )
mod8 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,1,10) & CODE %in% n2), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )


# mod6 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df3 %>% filter(between(Age,2,10) & CODE %in% lista2), control=lmeControl(msMaxIter=500,opt='optim',msVerbose = TRUE) )


df4 <- 
  df3 %>% 
  filter(between(Age,3,18)) %>% 
  mutate(res = resid(mod4, type="pearson"))

df3 %>% 
  filter(Age<5) %>% 
  count(Age)

df3 %>% 
  filter(Age>17.5) %>% 
  count(Age)
# summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

# broom.mixed::glance(mod1)
broom.mixed::glance(mod2)
broom.mixed::glance(mod3)
broom.mixed::glance(mod4)


# broom.mixed::glance(mod3)
broom.mixed::glance(mod4)

# df3$BMI_pred <- exp(predict(mod1,newdata = df3))
# df3$BMI_pred <- exp(predict(mod2,newdata = df3))
# df3$BMI_pred <- exp(predict(mod3,newdata = df3))

mod4 %>% fitted() %>% exp() %>% summary()
# mod4 %>% fitted(.,level=0) %>% exp() %>% summary()
mod4 %>% predict() %>% exp() %>% summary()
mod4 %>% predict(.,newdata=df3 %>% filter(between(Age,3,18))) %>% exp() %>% summary()
mod4 %>% predict(.,newdata=df3 %>% filter(between(Age,3,18))) %>% exp() %>% summary()

fixed.effects(mod4)
fixef(mod4)
2.9750555027+-0.0477977242-0.0679796008*3.00+0.0107127161*3.00^2-0.0003624384*3^3
2.9750555027+-0.0477977242-0.0679796008*3.00+0.0107127161*3.00^2-0.0003624384*3^3-0.0101951124 -0.0158460691*3 -0.0004255192*3^2+  5.569827e-05*3^3
df3 %>% filter(between(Age,3,18)) %>% filter(CODE==874 & time==4)

predict(mod4,newdata=df3 %>% filter(between(Age,3,18)))[1:10]
predict(mod4,newdata=df3 %>% filter(between(Age,3,18)),level=0)[1:10]

predict(mod4,newdata=df3,level=0)[df3$CODE==874 & df3$time==4]
predict(mod4,newdata=df3,level=0)[df3$CODE==874 & df3$time==4]

fitted(mod4)[1:10]
fitted(mod4, level=0)[1:10]
predict(mod4,newdata=df3 %>% filter(CODE %in% c(1,874) & time==4))

predict(mod4,newdata=df3 %>% filter(CODE %in% c(1,874) & time==4))
predict(mod4,newdata=df3 %>% filter(CODE %in% c(1,874) & time==4),level=0)

ranef(mod4) %>% names()
mod4$coefficients$random$CODE[670:680,]

df3 %>% filter(between(Age,3,18)) %>% select(BMI) %>% summary()





df3 %>% filter(CODE %in% sample(CODE,20)) %>% 
  ggplot() +
  aes(x = Age, y = BMI_pred,group=CODE) +
  geom_line() +
  stat_smooth(aes(group=1),se=F) +
  # geom_smooth(se = F) +
  coord_cartesian(xlim=c(0,20), ylim=c(12,35))

n <- ranef(mod1) %>% nrow()
# A <- ranef(mod1) +  matrix(1,n,1) %*% mod1$coef$fixed
# A <- ranef(mod4) +  matrix(1,n,1) %*% mod1$coef$fixed
mod4$coefficients$fixed
mod=mod3

efectos <- function(mod){
  f = fixef(mod)
  B1 = c(f[4],f[3])
  B2 = c(f[6],f[5])
  B3 = c(f[8],f[7])
  
  A = ranef(mod) 
  
  names(A) <- c('b0','b1','b2','b3')
  
  A <- 
    A %>% 
    mutate(
      AR1= (-(B2[1]+b2)+ sqrt((B2[1]+b2)^2-3*(B3[1]+b3)*(B1[1]+b1)))/(3*(B3[1]+b3)),
      AR2= (-(B2[2]+b2)+ sqrt((B2[2]+b2)^2-3*(B3[2]+b3)*(B1[2]+b1)))/(3*(B3[2]+b3)),
      s1 = sign(6*(B3[1]+b3)*AR1+2*(B2[1]+b2)),
      s2 = sign(6*(B3[2]+b3)*AR2+2*(B2[2]+b2))
    )
  
  B <- A %>% as.tbl
  B$CODE = A %>% row.names() %>% as.numeric()
  B <- B %>% left_join(df2 %>% distinct(CODE,Sex))
  B <- B %>% mutate(AR = ifelse(Sex=='Masculino',AR1,AR2), sgn = ifelse(Sex=='Masculino',s1,s2)) %>% select(-AR1,-AR2,-s1,-s2)
  
  return(B)
  
}
B1 = c(mod4$coefficients$fixed[4],mod4$coefficients$fixed[3])
B2 = c(mod4$coefficients$fixed[6],mod4$coefficients$fixed[5])
B3 = c(mod4$coefficients$fixed[8],mod4$coefficients$fixed[7])

A = ranef(mod4) 

names(A) <- c('b0','b1','b2','b3')

A <- 
  A %>% 
  mutate(
    AR1= (-(B2[1]+b2)+ sqrt((B2[1]+b2)^2-3*(B3[1]+b3)*(B1[1]+b1)))/(3*(B3[1]+b3)),
    AR2= (-(B2[2]+b2)+ sqrt((B2[2]+b2)^2-3*(B3[2]+b3)*(B1[2]+b1)))/(3*(B3[2]+b3))
  )
# (-(B2[i]+b2)+ sqrt((B2[i]+b2)^2-3*(B3[i]+b3)(B1[i]+b1)))/(3*(B3[i]+b3))


  

B <- A %>% as.tbl
B$CODE = A %>% row.names() %>% as.numeric()
B <- B %>% left_join(df2 %>% distinct(CODE,Sex))
B <- B %>% mutate(AR = ifelse(Sex=='Masculino',AR1,AR2)) 

summary(B$AR)
B %>% count(AR<3) %>% mutate(prop = n/nrow(B))



df2 %>%
  ggplot() +
  aes(x = Age, y = BMI) +
  geom_point(alpha = 1/30) +
  geom_smooth(se=F,lty='dotted', size=1.2,color='red') +
  geom_path(data = df2 %>% filter(CODE %in% c(2,5,15,10)) %>% mutate(CODE=factor(CODE)),aes(color=CODE),size = 0.8) +
  # geom_path(data = df2 %>% filter(CODE %in% c(1,9,16,17,25)) %>% mutate(CODE=factor(CODE)),aes(color=CODE),size = 0.8) +
  
  coord_cartesian(xlim=c(0,15), ylim=c(15,40))


B %>% filter(between(AR,5,7))
B %>% count(between(AR,4,7))%>% mutate(prop = n/nrow(B))

summary(mod4)


## Residuos

df4 <- 
  df3 %>% 
  filter(between(Age,3,18)) %>% 
  arrange(CODE,vis_date) %>% 
  mutate(res = resid(mod4, type="pearson")) %>% 
  group_by(CODE) %>% 
  mutate(res_max = max(abs(res)),
         res_mean = mean(abs(res))) %>% 
  ungroup()

df4 %>% 
  count(CODE,res_max) %>% 
  arrange(desc(res_max)) %>% 
  filter(CODE==17)
  # arrange(res_max)


df4 %>% 
  count(CODE,res_mean) %>% 
  arrange(desc(res_mean))
# arrange(res_max)

library(lattice)

plot(mod4, resid(., type = "pearson") ~ time | CODE)
bwplot(resid(mod4, type = "p") ~ time | CODE, panel = panel.bwxplot0, data = df4)
qqnorm(mod4, ~resid(.) | time)
qqnorm(mod4, ~ranef(.))


panel.bwxplot0 <- function(x,y, subscripts, ...)
 {
 panel.grid(h = -1)
 panel.stripplot(x, y, col = "grey", ...)
 panel.bwplot(x, y, pch = "|", ...)
}

install.packages('car')
install.packages(('influence.ME'))
# library(car)
library(influence.ME)
influ <- influence(mod3,group = 'CODE')
influence

dfbetas(influ,
        parameters=c(2,3))

plot(influ,
     which="dfbetas",
     parameters=c(2,3),
     xlab="DFbetaS",
     ylab="School ID")

cooks.distance(influ,
               parameter=3, sort=TRUE)

plot(influ, which="cook",
     cutoff=.17, sort=TRUE,
     xlab="Cook´s Distance",
     ylab="School ID")


bwplot(resid(mod4) ~ factor(time), # Fig. 12.2
       panel = panel.bwxplot0,
       ylab = "Residuals", data = df4)


plot(mod4, resid(., type = "pearson") ~ factor(df4$time))

bwplot(resid(mod4, type = "p") ~ factor(time), # Fig. 16.2
       panel = panel.bwxplot0, # User-defined panel (not shown)
       data = df4)
qqnorm(mod4, ~resid(.) | factor(time)) #

qqnorm(mod4, ~ranef(.))


id <- 0.05 # Argument for qnorm()
outliers.idx <-
   within(df4,
        {
              resid.p <- resid(mod4, type = "pearson") # Pearson resids.
              idx <- abs(resid.p) > -qnorm(id/2) # Indicator vector
              })
qnorm(id/2)
outliers <- subset(outliers.idx, idx) # Data with outliers
nrow(outliers) # Number of outliers

outliers$CODE # IDs of outliers

aug.Pred <- # augPred for M16.2
   augPred(mod4,
             primary = ~time, # Primary covariate
             level = 0:1, # Marginal(0) and subj.-spec.(1)
             length.out = 2)
 plot(aug.Pred, layout = c(4, 4, 1), # Fig. 16.5
        key = list(lines = list(lty = c(1,2)),
                     text = list(c("Marginal", "Subject-specific")),
                     columns = 2))
 
 rsd6p <- resid(fm18.6, type = "p")
 keep <- abs(rsd6p) < 3
 
 df4
 
qqDtx <- qqnorm(mod4, plot.it = FALSE)
xp1 <- xyplot(x ~ y, data.frame(qqDtx)) # Draft plot
update(xp1, # Plot updated
           ylab = "Quantiles of standard normal",
           xlab = "Standardized residuals",
           grid = TRUE)

qqnorm(df4$res ) #

rsd6p <- resid(mod4, type = "p")
keep <- abs(rsd6p) < 2
rsd6x <- rsd6p[keep]
qqnorm(rsd6x)

mod5 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df4 %>% filter(abs(res)<3), method = "ML",control=ctrl)
mod6 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df4 %>% filter(abs(res)<2), method = "ML",control=ctrl)
mod7 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df4 %>% filter(res_max<3), method = "ML",control=ctrl)
mod8 <- lme(logBMI ~ Sex + Sex:Age + Sex:I(Age^2) + Sex:I(Age^3), random= ~ 1+ Age +I(Age^2) +I(Age^3) | CODE, data=df4 %>% filter(res_max<2), method = "ML",control=ctrl)


qqnorm(mod4, ~resid(., type='p'))
qqnorm(mod5, ~resid(., type='p'))
qqnorm(mod6, ~resid(., type='p'))
qqnorm(mod7, ~resid(., type='p'))
qqnorm(mod8, ~resid(., type='p'))

a <- 3
b <- 8


(B <- efectos(mod4)) %>% count(between(AR,a,b))%>% mutate(prop = n/nrow(B))
(B <- efectos(mod5)) %>% count(between(AR,a,b))%>% mutate(prop = n/nrow(B))
(B <- efectos(mod6)) %>% count(between(AR,a,b))%>% mutate(prop = n/nrow(B))
(B <- efectos(mod7)) %>% count(between(AR,a,b))%>% mutate(prop = n/nrow(B))
(B <- efectos(mod8)) %>% count(between(AR,a,b))%>% mutate(prop = n/nrow(B))




efectos(mod4) %>% count(between(AR,3,8))
efectos(mod5) %>% count(between(AR,3,8))
efectos(mod6) %>% count(between(AR,3,8))
efectos(mod7) %>% count(between(AR,3,8))
efectos(mod8) %>% count(between(AR,3,8))

efectos(mod4) %>% nrow()
efectos(mod5) %>% nrow()
efectos(mod6) %>% nrow()
efectos(mod7) %>% nrow()
efectos(mod8) %>% nrow()

(B <- efectos(mod4)) %>% count(sgn)%>% mutate(prop = n/nrow(B))
(B <- efectos(mod5)) %>% count(sgn)%>% mutate(prop = n/nrow(B))
(B <- efectos(mod6)) %>% count(sgn)%>% mutate(prop = n/nrow(B))
(B <- efectos(mod7)) %>% count(sgn)%>% mutate(prop = n/nrow(B))
(B <- efectos(mod8)) %>% count(sgn)%>% mutate(prop = n/nrow(B))

B <- efectos(mod4)

B3 <- efectos(mod3)

B$AR %>% summary()
df3$BMI_pred <- exp(predict(mod4,newdata = df3))
df3$BMI_pred3 <- exp(predict(mod3,newdata = df3))

df2$BMI_pred <- exp(predict(mod4,newdata = df2))
df3 %>% 
  filter(between(Age,3,8)) %>% 
  summarise(n_distinct(CODE))

B2 %>% select(CODE,AR2=AR)

df3 <- 
  df3 %>% 
  left_join(B %>% select(CODE,AR))

df3 <- 
  df3 %>% 
  left_join(B3 %>% select(CODE,AR3=AR))

df2 <- 
  df2 %>% 
  left_join(B %>% select(CODE,Age_AR=AR))

df6 <- 
  df2 %>% 
  mutate(Age = Age_AR)


df6 %>% 
  arrange(CODE,time)

df2$BMI_AR <- exp(predict(mod4,newdata = df6))



df3 %>% 
  mutate(AR_cat = case_when(AR<3 ~ '<3',
                            AR<7 ~ '<7',
                            TRUE ~ '>=7'),
         AR2 =factor(round(AR,0))) %>% 
  group_by(CODE,AR2) %>% 
  sample_n(20) %>% 
  ungroup %>% 
  group_by(AR2) %>%
  summarise(n_distinct(CODE))

df3 %>% filter(AR<1) %>% summarise(n_distinct(CODE))

fwrite(
df3 %>% group_by(AR = trunc(AR)) %>% summarise(n_distinct(CODE)),
'distribucion AR.csv', sep=';')

(df3 %>% group_by(AR = trunc(AR)) %>% summarise(n_distinct(CODE)))$`n_distinct(CODE)` %>% sum()

set.seed(1169)
df4 <- 
  df3 %>% 
  mutate(AR = trunc(AR)) %>% 
  filter(between(AR,1,7)) %>% 
  distinct(CODE,AR) %>% 
  group_by(AR) %>% 
  sample_n(5) %>%
  ungroup()

df4 %>% count(AR)
df4 %>% count(AR2)
df5 <- 
  df3 %>% 
  mutate( AR2 =factor(trunc(AR)))
p <- 
    df5 %>% 
     ggplot() +
     aes(x = Age, y = BMI_pred,group=CODE,color=AR2) +
     geom_line(data=df5 %>% filter(CODE %in% df4$CODE)) +
      # geom_point(data=df5 %>% filter(CODE %in% df4$CODE)) +
     stat_smooth(aes(group=1),se=F,color='black',size=1,linetype = "dotted") +
     coord_cartesian(xlim=c(0,18), ylim=c(12,40)) +
     labs(color = 'Age AR', y=expression(hat(BMI)) ) +
      theme(legend.position="top")
p
    
ggsave('AR - curvas ajustadas.png',plot = p, scale=1.1,width=15)
fwrite(df2,'AR.csv',sep=';',dec=',')

df2 %>% summarise(n_distinct(CODE))
df2 %>% group_by(AR = trunc(Age_AR)) %>% summarise(n_distinct(CODE))

 df7 <- as.data.frame(B)
  
set.seed(123)
df7 <- 
  B %>% 
  filter(AR<1) %>% 
  sample_n(2)

df7 <- 
  B %>% 
  filter(AR<1) %>% 
  sample_n(2)


df3 %>% 
  filter(CODE %in% df7$CODE) %>% 
  mutate(CODE=factor(CODE)) %>% 
  ggplot() +
  aes(x = Age, y = BMI,group=CODE,color=CODE) +
  geom_point()+
  geom_path(linetype='dashed')+
  geom_line(aes(y=BMI_pred)) +
  geom_line(aes(y=BMI_pred3),color='green') +
  # geom_point(data=df5 %>% filter(CODE %in% df4$CODE)) +
  # stat_smooth(aes(group=1),se=F,color='black',size=1,linetype = "dotted") +
  coord_cartesian(xlim=c(0,18), ylim=c(12,40)) +
   theme(legend.position="top")
 
plotly::ggplotly()

df7 <- 
  B %>% 
  filter(between(AR,1,3)) %>% 
  sample_n(2)

df3 %>% 
  filter(CODE %in% df7$CODE) %>% 
  mutate(CODE=factor(CODE)) %>% 
  ggplot() +
  aes(x = Age, y = BMI,group=CODE,color=CODE) +
  geom_point()+
  geom_path(linetype='dashed')+
  geom_line(aes(y=BMI_pred)) +
  # geom_point(data=df5 %>% filter(CODE %in% df4$CODE)) +
  # stat_smooth(aes(group=1),se=F,color='black',size=1,linetype = "dotted") +
  theme(legend.position="top")

df7 <- 
  B %>% 
  filter(between(AR,3,5)) %>% 
  sample_n(2)

df3 %>% 
  filter(CODE %in% df7$CODE) %>% 
  mutate(CODE=factor(CODE)) %>% 
  ggplot() +
  aes(x = Age, y = BMI,group=CODE,color=CODE) +
  geom_point()+
  geom_path(linetype='dashed')+
  geom_line(aes(y=BMI_pred)) +
  # geom_point(data=df5 %>% filter(CODE %in% df4$CODE)) +
  # stat_smooth(aes(group=1),se=F,color='black',size=1,linetype = "dotted") +
  coord_cartesian(xlim=c(0,18), ylim=c(12,40)) +
  theme(legend.position="top")

df7 <- 
  B %>% 
  filter(between(AR,5,8)) %>% 
  sample_n(2)

df3 %>% 
  filter(CODE %in% df7$CODE) %>% 
  mutate(CODE=factor(CODE)) %>% 
  ggplot() +
  aes(x = Age, y = BMI,group=CODE,color=CODE) +
  geom_point()+
  geom_path(linetype='dashed')+
  geom_line(aes(y=BMI_pred)) +
  # geom_point(data=df5 %>% filter(CODE %in% df4$CODE)) +
  # stat_smooth(aes(group=1),se=F,color='black',size=1,linetype = "dotted") +
  coord_cartesian(xlim=c(0,18), ylim=c(12,40)) +
  theme(legend.position="top")

df7 <- 
  B %>% 
  filter(AR>8) %>% 
  sample_n(2)

df3 %>% 
  filter(CODE %in% df7$CODE) %>% 
  mutate(CODE=factor(CODE)) %>% 
  ggplot() +
  aes(x = Age, y = BMI,group=CODE,color=CODE) +
  geom_point()+
  geom_path(linetype='dashed')+
  geom_line(aes(y=BMI_pred)) +
  # geom_point(data=df5 %>% filter(CODE %in% df4$CODE)) +
  # stat_smooth(aes(group=1),se=F,color='black',size=1,linetype = "dotted") +
  coord_cartesian(xlim=c(0,18), ylim=c(12,40)) +
  theme(legend.position="top")

df3 %>% filter(CODE==290)
library(plotly)

B %>% filter