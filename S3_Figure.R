library(dplyr)
library(data.table)
library(graphics)
library(Hmisc)

edu = fread("/path/educ_materna_se.tsv", sep = "\t", header = T)
sig = fread("/path/significancia.tsv", header = T)
sig = data.frame(sig) %>% 
  filter(grepl("MPU",term)) %>%
  filter(p.value<0.05)

edu_eu = subset(edu,MPU==0, c(EDU,INT,fit,lwr,upr))
edu_eu_1 = subset(edu_eu, EDU==123, select = -EDU)
edu_eu_2 = subset(edu_eu, EDU==4, select = -EDU)
edu_eu_3 = subset(edu_eu, EDU==567, select = -EDU)
edu_eu_4 = subset(edu_eu, EDU==8910, select = -EDU)

edu_mpu = subset(edu,MPU==1, c(EDU,INT,fit,lwr,upr))
edu_mpu_1 = subset(edu_mpu, EDU==123, select = -EDU)
edu_mpu_2 = subset(edu_mpu, EDU==4, select = -EDU)
edu_mpu_3 = subset(edu_mpu, EDU==567, select = -EDU)
edu_mpu_4 = subset(edu_mpu, EDU==8910, select = -EDU)

pdf("/path/plot_agear_bmiar_ga_4EDU.pdf", width=5, height=5) 

layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
par(cex=0.6) # size of chart
par(mai=c(0.8,0.6,0.3,0.3)) 
par(oma=c(0.3,0.3,0.3,0.3))

x <- c(1,2,3,4)

for (i in x) {
  edu = cbind(get(paste0("edu_mpu_",i)),get(paste0("edu_eu_",i)))
  colnames(edu)=c("INT","fM","lM","uM","INT2","fE","lE","uE")   
  
  plot(edu$INT,edu$fM,
       type="p", 
       bty="l",
       col = "red",
       lwd=0.001,
       yaxt= 'n', 
       pch = 19,
       xlim =c(1,16),
       ylim=c(1,30),
       xaxt = 'n',
       las=2,
       mgp = c(3, 0.8, 0),
       frame.plot = F,
       cex=1.5,
       ann = FALSE) 
  
  # Error bars:
  errbar(edu$INT, edu$fM, yplus=edu$uM, yminus=edu$lM, cap=.015,  add=T, 
         lty=1, ylim=30, lwd=0.01, errbar.col="red",type="p")
  errbar(edu$INT, edu$fE, yplus=edu$uE, yminus=edu$lE, cap=.015,  add=T, 
         lty=1, ylim=30, lwd=0.01, errbar.col="deepskyblue2",type="p")
  
  lines(edu$INT,edu$fE, lwd=0.01, pch=19 , type="p", bty="l",col = "deepskyblue2", cex=1.5 )
  lines(edu$INT,edu$fM, lwd=0.01, pch=19 , type="p", bty="l",col = "red", cex=1.5 )
  
  #Add ticks to X axis
  axis(side=1, at = seq(1,16,1), 
       labels = c("0.5-1.5","1.5-2.5","2.5-3.5","3.5-4.5","4.5-5.5","5.5-6.5",
                  "6.5-7.5","7.5-8.5","8.5-9.5","9.5-10.5","10.5-11.5","11.5-12.5",
                  "12.5-13.5","13.5-14.5","14.5-15.5","15.5-16.5"), 
       las=1, cex.axis=0.9, font = 1,las =2,
       tck=-0.03, pos = 0,
       mgp = c(3, 0.8, 0))
  axis(2,cex.axis=1, las=2)
  
  text(x=6, y=21, labels = "**", cex = 1.0)
  text(x=7, y=21, labels = "**", cex = 1.0)
  text(x=8, y=22.5, labels = "*", cex = 1.0)
  text(x=9, y=23, labels = "**", cex = 1.0)
  text(x=10, y=24, labels = "**", cex = 1.0)
  text(x=11, y=26, labels = "**", cex = 1.0)
  text(x=12, y=27, labels = "**", cex = 1.0)
  text(x=13, y=28, labels = "**", cex = 1.0)
  text(x=14, y=29.5, labels = "**", cex = 1.0)
  text(x=15, y=30, labels = "**", cex = 1.0)
  text(x=16, y=30.5, labels = "**", cex = 1.0)

  # Legend in 4th plot:
  if (i == 4) {legend(7, 10, legend=c("Mapuche", "European"),
                      col=c("red", "deepskyblue2"), lty=1:2, cex=1.1,
                      bty = "n",pch = 19,lwd=0.01)} 
  # Label axes:
  if (i == 3 | i == 4){mtext("Age (years)", side=1, line=5, cex=1.1)}
  if (i == 1 | i == 3){mtext(expression(paste(plain("BMI (kg/m") ^ plain("2"),plain(")"),sep="")), 
                             side=2, line=2.5, cex=1.1)}
  # Label panels:
  mtext(side = 2, text = LETTERS[i], line = 1, 
        adj = 1.8, padj = -2.8, las = 2, cex=2) 
}

dev.off()
