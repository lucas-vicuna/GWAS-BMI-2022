library(dplyr)
librady(data.table)
library(graphics)
library(Hmisc)

# Format data:

# BMI x Age data:

edu = fread("/data/lucas/bmi/data_esteban/age_strata/educ_materna_se.tsv", header = T)
edu_eu = subset(edu,MPU==0, c(EDU,INT,fit,lwr,upr))
edu_eu_1 = subset(edu_eu, EDU==123, select = -EDU)

edu_mpu = subset(edu,MPU==1, c(EDU,INT,fit,lwr,upr))
edu_mpu_1 = subset(edu_mpu, EDU==123, select = -EDU)

edu = cbind(edu_mpu_1,edu_eu_1)
colnames(edu)=c("INT","fM","lM","uM","INT2","fE","lE","uE")  

# Age-AR & BMI-AR data:

data = read.table("/data/lucas/bmi/data_esteban/ageAR/data_AR_trans.tsv", header = T)
colnames(data)[1]="ID"

# convert classes from character to numeric:
data[,"ID"] <- sapply(data[,"ID"], as.numeric)

ga = read.table("/data/lucas/height/glob_anc/ancestria_global_gocs.csv", header = T,sep=",", stringsAsFactors = F)

# convert classes from character to numeric:
ga[,"ID"] <- sapply(ga[,"ID"], as.numeric)

data_ga = left_join(data,ga,"ID")

#####################

# PLOTS

pdf("/data/lucas/bmi/figures/ga_bmi_ageAR_bmiAR_mapuche_european2.pdf",  
    width=6, height=6)

layout(matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE))
par(cex=0.6) # size of chart
par(mai=c(0.8,0.6,0.3,0.3)) # adjusts thickness of margins (white space)
par(oma=c(0.3,0.3,0.3,0.3)) # adjusts margins of figure (subpanel)

# Plot 1

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

errbar(edu$INT, edu$fE, yplus=edu$uM, yminus=edu$lM, cap=.015,  add=T, 
       lty=1, ylim=36, lwd=0.001, errbar.col="red",type="p")
errbar(edu$INT, edu$fM, yplus=edu$uE, yminus=edu$lE, cap=.015,  add=T, 
       lty=1, ylim=36, lwd=0.001, errbar.col="deepskyblue2",type="p")

lines(edu$INT,edu$fE, lwd=0.001, pch=19 , type="p", 
      bty="l",col = "deepskyblue2", cex=1.5 )
lines(edu$INT,edu$fM, lwd=0.001, pch=19 , type="p", 
      bty="l",col = "red", cex=1.5 )

#Add ticks to X axis
axis(side=1, at = seq(1,16,1), 
     labels = c("0.5-1.5","1.5-2.5","2.5-3.5",
                "3.5-4.5","4.5-5.5","5.5-6.5",
                "6.5-7.5","7.5-8.5","8.5-9.5",
                "9.5-10.5","10.5-11.5","11.5-12.5",
                "12.5-13.5","13.5-14.5","14.5-15.5",
                "15.5-16.5"), 
     las=1, cex.axis=0.9, font = 1,las =2,
     tck=-0.03, pos = 0,
     mgp = c(3, 0.8, 0))

text(x=6, y=21, labels = "**", cex = 1.3)
text(x=7, y=21, labels = "**", cex = 1.3)
text(x=8, y=22.5, labels = "*", cex = 1.3)
text(x=9, y=23, labels = "**", cex = 1.3)
text(x=10, y=24, labels = "**", cex = 1.3)
text(x=11, y=26, labels = "**", cex = 1.3)
text(x=12, y=27, labels = "**", cex = 1.3)
text(x=13, y=28, labels = "**", cex = 1.3)
text(x=14, y=29.5, labels = "**", cex = 1.3)
text(x=15, y=30, labels = "**", cex = 1.3)
text(x=16, y=30.5, labels = "**", cex = 1.3)

axis(2,cex.axis=1, las=2)
mtext(expression(paste(plain("BMI (kg/m") ^ plain("2"),plain(")"),sep="")), 
      side=2, line=2.5, cex=1.1)
mtext("Age (years)", side=1, line=5, cex=1.1)

legend(12, 10, legend=c("Mapuche", "European"),
       col=c("red", "deepskyblue2"), lty=1:2, cex=1.5,
       bty = "n",pch = 19,lwd=0.01)

mtext(side = 2, text = "A", line = 1, adj = 2, padj = -3.3, las = 2, cex=2)  # to adjust distance of text from axis

# Plot 2:

plot(data_ga$MPU,data_ga$Age_AR,
     xlab = "", 
     xaxt='n',
     yaxt='n',
     pch = 19,
     cex=0.7, col="steelblue2", fg=NULL,
     ann = FALSE) # remove axis titles

axis(2,cex.axis=1.2)
mtext("Age-AR (years)", side=2, line=2.4, cex=1.1)
axis(1,cex.axis=1.2)
mtext("Mapuche ancestry proportion", side=1, line=2.8, cex=1.1)

abline(lm(Age_AR ~ MPU, data = data_ga), col = "red")

mtext(side = 2, text = "B", line = 1, adj = 2, padj = -3.3, las = 2, cex=2)

# Plot 3:

plot(data_ga$MPU,data_ga$BMI_AR,
     xlab = "", 
     ylab = "",
     xaxt='n',
     yaxt='n',
     pch = 19,
     cex=0.7, col="steelblue2", fg=NULL,
     ann = FALSE) # remove axis titles

axis(2,cex.axis=1.2)
mtext(expression(paste(plain("BMI-AR (kg/m") ^ plain("2"),plain(")"),sep="")), side=2, line=2.4, cex=1.1)

axis(1,cex.axis=1.2)
mtext("Mapuche ancestry proportion", side=1, line=2.8, cex=1.1)

abline(lm(BMI_AR ~ MPU, data = data_ga), col = "red")

mtext(side = 2, text = "C", line = 1, adj = 2, padj = -3.3, las = 2, cex=2)

dev.off()





