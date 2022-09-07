
library(graphics)
library(dplyr)
library(RColorBrewer)

for(i in seq(0.5,15.5,1)) {
  pdf(paste0("/data/lucas/bmi/figures/manhattan_bmi_INT",i,".pdf", sep=""),
      width=5.5, height=4.25)
  X = data1[-log10(data1$p_value) >= 2,]
  X2 = X[!is.na(X$SNP),]

  #create data
  x = X2$POS
  y = -log10(X2$p_value)
  z = data.frame(x,y)
  
  #cut in segments
  my_segments = c(1,48210,100848,145702,186543,223330,
                  266721,300915,333096,359821,390431,
                  420320,448935,470867,490494,509175,
                  528626,546438,564732,577262,592132,
                  600843,609122)
  
  gr <- cut(z$x, my_segments,labels = FALSE, right = T)
  gr[is.na(gr)] <- 0
  
  # create color vector with 1 == red, and 2 == deepskyblue2
  z$color <- ifelse(gr %% 2 == 0, "deepskyblue4", "deepskyblue2")
  
  plot(z$x, z$y, type="p", cex = 0.6, pch = 16,
       col = z$color,
       lwd=0.1,
       frame.plot = F,
       xaxt = 'n', # removes x labels,
       xlim = c(1,609122),
       ylim = c(2, 7),
       las = 2,
       cex.lab=1, # size of axis labels
       ann = FALSE, # remove axis titles
       mgp = c(3, 0.8, 0)) 
  
  # adjust y axis label size
  par(cex.axis= 1.1, tck=-0.03)
  
  mtext(side = 1, text = "Chromosome", line = 1, cex = 1) 
  mtext(side = 2, text = expression(paste("-log"[10]," ",italic("P"),"-value",sep="")),
        line = 1.7, adj = 0.5, cex = 1)  
  mtext(side = 3, text = "A", line = 1, adj = -0.15, padj = 1, las = 1, cex=1.8) # to adjust distance of text from axis
   
  options(scipen = 999) # unables scientific notation por axes
  
  # Add ticks to X axis
  axis(side=1, at = c(1,48210,100848,145702,186543,223330,
                      266721,300915,333096,359821,390431,
                      420320,448935,470867,490494,509175,
                      528626,546438,564732,577262,592132,
                      600843,609122), 
       labels = F, las=1, cex.axis=0.6, font = 2,
       tck=-0.03, pos = 1.85)
  
  # Add labels between ticks
  axis(side=1, at = c(24106,74530,123276,166124,204938,
                      245027,283819,317007,346460,375127,
                      405377,434629,459902,480682,499836,
                      518902,537533,555586,570998,584698,
                      596489,604984), 
       labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22), 
       las=2, cex.axis=0.7, font = 1, tck= F, lwd=0, mgp = c(3, 0.12, 0))
       
  dev.off()
}


