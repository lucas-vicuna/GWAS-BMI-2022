library(graphics)

k = c(2,3,4,5)
cv = c(0.55651,0.54447,0.54281,0.54239)

kcv = cbind(k,cv)

#par(las=2)
par(mar=c(5,6,1,1))

pdf("/data/lucas/bmi/figures/plot_cv_error.pdf", width=4, height=4.5)

plot(kcv,
     pch = 19,
     type = "b",
     ylim = c(0.54,0.56),
     xlim = c(2,6),
     xaxt='n',
     frame.plot = F,
     ann = FALSE, # remove axis titles
     mgp = c(3, 0.5, 0)) # separates y-axis labels from axis.)

#Add ticks to X axis
axis(side=1, at = c(2,3,4,5), labels = F, las=1, cex.axis=0.6, 
     font = 2, tck=-0.03, pos = 0.540)

# Add labels between ticks
axis(side=1, at = c(2,3,4,5), labels = c("2","3","4","5"), 
     las=1, cex.axis=0.8, font = 1, tck= F, lwd=0, mgp = c(3, 0.7, 0))

mtext(side = 1, text = "K", line = 2.2, adj = 0.4, cex = 1.2, las = 1) 

# Add labels between ticks
axis(2, at = c(540,545,550,555), labels = F, 
     las=1, cex.axis=0.5, font = 1, tck= F, lwd=0, mgp = c(3, 0.7, 0))

mtext(side = 2, 
      text = expression("CV error"),
      line = 2.4, adj = 0.5, cex =1.2, las = 0) 

dev.off()
