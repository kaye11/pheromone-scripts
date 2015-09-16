library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
library(plyr)
source("vif.R")
source ("AED.R")
source("lang.R")
source("distmmmarySE.R")
source("lang.R")
s=mgcv:::s

write.table (distall, "d:/Karen's/PhD/R program/Pheromone data/distall.csv", 
             sep=";", col.names=T, row.names=F)

distall <- read.csv("D:/Karen's/PhD/R program/Pheromone data/distall.csv", sep=";")


##check for collinearity and correlation, this only applies to the explanatory variables!
exp=as.data.frame(data.table(cbind(cond=distall$cond, T=distall$T)))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(1,2))
boxplot(distmm~cond, data=distall)
boxplot(distmm~time, data=distall)

#levene
library(lawstat)
levene.test(distall$distmm, group=distall$cond, location="mean") #unequal
levene.test(distall$distmm, group=distall$T, location="mean") #unequal

#distall
PA <- gam (distmm~s(T, by=cond, bs="fs"), method="REML", data = distall) #best
PA1 <- gam (distmm~s(T, by=cond, bs="fs", xt="cr"), method="REML", data = distall) 
PA2 <- gam (distmm~s(T, by=cond, bs="fs", xt="cs"), method="REML", data = distall) 
PA3 <- gamm (distmm~s(T, by=cond, bs="fs", xt="cs"), method="REML", weights=varIdent(form=~1| cond), data = distall) 
#PA4 <- gamm (distmm~s(T, by=cond, bs="fs", xt="cs"), method="REML", weights=varIdent(form=~1| time), data = distall) 
PA5 <- gamm (distmm~s(T, by=cond, bs="fs", xt="cs"), method="REML", correlation=corAR1(), data = distall) 



AIC(PA, PA1, PA2)

summary(PA2)
gam.check(PA2)

dev.new(width=6, height=9)
source("resizewin.R")
resize.win(8,10)

op=par(mfrow=c(2,1), mar=c(5.1,6.1,4.1,2.1))
plot(PA1, cex.lab=1.7, cex.axis=1.7, xlab="Time (s)")
