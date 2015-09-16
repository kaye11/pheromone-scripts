source("summarySE.R")

#merging coordinates and data set (pherdist)
pherdata=merge(pherdist, coordinates, by="rep")

#calculate distance between cell position and bead
pherdata$dist=sqrt(((pherdata$bead.X-pherdata$X)^2)+((pherdata$bead.Y-pherdata$Y)^2))

distsum<- summarySE(pherdata, measurevar="dist", groupvars=c("cond", "time"), na.rm=TRUE)

write.table (pherdata, "d:/Karen's/PhD/R program/Pheromone data/pherdatacomp.csv", 
             sep=";", col.names=T, row.names=F)

pher=pherdata

#distance gamm

library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
source("vif.R")
source ("AED.R")
source("lang.R")
source("summarySE.R")
source("lang.R")
s=mgcv:::s

distall <- ddply(pher, c("time", "cond"), summarise,
                 N    = length(dist),
                 mean = mean(dist, na.rm=TRUE),
                 sumdist= sum(dist, na.rm=TRUE), 
                 sd   = sd(dist, na.rm=TRUE),
                 se   = sd / sqrt(N), na.rm=TRUE,
                 angs = mean(angs, na.rm=TRUE))

distall$distmm <- distall$sumdist/1000


##check for collinearity and correlation, this only applies to the explanatory variables!
exp=as.data.frame(data.table(cbind(cond=pher$cond, T=pher$T, A=pher$ARF)))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)
corvif(exp)
pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(2,2))
boxplot(dist~cond, data=pher)
boxplot(dist~ARF, data=pher)
boxplot (dist~T, data=pher)
boxplot (dist~time, data=pher)

#summaries
source("summarySE.R")
distsum <- summarySE(pher, measurevar="dist", groupvars=c("ARF"))
distsumT <- summarySE(pher, measurevar="dist", groupvars=c("time"))
distsumCond <- summarySE(pher, measurevar="dist", groupvars=c("cond"))
distsumCondT <- summarySE(pher, measurevar="dist", groupvars=c("cond", "time"))

#levene test
library(lawstat)
levene.test(pher$dist, group=pher$cond, location="mean") #unequal
levene.test(pher$dist, group=pher$ARF, location="mean") #unequal
levene.test(pher$dist, group=pher$T, location="mean") #unequal
levene.test(pher$dist, group=pher$time, location="mean") #equal

#models
#gam formula
PA <- gamm (dist~s(time, by=cond, bs="fs", xt="cr"), method="REML", data = pher) #best
PA1 <- gamm (dist~s(time, by=cond, bs="fs", xt="cs"), method="REML", data = pher) 

anova(PA$lme, PA1$lme)

f1 <- dist~s(time, by=cond, bs="fs", xt="cr")

#add random structure
PA2 <- gamm (f1, method="REML",  random=list(ARF=~1), data = pher) 

#add correlation structure
PA3 <- gamm (f1, method="REML", random=list(ARF=~1), correlation= corAR1 (form=~1|cond/ARF), data = pher) 
PA4 <- gamm (f1, method="REML", random=list(ARF=~1), correlation= corAR1 (), data = pher) #same with PA3

anova(PA$lme, PA1$lme, PA2$lme, PA3$lme, PA4$lme)

#add variance structure
#PA5 <- gamm (f1, method="REML", random=list(ARF=~1), correlation= corAR1 (form=~1|cond/ARF), weights = varIdent(form=~1| ARF), data = pher) #no convergence

PA6 <- gamm (f1, method="REML", random=list(ARF=~1), correlation= corAR1 (form=~1|cond/ARF), 
             weights = varIdent(form=~1| cond), data = pher) #best

#PA6.1 <- gamm (f1, method="REML", random=list(ARF=~1), correlation= corAR1 (form=~1|cond/ARF), weights = varIdent(form=~1| cond*time), data = pher) #no convergence

#PA7 <- gamm (f1, method="REML", random=list(ARF=~1), correlation= corAR1 (form=~1|cond/ARF), weights = varIdent(form=~1|time), data = pher) #no convergence

#best model is PA6

gam.check (PA6$gam)

ggplot(data = pher, aes(x=T,y=dist, color=cond))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE)


op=par(mfrow=c(2,1))
plot(PA6$gam)
plot(PA6$lme)

#best model is PA6

gam.check (PA6$gam)

ggplot(data = pher, aes(x=T,y=Vlog, color=cond))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE)

op=par(mfrow=c(2,1))
plot(PA6$gam)
plot(PA6$lme)

summary(PA6$gam)
summary(PA6$lme)





