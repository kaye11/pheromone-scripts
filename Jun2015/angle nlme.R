library(ggplot2)
library(mgcv)
library(grid)
library(data.table)

source("summarySE.R")

phersum <- summarySE(pherspeed, measurevar="angs", groupvars=c("cond", "T", "ARF"), na.rm=TRUE)
phersum$cond2 <- factor(phersum$cond, levels=c("HLB", "DPR"), labels=c("Control", "DPR"))

phersum=na.omit(phersum)

##check for collinearity and correlation, this only applies to the explanatory variables!
source("vif.R")
source ("AED.R")
exp=as.data.frame(data.table(cbind(cond2=phersum$cond2, T=phersum$T, ID=phersum$ARF)))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

source("summarySE.R")
dfc2 <- summarySE(phersum, measurevar="angs", groupvars=c("ARF"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(angs~cond2, data=phersum)
boxplot(angs~ARF, data=phersum)
boxplot (angs~T, data=phersum)

#levene test
library(lawstat)
levene.test(phersum$angs, group=phersum$cond2, location="mean") #equal
levene.test(phersum$angs, group=phersum$ARF, location="mean") #unequal
levene.test(phersum$angs, group=phersum$T, location="mean") #equal

#nlme

#lm model
PM.lm <- lm(angs ~cond2*T, data=phersum)

#barplots again
op=par(mfrow=c(2,2))
plot(PM.lm)

PM.E=rstandard(PM.lm)
boxplot(PM.E ~ cond2, data=phersum, axes=TRUE)
abline(0,0); axis=(2)
boxplot(PM.E ~ T, data=phersum, axes=TRUE)
abline(0,0); axis=(2)
boxplot(PM.E ~ ARF, data=phersum, axes=TRUE)
abline(0,0); axis=(2)
par(op)

#fit a gls
Form <- formula (angs ~ cond2*T)
PM.gls<- gls(Form, data=phersum)



#nlme model
PM1.lme <- lme (Form, random = ~1|ARF, method="REML", data=phersum)

#PM2.lme <- lme (Form, random = ~1|ARF,  weights=varIdent(form=~1|ARF), method="REML", data=phersum) 

#PM3.lme <- lme (angs ~ cond2*T, random = ~1|ARF,  weights=varIdent(form=~1|ARF), correlation=corAR1 (form=~1|ARF/cond2), method="REML", data=phersum) 

#PM4.lme <- lme (Form, random = ~1|ARF,  weights=varIdent(form=~1|T), correlation=corAR1 (), method="REML", data=phersum) 

PM5.lme <- lme (Form, random = ~1|ARF,  weights=varIdent(form=~1|cond2), correlation=corAR1 (), method="REML", data=phersum) 

PM6.lme <- lme (angs ~ cond2*T, random = ~1|ARF,  weights=varIdent(form=~1|cond2), correlation=corAR1 (form=~1|ARF/cond2), 
                method="REML", data=phersum) #same with PM5.lme



anova(PM.gls, PM1.lme, PM5.lme, PM6.lme)

summary(PM6.lme)
anova(PM6.lme)

#residuals
PM.E2<-resid(PM6.lme,type="normalized")
PM.F2<-fitted(PM6.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=PM.F2,y=PM.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(PM.E2~cond2,data=phersum, main="cond2",ylab=MyYlab)
plot(x=phersum$T,y=PM.E,main="T",ylab=MyYlab,xlab="T (sec)")
par(op)

library(lattice)
xyplot (PM.E2 ~ T| cond2, data=phersum, ylab="Residuals", xlab="T (sec)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})
