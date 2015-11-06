#summaries and removing data points with contribution of just 1 track
library(plyr)

pher <- read.csv("D:/Karen's/PhD/R program/pheromone data/Processed_data/pherspeed.csv", sep=";")
pherdata <- na.omit(pher)

angssumall <- ddply(pherdata, c("cond", "time"), summarise,
                     N    = length(angs),
                     ID   = length(unique(ID)),
                     mean = mean(angs, na.rm=TRUE),
                     sd   = sd(angs, na.rm=TRUE),
                     se   = sd / sqrt(N))


library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)

source("AED.R")
source("vif.R")
source("tsDiagGamm.R")
source("summarySE.R")
source("resizewin.R")

resize.win(9,6)

qplot(timef, angs, color = cond, data = pherdata,  geom = "boxplot") + facet_wrap(~cond, scales="free") 

exp=as.data.frame(data.table(cbind(cond=as.numeric(as.factor(pherdata$cond)), T=pherdata$time, ID=as.numeric(as.factor(pherdata$ID)))))

cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)


#boxplots
op=par(mfrow=c(2,2))
boxplot(angs~cond, data=pherdata)
boxplot(angs~ID, data=pherdata)
boxplot (angs~time, data=pherdata)

#levene
library(lawstat)
levene.test(pherdata$angs, group=pherdata$ID, location="mean") #unequal
levene.test(pherdata$angs, group=pherdata$time, location="mean") #unequal
levene.test(pherdata$angs, group=pherdata$cond, location="mean") #unequal


#gamm
pher0 <- gamm (angs~s(time, by=cond, bs="fs"), method="REML", data = pherdata)
pher1 <- gamm (angs~s(time, by=cond, bs="fs", xt="cr"), method="REML", data = pherdata) #best
pher2 <- gamm (angs~s(time, by=cond, bs="fs", xt="cs"), method="REML", data = pherdata) 

anova(pher0$lme, pher1$lme, pher2$lme)

#make random factor and correlations
form <- angs~s(time, by=cond, bs="fs", xt="cr")

pher3 <- gamm (form, method="REML",  random=list(ID=~1), data = pherdata) 
pher4 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), data = pherdata) #BEST
pher5 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (), data = pherdata) #same with pher4

anova(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme)

#make variance structures
#pher6 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varIdent(form=~1| time), data = pherdata) #no convergence

#pher7 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varIdent(form=~1| ID), data = pherdata) #no convergence

pher8 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), 
               weights = varIdent(form=~1| cond), data = pherdata) 

pher9 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varExp(form=~fitted(.)), data = pherdata) 
#pher9 worked but kind of weirdly


anova(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme, pher8$lme)

AIC(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme, pher8$lme, pher9$lme)

with(pherdata, tsDiagGamm(pher8, timevar=time, observed=angs))

#best model is pher8
summary(pher8$gam)
anova(pher8$gam)
gam.check (pher8$gam)

op=par(mfrow=c(2,1), mar=c(5.1,5.1,3.5,2.1))
plot(pher8$gam, cex.lab=1.3, cex.axis=1.3)


summary(pher8$gam)
anova(pher8$gam)

plot(pher8$lme)
anova(pher8$lme)
