#summaries and removing data points with contribution of just 1 track
library(plyr)

pher <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/pherspeed.csv", sep=";")

speedsumall <- ddply(pher, c("cond", "time"), summarise,
                     N    = length(V),
                     ID   = length(unique(ID)),
                     mean = mean(V, na.rm=TRUE),
                     sd   = sd(V, na.rm=TRUE),
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

qplot(timef, Vlog, color = cond, data = pher,  geom = "boxplot") + facet_wrap(~cond, scales="free") 

exp=as.data.frame(data.table(cbind(cond=as.numeric(as.factor(pher$cond)), T=pher$time, ID=as.numeric(as.factor(pher$ID)))))

cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#summaries
speedpher.sum <- summarySE(pher, measurevar="V", groupvars=c("cond", "time"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~cond, data=pher)
boxplot(Vlog~ID, data=pher)
boxplot (Vlog~time, data=pher)

#levene
library(lawstat)
levene.test(pher$Vlog, group=pher$ID, location="mean") #unequal
levene.test(pher$Vlog, group=pher$time, location="mean") #unequal
levene.test(pher$Vlog, group=pher$cond, location="mean") #unequal

#gamm
pher0 <- gamm (Vlog~s(time, by=cond, bs="fs"), method="REML", data = pher)
pher1 <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cr"), method="REML", data = pher) #best
pher2 <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cs"), method="REML", data = pher) 

anova(pher0$lme, pher1$lme, pher2$lme)

#make random factor and correlations
form <- Vlog~s(time, by=cond, bs="fs", xt="cr")

pher3 <- gamm (form, method="REML",  random=list(ID=~1), data = pher) 
pher4 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), data = pher) #BEST
pher5 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (), data = pher) #same with pher4

anova(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme)

#make variance structures
#pher6 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varIdent(form=~1| time), data = pher) #no convergence

#pher7 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varIdent(form=~1| ID), data = pher) #no convergence

pher8 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), 
               weights = varIdent(form=~1| cond), data = pher) 

pher9 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varExp(form=~fitted(.)), data = pher) 
#pher9 worked but kind of weirdly
 

anova(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme, pher8$lme)

AIC(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme, pher8$lme, pher9$lme)

with(pher, tsDiagGamm(pher8, timevar=time, observed=Vlog))

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


#plotting

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 
source("resizewin.R")

resize.win(9,6)

ggplot(data=speedsumall, aes(x=time, y=mean, shape=cond, color=cond)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=20, size=1) +
  scale_colour_manual(values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") +
  labs(list(x = "Time (s)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 