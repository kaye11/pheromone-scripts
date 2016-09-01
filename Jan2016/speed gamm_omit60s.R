#summaries and removing data points with contribution of just 1 track
library(plyr)

pher1 <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/pherspeed.csv", sep=";")
pher <-  subset (pher1, pher1$time > 60, ) #exclude the first 60s


speedsumall <- ddply(pher, c("cond", "time"), summarise,
                     N    = length(V),
                     ID   = length(unique(ID)),
                     mean = mean(V, na.rm=TRUE),
                     sd   = sd(V, na.rm=TRUE),
                     se   = sd / sqrt(N))

speedsumall1 <- ddply(pher1, c("cond", "time"), summarise,
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
s=mgcv:::s


source("AED.R")
source("vif.R")
source("tsDiagGamm.R")
source("summarySE.R")
source("resizewin.R")

resize.win(8,6)

pher$timef=as.factor(pher$timef)
qplot(timef, Vlog, color = cond, data = pher,  geom = "boxplot") + facet_wrap(~cond, scales="free") 

exp=as.data.frame(data.table(cbind(cond=as.numeric(as.factor(pher$cond)), T=pher$time, ID=as.numeric(as.factor(pher$ID)))))

cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#summaries
speedpher.sum <- summarySE(pher, measurevar="V", groupvars=c("cond", "time"))
speedpher.sum.Vlog <- summarySE(pher, measurevar="Vlog", groupvars=c("cond", "time"))
speedpher.sum.condV <- summarySE(pher,measurevar="V", groupvars=c("cond"))

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
#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
pher1 <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cr"), method="REML", data = pher) #best
pher2 <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cs"), method="REML", data = pher)
pher2.1 <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cc"), method="REML", data = pher)

anova(pher0$lme, pher1$lme, pher2$lme, pher2.1$lme) #checking which smoothing spline is the best. you can go further by specifying 
#the number of knots on the formula but I just trust the optimum knot runs by the model

#make random factor and correlations
form <- Vlog~s(time, by=cond, bs="fs", xt="cr")

pher3 <- gamm (form, method="REML",  random=list(ID=~1), data = pher) 
pher4 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), data = pher) #BEST
pher5 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (), data = pher) #same with pher4

anova(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme)

#make variance structures #those with # did not converge
#i normally use varIdent as my data does not vary that much to use exponential or power form of the variance structure,
#see Zuur book for each variance structure formulation
#pher6 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varIdent(form=~1| time), data = pher) #no convergence

#pher7 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varIdent(form=~1| ID), data = pher) #no convergence

pher8 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), 
               weights = varIdent(form=~1| cond), data = pher) 

pher9 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights = varExp(form=~fitted(.)), data = pher) 
#pher9 worked but kind of weirdly

#pher10 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights=varComb(varIdent(form=~1|cond), varIdent (form=~1|ID)), data = pher) 

#pher11 <- gamm (form, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|cond/ID), weights=varComb(varIdent(form=~1|cond), varIdent (form=~1|time)), data = pher) 

anova(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme, pher8$lme)

AIC(pher0$lme, pher1$lme, pher2$lme, pher3$lme, pher4$lme, pher5$lme, pher8$lme, pher9$lme)

with(pher, tsDiagGamm(pher8, timevar=time, observed=Vlog)) #visual checking of residuals

#best model is pher8
summary(pher8$gam)
anova(pher8$gam)


resize.win(6,6)
#checking the model produced by GAMM
op=par(mfrow=c(2,1), mar=c(4.5,4.5,1.5,1.5))
plot(pher8$gam, cex.lab=1.1, cex.axis=1.1, xlab ="Time (s)")

#GAMM has two components: GAM and LME. I normally report the GAM side but it doesn't hurt to look at the LME as well
#most of the time, their trends (sig p value) is the same
summary(pher8$gam)
anova(pher8$gam)

plot(pher8$lme)
anova(pher8$lme)


#plotting (colored and b/w)

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 
source("resizewin.R")

resize.win(8,6)

#colored
ggplot(data=speedsumall1, aes(x=time, y=mean, shape=cond, color=cond)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=20, size=1) +
  scale_colour_manual(values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") +
  labs(list(x = "Time (s)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

#bw

ggplot(data=speedsumall1, aes(x=time, y=mean, shape=cond)) + 
  geom_point(size=5)+ facet_grid(cond~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=20, size=1) +
  scale_shape_discrete (name="Treatment") +
  labs(list(x = "Time (s)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit (0.5, "lines"),
        legend.title=element_blank(),legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

ggplot(data=speedsumall1, aes(x=time, y=mean, shape=cond)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=20, size=1) +
  scale_shape_discrete (name="Treatment") +
  labs(list(x = "Time (s)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit (0.5, "lines"),
        legend.title=element_blank(),legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

SpeedsumCondV <- summarySE(pher, measurevar="V", groupvars=c("cond"))
resize.win (6, 4)
ggplot(data=SpeedsumCondV, aes(x=cond, y=V, fill=cond)) + 
  geom_bar(position = "dodge", stat="identity", color="black") +
  geom_errorbar(width=0.2, size=0.8, aes(ymin=V-se, ymax=V+se))+ 
  scale_fill_manual(values = c(Control="white", Diproline="gray"), name="Treatment")+
  labs(y = "Mean speed (µm/s)") + 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(), 
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        legend.title=element_blank(),legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

#colored
ggplot(data=SpeedsumCondV, aes(x=cond, y=V, fill=cond)) + 
  geom_bar(position = "dodge", stat="identity", color="black") +
  geom_errorbar(width=0.2, size=0.8, aes(ymin=V-se, ymax=V+se))+ 
  scale_fill_manual(values =  c(Control="#f1a340", Diproline="#998ec3"), name="Treatment")+
  labs(y = "Mean speed (µm/s)") + 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(), 
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        legend.title=element_blank(),legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))


#SPEED (Overall)

levene.test(pher$V, group=pher$cond, location="mean") #homogenous variances

wilcox.test(V~cond, data=pher, paired=FALSE) #significantly different

levene.test(pher$Vlog, group=pher$cond, location="mean") #homogenous variances

wilcox.test(Vlog~cond, data=pher, paired=FALSE) #significantly different
