library(ggplot2)
library(mgcv)
library(grid)
library(data.table)

Pher$Time=Pher$T*30

Pher$treat2 <- factor(Pher$Treatment, levels=c("HLB", "DPR"), labels=c("Control", "DPR"))
#check plot
qplot(Time,Count, data = Pher, colour = treat2) + facet_grid(~treat2) + 
  geom_smooth(method="loess")

#normalized for area
Pher$CountC=Pher$Count/(pi*Pher$Rad^2)

qplot(Time,CountC, data = Pher, colour = treat2) + facet_grid(~treat2) + 
  geom_smooth(method="loess")

#Normalization of data per treat2
Pher$CountN = NA
AB<-split(Pher, Pher$treat2)
CountNorm <- lapply( AB , function(x) scale(x[,c("CountC")], center=T, scale=T))
Pher$CountN<- unsplit(CountNorm, Pher$treat2)

qplot(Time,CountN, data = Pher, colour = treat2) + facet_grid(~treat2) + 
  geom_smooth(method="loess")

qplot(Time,CountN, color = treat2, data = Pher,  geom = "boxplot") + facet_wrap(~treat2)

#make new variable as ID (treat2, ID2, bead)

Pher$ID2 <- paste(Pher$Replicate, Pher$Bead, sep = "-")


##check for collinearity and correlation, this only applies to the explanatory variables!
source("vif.R")
source ("AED.R")
exp=as.data.frame(data.table(cbind(treat2=Pher$treat2, T=Pher$Time, ID=Pher$ID2)))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

source("summarySE.R")
dfc2 <- summarySE(Pher, measurevar="CountN", groupvars=c("ID2"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(CountN~treat2, data=Pher)
boxplot(CountN~ID2, data=Pher)
boxplot (CountN~Time, data=Pher)

shapiro.test(Pher$Count)
shapiro.test(Pher$CountC)
shapiro.test(Pher$CountN)


#levene test
library(lawstat)
levene.test(Pher$CountN, group=Pher$treat2, location="mean") #equal
levene.test(Pher$CountN, group=Pher$ID2, location="mean") #unequal
levene.test(Pher$CountN, group=Pher$Time, location="mean") #equal

#nlme

#lm model
PM.lm <- lm(CountN ~treat2*Time, data=Pher)

#barplots again
op=par(mfrow=c(2,2))
plot(PM.lm)

PM.E=rstandard(PM.lm)
boxplot(PM.E ~ treat2, data=Pher, axes=TRUE)
abline(0,0); axis=(2)
boxplot(PM.E ~ Time, data=Pher, axes=TRUE)
abline(0,0); axis=(2)
boxplot(PM.E ~ ID2, data=Pher, axes=TRUE)
abline(0,0); axis=(2)
par(op)

#fit a gls
Form <- formula (CountN ~ treat2*Time)
PM.gls<- gls(Form, data=Pher)



#nlme model
PM1.lme <- lme (Form, random = ~1|ID2, method="REML", data=Pher)

#PM2.lme <- lme (Form, random = ~1|ID2,  weights=varIdent(form=~1|ID2), method="REML", data=Pher) 

#PM3.lme <- lme (CountN ~ treat2*Time, random = ~1|ID2,  weights=varIdent(form=~1|ID2), correlation=corAR1 (form=~1|ID2/treat2), method="REML", data=Pher) 

#PM4.lme <- lme (Form, random = ~1|ID2,  weights=varIdent(form=~1|T), correlation=corAR1 (), method="REML", data=Pher) 

PM5.lme <- lme (Form, random = ~1|ID2,  weights=varIdent(form=~1|treat2), correlation=corAR1 (), method="REML", data=Pher) 

PM6.lme <- lme (CountN ~ treat2*Time, random = ~1|ID2,  weights=varIdent(form=~1|treat2), correlation=corAR1 (form=~1|ID2/treat2), 
                method="REML", data=Pher) #same with PM5.lme



anova(PM.gls, PM1.lme, PM5.lme, PM6.lme)

summary(PM6.lme)
anova(PM6.lme)

#residuals
PM.E2<-resid(PM6.lme,type="normalized")
PM.F2<-fitted(PM6.lme)

op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"
plot(x=PM.F2,y=PM.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(PM.E2~treat2,data=Pher, main="treat2",ylab=MyYlab)
plot(x=Pher$Time,y=PM.E,main="Time",ylab=MyYlab,xlab="Time (sec)")
boxplot(PM.E2~Time, xlab="Time", ylab=MyYlab)
par(op)

library(lattice)
xyplot (PM.E2 ~ T| treat2, data=Pher, ylab="Residuals", xlab="Time (sec)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.75,col=1,lwd=2)})


#gamm

PM.gamm <- gamm (CountN ~ s (Time, by=treat2, bs="fs"), random=list(ID2=~1), weights=varIdent(form=~1|treat2), 
                 correlation=corAR1 (form=~1|ID2/treat2), data=Pher)


summary(PM.gamm$gam) #DPR sig
anova(PM.gamm$gam) #DPR sig
plot(PM.gamm$gam) #things look linear, more support for using nlme
plot(PM.gamm$lme) # patterns seen on the residuals
summary(PM.gamm$lme)


##plotting
grid.newpage()
text <- element_text(size = 18, face="bold") #change the size of the axes
theme_set(theme_bw())

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="treat2") { 
    value[value=="DPR"] <- "Diproline"
    value[value=="HLB"]   <- "Control"
  }
  return(value)
}

Pher$cond <- factor(Pher$Treatment, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))

ggplot(data = Pher, aes(x=Time,y=CountN, linetype=cond))+ 
  stat_smooth(method="loess", formula=y~x, size=1.5, color="black") +
  labs (list(x = "Time (s)", y = "Normalized cell count"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size =20, face="bold"),legend.position=c(0.5,0.9), legend.direction="horizontal", 
        legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"), 
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


library(multcomp)
summary(glht(PM6.lme, covariate_average=TRUE))

K <- diag(length(coef(PM6.lme)))[-1,]
rownames(K) <- names(coef(PM6.lme))[-1]
summary(glht(PM6.lme,covariate_average=TRUE, linfct=K))
