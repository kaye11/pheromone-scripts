#summaries and removing data points with contribution of just 1 track
library(plyr)

pher <- read.csv("D:/Karen's/PhD/R program/pheromone data/Processed_data/pherspeed.csv", sep=";")
pher1 <- na.omit(pher)
pher2 <- subset (pher1, pher1$time > 60, ) #exclude the first 60s


angssumall <- ddply(pher2, c("cond", "time"), summarise,
                     N    = length(angs),
                     ID   = length(unique(ID)),
                     mean = mean(angs, na.rm=TRUE),
                     sd   = sd(angs, na.rm=TRUE),
                     se   = sd / sqrt(N))


angssumall1 <- ddply(pher1, c("cond", "time"), summarise,
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

pher2$timef <- as.factor(pher2$time)

qplot(timef, angs, color = cond, data = pher2,  geom = "boxplot") + facet_wrap(~cond, scales="free") 

exp=as.data.frame(data.table(cbind(cond=as.numeric(as.factor(pher2$cond)), T=pher2$time, ID=as.numeric(as.factor(pher2$ID)))))

cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)


#boxplots
op=par(mfrow=c(2,2))
boxplot(angs~cond, data=pher2)
boxplot(angs~ID, data=pher2)
boxplot (angs~time, data=pher2)

#levene
library(lawstat)
levene.test(pher2$angs, group=pher2$ID, location="mean") #unequal
levene.test(pher2$angs, group=pher2$time, location="mean") #unequal
levene.test(pher2$angs, group=pher2$cond, location="mean") #unequal


#fit a gls
Form <- formula (angs ~ cond*time)
pher.gls<- gls(Form, pher2, na.action=na.omit)


#nlme model
pher1.lme <- lme (Form, random = ~1|ID, method="REML", pher2, na.action=na.omit)

pher2.lme <- lme (Form, random = ~1|ID, method="REML", correlation=corAR1(), pher2, na.action=na.omit)

#pher3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", pher2, na.action=na.omit)

#pher4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), correlation=corAR1 (form=~1|ID/cond), method="REML", pher2, na.action=na.omit)

#pher5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|time), correlation=corAR1 (form=~1|ID/cond), method="REML", pher2, na.action=na.omit) 

pher6.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|cond), 
                  correlation=corAR1 (), method="REML", data=pher2, na.action=na.omit) 

#pher8.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|cond), varIdent (form=~1|ID)),  correlation=corAR1 (), method="REML", data=pher2, na.action=na.omit) 

#pher9.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|cond), varIdent (form=~1|time)),  correlation=corAR1 (), method="REML", data=pher2, na.action=na.omit) 


pher7.lme <- lme (Form, random = ~1|ID,  weights=varExp(form=~fitted(.)), 
                  correlation=corAR1 (), method="REML", data=pher2, na.action=na.omit)


anova(pher.gls, pher1.lme, pher2.lme, pher6.lme, pher7.lme)

#choose pher6.lme
summary(pher6.lme)
anova(pher6.lme)


#residuals
pher.E2<-resid(pher6.lme,type="normalized")
pher.F2<-fitted(pher6.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=pher.F2,y=pher.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(pher.E2~cond,data=pher2, main="Treatment",ylab=MyYlab)
plot(x=pher2$time,y=pher.E2,main="Time",ylab=MyYlab,xlab="Time (sec)")
par(op)

xyplot (pher.E2 ~ time| cond, data=pher, ylab="Residuals", xlab="Time (sec)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#let's plot this!

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 
library (AICcmodavg)


#pherangs fit

pherangs.fit <- as.data.frame(predictSE.lme(pher6.lme, pher2, se.fit = TRUE, level = 0,
                                            print.matrix = FALSE))

pherangs.fit$upr <- pherangs.fit$fit + (1.96 * pherangs.fit$se)
pherangs.fit$lwr <- pherangs.fit$fit - (1.96 * pherangs.fit$se)

pherangs.fit.combdata <- cbind(pher2, pherangs.fit)


ggplot(data=angssumall1, aes(x=time, y=mean, shape=cond, color=cond)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=15, size=1) + 
  geom_smooth(data=pherangs.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=cond), method="lm", stat="identity", alpha=0.2)+ 
  scale_colour_manual(values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment") +
  scale_shape_discrete (name="Treatment") + geom_hline(yintercept=0, size=1)+ 
  scale_fill_manual (values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment")+ 
  labs(list(x = "Time (s)", y = "Mean sine angle"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=element_blank(), 
        legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

#bw

ggplot(data=angssumall1, aes(x=time, y=mean, shape=cond)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=20, size=1) + 
  geom_smooth(data=pherangs.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), 
              color="black", method="lm", stat="identity", alpha=0.2)+ facet_grid(~cond)+
  scale_shape_discrete (name="Treatment") + geom_hline(yintercept=0, linetype="dotted", size=1)+ 
  scale_linetype_manual(values = c(Control="dashed", Diproline="solid"), name="Treatment")+
  labs(list(x = "Time (s)", y = "Mean sine angle"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = element_text(size=20), strip.text.y = text, legend.title=text, legend.text=text, 
        legend.title=element_blank(),legend.key.width=unit(2,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.margin=unit (0.5, "lines"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

