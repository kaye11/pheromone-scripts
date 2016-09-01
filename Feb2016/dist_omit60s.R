distall1 <- read.csv("d:/Karen's/PhD/R program/Pheromone data/Processed_data/distall.csv", sep=";")
distall <- subset(distall1, distall1$time> 60, )

library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(data.table)
library(plyr)
library(mgcv)
source("vif.R")
source ("AED.R")
source("lang.R")
source("summarySE.R")
source("lang.R")
source("tsDiagGamm.R")
source("resizewin.R")
resize.win(8,6)


qplot(x=time, y=sum, data=distall)+geom_line()+facet_grid(~cond, scales="free")

distall$sum_mm <- distall$sum/1000
distall1$sum_mm <- distall1$sum/1000

qplot(x=time, y=sum_mm, data=distall)+geom_line()+facet_grid(~cond, scales="free")


##check for collinearity and correlation, this only applies to the explanatory variables!
exp=as.data.frame(data.table(cbind(cond=distall$cond, T=distall$time)))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)


pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(1,2))
boxplot(sum_mm~cond, data=distall)
boxplot(sum_mm~time, data=distall)

#levene
library(lawstat)
levene.test(distall$sum_mm, group=distall$cond, location="mean") #unequal
levene.test(distall$sum_mm, group=distall$time, location="mean") #unequal

#try linear model

pherdist.gls <- gls (sum_mm~time*cond, method="REML", data=distall)
pherdist.lm <- lm (sum_mm~time*cond,  data=distall) 
#pherdist.lme <- lme (sum_mm~time*cond, method="REML", correlation=corAR1(), data=distall ) #lme doesn't work
#residuals of pherdist.lm shows nonlinearity but normal qq plots so choose lm

pherdist.gls1 <- gls (sum_mm~time*cond, method="REML", correlation=corAR1(), data=distall)
pherdist.gls2 <- gls (sum_mm~time*cond, method="REML", correlation= corAR1 (form=~1|cond), data=distall)

anova(pherdist.gls, pherdist.lm, pherdist.gls1, pherdist.gls2)
#pherdist2 produced lowest AIC values but there is a problem with the residuals (biased and heteroscedastic)


summary(pherdist.lm)
anova(pherdist.lm)

#residuals
distall.E2<-resid(pherdist.lm,type="response")
distall.F2<-fitted(pherdist.lm)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=distall.F2,y=distall.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(distall.E2~cond,data=distall, main="cond",ylab=MyYlab)
plot(x=distall$time,y=distall.E2,main="Time",ylab=MyYlab,xlab="Time (sec)")
par(op)

xyplot (distall.E2 ~ time| cond, data=distall, ylab="Residuals", xlab="Time (sec)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

library(lsmeans)
lsmeans(pherdist.lm, pairwise~cond*time, adjust="tukey")

#let's plot this!

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 



#pherdist fit

pherdist.fit <- as.data.frame(predict(pherdist.lm, distall, se.fit = TRUE, level = 0))

pherdist.fit$upr <- pherdist.fit$fit + (1.96 * pherdist.fit$se)
pherdist.fit$lwr <- pherdist.fit$fit - (1.96 * pherdist.fit$se)

pherdist.fit.combdata <- cbind(distall, pherdist.fit)


ggplot(data=distall1, aes(x=time, y=sum_mm, shape=cond, color=cond)) + geom_point(size=5)+ 
  geom_smooth(data=pherdist.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=cond), stat="identity", alpha=0.2)+ 
  scale_colour_manual(values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_manual(values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment") + 
  labs(list(x = "Time (s)", y = "Sum distance \n from the bead (mm)"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom", legend.direction="horizontal",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

#bw

ggplot(data=distall1, aes(x=time, y=sum_mm, shape=cond)) + geom_point(size=5)+ 
  geom_smooth(data=pherdist.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, linetype=cond), 
              color="black", stat="identity", alpha=0.2)+ 
  scale_shape_discrete (name="Treatment") +
  scale_linetype_manual(values = c(Control="dashed", Diproline="solid"), name="Treatment")+
  labs(list(x = "Time (s)", y = "Sum distance \n from the bead (mm)"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom", legend.direction="horizontal",
        legend.title=element_blank(),legend.key.width=unit(1.8,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 


