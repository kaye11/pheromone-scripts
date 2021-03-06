pher <- read.csv("D:/Karen's/PhD/R program/pheromone data/Processed_data/pherspeed.csv", sep=";")
coordinates <- read.csv("D:/Karen's/PhD/R program/Pheromone data/coordinates.csv", sep=";")

#merging coordinates and data set (pherdist)
pherdata=merge(pher, coordinates, by="rep")

pherdata$dist=sqrt(((pherdata$bead.X-pherdata$X)^2)+((pherdata$bead.Y-pherdata$Y)^2))

library(plyr)
distall <- ddply(pherdata, c("time", "cond"), summarise,
                 N    = length(dist),
                 mean = mean(dist, na.rm=TRUE),
                 sum= sum(dist, na.rm=TRUE), 
                 sd   = sd(dist, na.rm=TRUE),
                 se   = sd / sqrt(N))

write.table (pherdata, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/pherdatadist.csv", 
             sep=";", col.names=T, row.names=F)

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
resize.win(9,6)


qplot(x=time, y=sum, data=distall)+geom_line()+facet_grid(~cond, scales="free")

distall$sum_mm <- distall$sum/1000

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

#gamm
pherdist <- gam (sum_mm~s(time, by=cond, bs="fs"), method="REML", data = distall)
pherdist1 <- gam (sum_mm~s(time, by=cond, bs="fs", xt="cr"), method="REML", data = distall) 
pherdist2 <- gam (sum_mm~s(time, by=cond, bs="fs", xt="cs", k=5), method="REML", data = distall) 
pherdist3 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", weights=varIdent(form=~1| cond), data = distall) 
#pherdist4 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", weights=varIdent(form=~1| time), data = distall) 
pherdist5 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", correlation=corAR1(), data = distall) 
pherdist6 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", correlation= corAR1 (form=~1|cond), data = distall) 
pherdist7 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", correlation= corAR1 (form=~1|time), data = distall) 
#pherdist8 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", correlation= corAR1 (form=~1|cond), weights= varIdent(form=~1| time), data = distall) #best
pherdist9 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", correlation= corAR1 (form=~1|cond), 
                   weights= varIdent(form=~1|cond), data = distall) #best
pherdist10 <- gamm (sum_mm~s(time, by=cond, bs="fs", xt="cs"), method="REML", correlation= corAR1 (form=~1|cond), 
                    weights= varIdent(form=~1|cond), data = distall) #same with pherdist9

anova(pherdist3$lme, pherdist5$lme, pherdist6$lme, pherdist7$lme, pherdist9$lme, pherdist10$lme)
#best is pherdist6

with(distall, tsDiagGamm(pherdist6, timevar=time, observed=sum_mm))

#best model is pherdist6
summary(pherdist6$gam)
anova(pherdist6$gam)
gam.check (pherdist6$gam)

resize.win(6,6)
op=par(mfrow=c(2,1), mar=c(4.5,4.5,1.5,1.5))
plot(pherdist6$gam, cex.lab=1.1, cex.axis=1.1, xlab ="Time (s)")


plot(pherdist6$lme)
anova(pherdist6$lme)

#lme
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


ggplot(data=distall, aes(x=time, y=sum_mm, shape=cond, color=cond)) + geom_point(size=5)+ 
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

ggplot(data=distall, aes(x=time, y=sum_mm, shape=cond)) + geom_point(size=5)+ 
  geom_smooth(data=pherdist.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, linetype=cond), 
              color="black", stat="identity", alpha=0.2)+ 
  scale_shape_discrete (name="Treatment") +
  scale_linetype_manual(values = c(Control="dashed", Diproline="solid"), name="Treatment")+
  labs(list(x = "Time (s)", y = "Sum distance from the bead (mm)"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom", legend.direction="horizontal",
        legend.title=element_blank(),legend.key.width=unit(1.8,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 


ggplot(data=distall, aes(x=time, y=sum_mm, shape=cond)) + geom_point(size=5)+ 
  scale_shape_discrete (name="Treatment") +
  labs(list(x = "Time (s)", y = "Sum distance from the bead (mm)"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom", legend.direction="horizontal",
        legend.title=element_blank(),legend.key.width=unit(1.8,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 
