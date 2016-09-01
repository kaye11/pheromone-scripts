pher <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Raw data_counts/Phercount.csv", sep=";")

library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
library(nlme)

source("vif.R")
source("AED.R")
source("lang.R")

dev.new(width=6, height=9)
source("resizewin.R")
resize.win(9,6)

#correct Time variable
pher$T= pher$time*30
pher$T.factor= as.factor(pher$T)

#rename treatment levels
pher$treatment2 <- factor(pher$Treatment, levels=c("HLB", "DPR"), labels=c("Control", "Diproline"))


#make new variable as ID (replicate, bead)
pher$ID <- paste(pher$Replicate, pher$Bead, sep = "-")

#check initial plot
qplot(T.factor,Cells, color = Treatment, data = pher,  geom = "boxplot") + facet_wrap(~Treatment, scales="free") 

#standardization

pher$CellsS=NA
k=split(pher, pher$treatment2)
pherstd <- lapply(k, function (x) scale(x[,c("Cells")], center=T, scale=T))
pher$CellsS=unsplit(pherstd, pher$treatment2)

qplot(T.factor,CellsS, color = treatment2, data = pher,  geom = "boxplot") + facet_wrap(~treatment2, scales="free") 


#baselining to 0 at time point 0
NT<-data.table(pher, key=c("ID"))

t1=NT[,list(treatment=treatment2, T=T, T.factor=T.factor, ID=ID, Cells=Cells, CellsS=CellsS,
            CellsBase=(CellsS-CellsS[1])), by=c("ID")]

pherbase1 <- t1 #DATA IS NOW CALLED pherBASE
pherbase <- subset(pherbase, pherbase$T>60, )


qplot(T.factor,CellsBase, color = treatment, data = pherbase,  geom = "boxplot") + facet_wrap(~treatment) +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))



#summary
source("summarySE.R")

pherbase.sum <- summarySE(pherbase, measurevar="CellsBase", groupvars=c("T", "treatment"))
pherbase1.sum <- summarySE(pherbase1, measurevar="CellsBase", groupvars=c("T", "treatment"))

ggplot(data=pherbase.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=10, size=1) + facet_grid(~treatment)

#model

exp=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(pherbase$treatment)), T=pherbase$T, 
                                   ID=as.numeric(as.factor(pherbase$ID)))))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)

pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

library(lawstat)
levene.test(pherbase$CellsBase, group=pherbase$ID, location="mean") #significant
levene.test(pherbase$CellsBase, group=pherbase$T, location="mean") # significant
levene.test(pherbase$CellsBase, group=pherbase$treatment, location="mean") #not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(CellsBase~treatment, data=pherbase)
boxplot(CellsBase~ID, data=pherbase)
boxplot (CellsBase~T, data=pherbase)

#fit a gls
Form <- formula (CellsBase ~ treatment*T)
pherbase.gls<- gls(Form, data=pherbase)


#nlme model
pherbase1.lme <- lme (Form, random = ~1|ID, method="REML", data=pherbase)

pherbase2.lme <- lme (Form, random = ~1|ID, correlation=corAR1(), method="REML", data=pherbase)

#pherbase3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), correlation=corAR1 (), method="REML", data=pherbase) 

#pherbase4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|T), correlation=corAR1 (), method="REML", data=pherbase) 

pherbase5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), 
                      correlation=corAR1(), method="REML", data=pherbase) #best is A5 

#pherbase6.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|treatment), varIdent (form=~1|ID)), correlation=corAR1(), method="REML", data=pherbase) 

#pherbase7.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|treatment), varIdent (form=~1|T)), correlation=corAR1(), method="REML", data=pherbase) 

anova(pherbase.gls, pherbase1.lme, pherbase2.lme, pherbase5.lme) #best is pherbase2.lme


summary(pherbase2.lme)
anova(pherbase2.lme)


#residuals
pherbase.E2<-resid(pherbase2.lme,type="normalized")
pherbase.F2<-fitted(pherbase2.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=pherbase.F2,y=pherbase.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(pherbase.E2~treatment,data=pherbase, main="Treatment",ylab=MyYlab)
plot(x=pherbase$T,y=pherbase.E2,main="Time",ylab=MyYlab,xlab="Time (sec)")
par(op)

xyplot (pherbase.E2 ~ T| treatment, data=pherbase, ylab="Residuals", xlab="Time (sec)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#let's plot this!

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 
library (AICcmodavg)

source("resizewin.R")
resize.win(9, 6)


#pherbase fit

pherbase.fit <- as.data.frame(predictSE.lme(pherbase2.lme, pherbase, se.fit = TRUE, level = 0,
                                            print.matrix = FALSE))

pherbase.fit$upr <- pherbase.fit$fit + (1.96 * pherbase.fit$se)
pherbase.fit$lwr <- pherbase.fit$fit - (1.96 * pherbase.fit$se)

pherbase.fit.combdata <- cbind(pherbase, pherbase.fit)


ggplot(data=pherbase1.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=15, size=1) + 
  geom_smooth(data=pherbase.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), method="lm", stat="identity", alpha=0.2)+ 
  scale_colour_manual(values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_manual (values = c(Control="#f1a340", Diproline="#998ec3"), name="Treatment")+ 
  labs(list(x = "Time (s)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=element_blank(), 
        legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_continuous (breaks=c(200, 400, 600)) 

#bw
resize.win(8,6)
ggplot(data=pherbase1.sum, aes(x=T, y=CellsBase, shape=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=15, size=1) + 
  geom_smooth(data=pherbase.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, linetype=treatment), 
              color="black", method="lm", stat="identity", alpha=0.2)+ 
  scale_linetype_manual(values = c(Control="dashed", Diproline="solid"), name="Treatment")+
  scale_shape_discrete(name="Treatment")+
  labs(list(x = "Time (s)", y = "Normalized cell count"))+ 
  theme(axis.text=text, axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5), legend.title=element_blank(), 
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(), legend.key.width=unit(1.8,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 
