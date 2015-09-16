library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
s=mgcv:::s

pher$cond <- factor(pher$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))


##check for collinearity and correlation, this only applies to the explanatory variables!
source("vif.R")
source ("AED.R")
exp=as.data.frame(data.table(cbind(cond=pher$cond, T=pher$time, A=pher$ARF)))
cor(exp, method = "spearman")

vif_func(in_frame=exp,thresh=5,trace=T)
corvif(exp)
pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#transform V to Vlog+1
pher$Vlog=log(pher$V+1)

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~cond, data=pher)
boxplot(Vlog~ARF, data=pher)
boxplot (Vlog~time, data=pher)
boxplot (Vlog~T, data=pher)


#summaries
source("summarySE.R")
Speedsum <- summarySE(pher, measurevar="Vlog", groupvars=c("ARF"))
SpeedsumT <- summarySE(pher, measurevar="Vlog", groupvars=c("T"))
SpeedsumCond <- summarySE(pher, measurevar="Vlog", groupvars=c("cond"))
SpeedsumCondV <- summarySE(pher, measurevar="V", groupvars=c("cond"))
SpeedsumCondTcond <- summarySE(pher, measurevar="V", groupvars=c("cond", "time"))

#levene test
library(lawstat)
levene.test(pher$Vlog, group=pher$cond, location="mean") #unequal
levene.test(pher$Vlog, group=pher$ARF, location="mean") #unequal
levene.test(pher$Vlog, group=pher$T, location="mean") #unequal
levene.test(pher$Vlog, group=pher$time, location="mean") #unequal


#models
#gam formula
PA <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cr"), method="REML", data = pher) #best
PA1 <- gamm (Vlog~s(time, by=cond, bs="fs", xt="cs"), method="REML", data = pher) 

anova(PA$lme, PA1$lme)

f1 <- Vlog~s(time, by=cond, bs="fs", xt="cr")

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

ggplot(data = pher, aes(x=T,y=Vlog, color=cond))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE)

dev.new(width=6, height=9)
source("resizewin.R")
resize.win(8,10)


op=par(mfrow=c(2,1), mar=c(5.1,6.1,4.1,2.1))
plot(PA6$gam, cex.lab=1.7, cex.axis=1.7, xlab ="Time (s)")



plot(PA6$lme, cex.lab=1.7, cex.axis=1.7)


summary(PA6$gam)
summary(PA6$lme)

##plotting
grid.newpage()
text <- element_text(size = 18, face="bold") #change the size of the axes
theme_set(theme_bw())

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="cond") { 
    value[value=="DPR"] <- "Diproline"
    value[value=="HLB"]   <- "Control"
  }
  return(value)
}


ggplot(data = pher, aes(x=T,y=V))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE, color="black") +
  facet_grid(cond~., labeller=mf_labeller)+ 
  labs(list(x = "Time (s)", y = "Speed (µm/s)"))+ labs (color="Experimental condition")+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#extract estimates of gam
summary_model <- summary(PA6$gam)
summary_model$p.table
summary_model$s.table

p_table <- data.frame(summary_model$p.table)
p_table <- within(p_table, {lci <- Estimate - qnorm(0.975) * Std..Error
                            uci <- Estimate + qnorm(0.975) * Std..Error})
p_table

summary(glht(PA6$gam, linfct=mcp(cond = "Tukey")))

#speedsum

grid.newpage()
text <- element_text(size = 40) #change the size of the axes
theme_set(theme_bw())

ggplot(data=SpeedsumCondV, aes(x=cond, y=V, width=0.75)) + 
  geom_bar(aes(fill = cond), position = "dodge", stat="identity") +
  geom_errorbar(width=.1, size=1, aes(ymin=V-se, ymax=V+se))+
  labs(y = "Mean Speed (µm/s)") + 
  scale_x_discrete(breaks=c("DPR", "HLB"),
                   labels=c("Diproline", "Control")) +
  scale_fill_manual("Cond",values = c("#CC1100", "#00008B"))+
  theme(axis.text.y=element_text(size=30), axis.title.y=element_text(size=40,face="bold", vjust=-0.05), 
        plot.title = element_text(size =40, face="bold"), axis.title.x = element_blank(), 
        axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,2), "cm"), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +  
  scale_y_continuous(expand=c(0,0)) 

shapiro.test(pher$V) # not normal
library(lawstat)
levene.test(pher$V, group=pher$cond, location="mean") #homogenous variances

wilcox.test(V~cond, data=pher, paired=FALSE) #significantly different
