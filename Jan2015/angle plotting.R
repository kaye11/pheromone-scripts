library(gridExtra)
library(plyr)
library(ggplot2)
library(ggthemes)


source("summarySE.R")

phersum <- summarySE(pherspeed, measurevar="angs", groupvars=c("cond", "T", "ARF"), na.rm=PASS)
phersum$cond2 <- factor(phersum$cond, levels=c("HLB", "DPR"), labels=c("Control", "DPR"))

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

ggplot(data=dfc2, aes(x=T, y=angs, shape=cond)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=angs-se, ymax=angs+se), width=10, size=1) + facet_grid(cond~., labeller=mf_labeller) + geom_hline(yintercept=0)+
  labs(list(x = "Time (s)", y = "Sine angle"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#normality test
by(dfc2$angs, dfc2$cond, shapiro.test) 
shapiro.test(dfc2$angs)

#If p values from normality tests are all greater 0.05 then it is acceptable. 

##histogram checking
ggplot(data = dfc2[dfc2$cond == "HLB",] , aes(x = angs))+ geom_histogram(binwidth = diff(range(dfc2$angs))/10)+ 
  facet_wrap(~cond, scale="free")

##histogram checking
ggplot(data = dfc2[dfc2$cond == "DPR",] , aes(x = angs))+ geom_histogram(binwidth = diff(range(dfc2$angs))/10)+ 
  facet_wrap(~cond, scale="free")

hist(dfc2$angs)


#homogeneity of variance
#  the null hypothesis is that all populations variances are equal; 
# the alternative hypothesis is that at least two of them differ.
qplot(cond, angs, data =dfc2,  geom = "boxplot")
bartlett.test(angs ~ cond, data=dfc2) #for normal data set
library(lawstat)
levene.test(dfc2$angs, group=dfc2$cond, location="mean")

#normal and homogenous, use t-test

t.test(angs~cond, data=dfc2) # no significant difference


#nlme

