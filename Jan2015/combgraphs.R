library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)

source("summarySE.R")

countsum <- summarySE(phercount, measurevar="CountN", groupvars=c("Treatment","Time"))

grid.newpage()
text <- element_text(size = 35, face="bold") #change the size of the axes
theme_set(theme_bw()) 

countplot=ggplot(data = phercount, aes(x=Time,y=CountN, color=Treatment))+
  stat_smooth(method="loess", formula=y~x, size=1.5, span=0.7) + scale_colour_manual(values = c("#CC1100", "#00008B")) +
  labs(list(x = "Time (s)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold", vjust=0.20), 
        plot.title = element_text(size =35, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank())

#pherspeed$cond <- factor(pherspeed$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))

velplot=ggplot(data = pherspeed, aes(x=T,y=V, color=cond))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE) + scale_colour_manual(values = c("#CC1100", "#00008B")) +
  labs(list(x = "Time (s)", y = "Speed (µm/s)"))+ theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold", vjust=0.20), 
        plot.title = element_text(size =35, face="bold"),legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.text=element_text(size=35, face="bold"), panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank())


grid.arrange (countplot, velplot, ncol=1)

angsum <- summarySE(pherspeed, measurevar="angs", groupvars=c("cond", "T"), na.rm=TRUE)

#angsum$cond <- factor(angsum$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))

angplot=ggplot(data=angsum, aes(x=T, y=angs, shape=cond, color=cond)) + geom_point(size=8)+ 
  geom_errorbar(aes(ymin=angs-se, ymax=angs+se), width=10, size=1) + geom_hline(yintercept=0)+ 
  scale_colour_manual(values = c("#CC1100", "#00008B")) + 
  labs(list(x = "Time (s)", y = "Sine angle"))+ theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold", vjust=0.20), 
        plot.title = element_text(size =35, face="bold"),legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.text=element_text(size=35, face="bold"), panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank())

grid.arrange (countplot, velplot, angplot, ncol=1)

grid.arrange (countplot, velplot, ncol=1)


#distance

distplot=ggplot(data=pherdatacomp, aes(x=time, y=dist, color=cond)) +  
  stat_smooth(aes(group=cond), method="loess", size=2, se=TRUE)+
  labs(list(x = "Time (s)", y = "Sum distance \n from the bead (mm)")) +
  scale_colour_manual(values = c("lightcoral", "steelblue2"), labels=c("Control", "dSi")) + 
  facet_grid(.~bin, label=mf_labeller2, scales="free_y")+
  labs (color="Experimental condition")+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=0.50), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.01),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(0.5,1,1,0), "cm")) + 
  scale_x_discrete (breaks=c("0", "200", "400", "600"))

