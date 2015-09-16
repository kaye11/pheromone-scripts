library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)

source("summarySE.R")



dev.new(width=6, height=9)
source("resizewin.R")
resize.win(9,6)


grid.newpage()
text <- element_text(size = 18, face="bold") #change the size of the axes
theme_set(theme_bw()) 

#COUNT
countsum <- summarySE(phercount, measurevar="CountN", groupvars=c("Treatment","Time"))
phercount$cond <- factor(phercount$Treatment, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))


ggplot(data = phercount, aes(x=Time,y=CountN, linetype=cond))+ 
  stat_smooth(method="loess", formula=y~x, size=1.5, color="black") +
  labs (list(x = "Time (s)", y = "Normalized cell count"))+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"),legend.position=c(0.90,0.10), 
        legend.key.width=unit(1.2,"cm"),legend.key.height=unit(0.8,"cm"), 
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), panel.border=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())+scale_x_continuous (breaks=c(200, 400, 600)) 

#SPEED
#pherspeed$cond <- factor(pherspeed$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))


mainplot=ggplot(data = pherspeed, aes(x=time,y=V, linetype=cond))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE, color="black") +
  labs(list(x = "Time (s)", y = "Speed (µm/s)"))+ theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"),legend.position=c(0.15,0.10), 
        legend.key.width=unit(1.2,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), plot.margin = unit(c(3,6,0.5,0.5), "lines"), panel.border=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())+scale_x_continuous (breaks=c(200, 400, 600)) 

SpeedsumCondV <- summarySE(pherspeed, measurevar="V", groupvars=c("cond"))

subplot=ggplot(data=SpeedsumCondV, aes(x=cond, y=V, width=0.75)) + 
  geom_bar(fill="white", position = "dodge", stat="identity", color="black") +
  geom_errorbar(width=.1, size=1, aes(ymin=V-se, ymax=V+se))+ 
  labs(y = "Mean speed (µm/s)") + 
  theme(axis.text.y=element_text(size=12, face="bold"), axis.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=12,face="bold"), axis.title.x = element_blank(),
        plot.title = element_text(size =5, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

vp <- viewport(.85, .82, .30, .37)

full <- function() {
  print(mainplot)
  theme_set(theme_bw(base_size = 12))
  theme_bw()
  print(subplot, vp = vp)
  theme_set(theme_bw())
}

full()

#ANGLE
angsum <- summarySE(pherspeed, measurevar="angs", groupvars=c("cond", "T"), na.rm=TRUE)

#angsum$cond <- factor(angsum$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))

ggplot(data=angsum, aes(x=T, y=angs, shape=cond, linetype=cond)) + geom_point(size=8)+ 
  geom_errorbar(aes(ymin=angs-se, ymax=angs+se), width=10, size=1) + geom_hline(yintercept=0)+ 
  labs(list(x = "Time (s)", y = "Sine angle"))+ theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"),legend.position=c(0.15,0.10), 
        legend.key.width=unit(2,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), panel.border=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())+ scale_x_continuous (breaks=c(200, 400, 600)) 


#distance

distall <- ddply(pherdist, c("time", "cond"), summarise,
                 N    = length(dist),
                 mean = mean(dist, na.rm=TRUE),
                 sumdist= sum(dist, na.rm=TRUE), 
                 sd   = sd(dist, na.rm=TRUE),
                 se   = sd / sqrt(N), na.rm=TRUE,
                 angs = mean(angs, na.rm=TRUE))

distall$distmm <- distall$sumdist/1000

distall$cond <- factor(distall$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))


ggplot(data = distall, aes(x=time,y=distmm, linetype=cond))+ 
  stat_smooth(method="gam", formula=y~s(x, k=6), size=2, se=TRUE, color="black") +
  labs(list(x = "Time (s)",y = "Sum distance \n from the bead (mm)"))+ theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"),legend.position=c(0.15,0.10), 
        legend.key.width=unit(1.2,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), panel.border=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank()) +scale_x_continuous (breaks=c(200, 400, 600)) 

#trajectories

qplot(X, Y, data = pherspeed [pherspeed$rep=="DPR018", ], color = factor(ARF), group = factor(ARF))+
  guides(col = guide_legend(nrow = 25)) + scale_y_reverse() + geom_path(aes(group=factor(ARF)))
