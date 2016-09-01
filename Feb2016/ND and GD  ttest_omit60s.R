library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)

pherspeed <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/pherspeed.csv", sep=";")

pher410 <- subset (pherspeed, pherspeed$T>60 & pherspeed$T<410, )

#GD
NT = data.table(pher410, key="ID")

GDtable <- NT[, list(GD=max(GD), T=max(T)), by=c("ID", "cond")] 

qplot(cond,GD, data = GDtable,  geom = "boxplot") + labs(list(x = "Experimental Condition", y = "Gross distance traveled (µm)"))

shapiro.test(GDtable$GD) #normal
bartlett.test(GD ~ cond, data=GDtable) #homogenous

t.test(GD~cond, data=GDtable) # significant difference in GD


grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 



ggplot(data=GDtable, aes(x=cond, y=GD, color=cond)) +  geom_boxplot() +
  labs(list(x = "Experimental condition", y = "Gross distance traveled (µm)")) + scale_colour_manual(values = c("blue", "red"))+
  theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20, face="bold"),
        axis.title=element_text(size=20,face="bold"),  axis.title.x = element_blank(),
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ND

## Split the data
phersub <- subset (pherspeed, pherspeed$T>60 & pherspeed$T<410,)
dfs <- split(pherspeed,pherspeed$ID)

first <- NT [, head (.SD, 1), by=c("ID")]
last <- NT [, tail (.SD, 1), by=c("ID")]

fl<- rbind (first, last)

first2 <- subset (first, select=c("ID", "X", "Y", "time", "cond" ))
last2 <-  subset (last, select=c("ID", "X", "Y", "time" , "cond"))

fl2 <- rbind(first2, last2)

NT5=data.table(fl2, key="ID")
NT6=NT5[, ND:=sqrt(diff(X)^2 + diff(Y)^2), by=c("ID") ]

NTlast <- NT6 [, tail (.SD, 1), by=c("ID")]

ND410=NTlast

shapiro.test(ND410$ND) # not normal
library(lawstat)
levene.test(ND410$ND, group=ND410$cond, location="mean") #not homogenous variances

wilcox.test(ND~cond, data=ND410, paired=FALSE) #significantly different


ggplot(data=ND410, aes(x=cond, y=ND, color=cond)) +  geom_boxplot() +
  labs(list(y = "Net distance traveled (µm)")) + scale_colour_manual(values = c("blue", "red"))+
  theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20, face="bold"),
        axis.title=element_text(size=20,face="bold"), axis.title.x = element_blank(),
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggplot(data=ND410, aes(x=cond, y=ND)) +  geom_boxplot() +
  labs(list(y = "Net distance traveled (µm)")) +
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=element_blank(), legend.text=element_text(size=20), 
        legend.title=element_blank(),legend.key.width=unit(2,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.margin=unit (0.5, "lines"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

source("summarySE.R")
source("resizewin.R")
NDsum <- summarySE(ND410, measurevar="ND", groupvars=c("cond"))


resize.win (6,4)
ggplot(data=NDsum, aes(x=cond, y=ND, fill=cond)) + 
  geom_bar(position = "dodge", stat="identity", color="black") +
  geom_errorbar(width=0.2, size=0.8, aes(ymin=ND-se, ymax=ND+se))+ 
  scale_fill_manual(values = c(Control="white", Diproline="gray"), name="Treatment")+
  labs(y = "Net distance traveled (µm)") + 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(), 
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        legend.title=element_blank(),legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

#colored
ggplot(data=NDsum, aes(x=cond, y=ND, fill=cond)) + 
  geom_bar(position = "dodge", stat="identity", color="black") +
  geom_errorbar(width=0.2, size=0.8, aes(ymin=ND-se, ymax=ND+se))+ 
  scale_fill_manual(values = c(Control="#f1a340", Diproline="#998ec3"))+
  labs(y = "Net distance traveled (µm)") + 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(), 
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        legend.title=element_blank(),legend.key.width=unit(1,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))


