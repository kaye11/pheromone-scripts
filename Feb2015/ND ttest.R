library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)

pher360 <- subset (pherspeed, pherspeed$time<361)

#GD
NT = data.table(pher360, key="ARF")

GDtable <- NT[, list(GD=max(GD), T=max(T)), by=c("ARF", "cond")] 

qplot(cond,GD, data = GDtable,  geom = "boxplot") + labs(list(x = "Experimental Condition", y = "Gross distance traveled (µm)"))

shapiro.test(GDtable$GD) #normal
bartlett.test(GD ~ cond, data=GDtable) #homogenous

t.test(GD~cond, data=GDtable) # significant difference in GD


grid.newpage()
text <- element_text(size = 30, face="bold") #change the size of the axes
theme_set(theme_bw())


ggplot(data=GDtable, aes(x=cond, y=GD, color=cond)) +  geom_boxplot() +
  labs(list(x = "Experimental condition", y = "Gross distance traveled (µm)")) + scale_colour_manual(values = c("blue", "red"))+
  theme(axis.text.y=element_text(size=30), axis.text.x=element_text(size=30, face="bold"),
        axis.title=element_text(size=35,face="bold"),  axis.title.x = element_blank(),
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ND

## Split the data
phersub <- subset (pherspeed, pherspeed$time<361)
dfs <- split(pherspeed,pherspeed$ARF)

first <- NT [, head (.SD, 1), by=c("ARF")]
last <- NT [, tail (.SD, 1), by=c("ARF")]

fl<- rbind (first, last)

first2 <- subset (first, select=c("ARF", "X", "Y", "time", "cond" ))
last2 <-  subset (last, select=c("ARF", "X", "Y", "time" , "cond"))

fl2 <- rbind(first2, last2)

NT5=data.table(fl2, key="ARF")
NT6=NT5[, ND:=sqrt(diff(X)^2 + diff(Y)^2), by=c("ARF") ]

NTlast <- NT6 [, tail (.SD, 1), by=c("ARF")]

ND360=NTlast

shapiro.test(ND360$ND) # not normal
library(lawstat)
levene.test(ND360$ND, group=ND360$cond, location="mean") #homogenous variances

wilcox.test(ND~cond, data=ND360, paired=FALSE) #significantly different

grid.newpage()
text <- element_text(size = 30, face="bold") #change the size of the axes
theme_set(theme_bw())


ggplot(data=ND360, aes(x=cond, y=ND, color=cond)) +  geom_boxplot() +
  labs(list(y = "Net distance traveled (µm)")) + scale_colour_manual(values = c("blue", "red"))+
  theme(axis.text.y=element_text(size=30), axis.text.x=element_text(size=30, face="bold"),
        axis.title=element_text(size=35,face="bold"), axis.title.x = element_blank(),
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

