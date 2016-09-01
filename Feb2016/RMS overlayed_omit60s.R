HLB.RMS <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/RMS data/ control_omit60s_410s .csv", sep=";")
DPR.RMS <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/RMS data/ dpr_omit60s_410s .csv", sep=";")

library(ggplot2)
library(grid)
library(ggthemes)
library(gridExtra)
library(reshape2)

dev.new(width=6, height=9)
source("resizewin.R")
resize.win(8,6)

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

##RMS plotting overlayed
df <- data.frame(HLB.RMS, DPR.RMS)
df$DPR=df$MF.1
df$HLB=df$MF


mainplot=ggplot(df, aes(time, y = value)) + 
  geom_line(aes(y = HLB, linetype = "Control"), size=2) + 
  geom_line(aes(y = DPR, linetype = "Diproline"), size=2)+ labs(list(x = "Time (s)", y = "RMS (µm)")) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"),legend.position=c(0.85,0.10), 
        legend.key.width=unit(2.5,"cm"),legend.key.height=unit(0.8,"cm"), 
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), panel.border=element_blank(), 
        axis.line = element_line(colour = "black"), plot.margin = unit(c(3,0.5,0.5,0.5), "lines"),
        panel.grid.minor = element_blank()) + scale_linetype_manual(values=c("dashed", "solid")) 

pherspeed <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/pherspeed.csv", sep=";")
phersub <- subset (pherspeed, pherspeed$time<361)
ND360$cond2 <- factor(ND360$cond, levels=c("DPR", "HLB"), labels=c("Diproline", "Control"))

subplot=ggplot(data=ND360, aes(x=cond2, y=ND)) +  geom_boxplot() +
  labs(list(y = "Net distance \n traveled (µm)")) +
  theme(axis.text.y=element_text(size=15, face="bold"), axis.text.x=element_text(size=12, face="bold"),
        axis.title=element_text(size=15,face="bold"), axis.title.x = element_blank(),
        plot.title = element_text(size =5, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

vp <- viewport(.3, .79, .38, .45)

full <- function() {
  print(mainplot)
  theme_set(theme_bw(base_size = 12))
  theme_bw()
  print(subplot, vp = vp)
  theme_set(theme_bw())
}

full()


#bw

ggplot(df, aes(time, y = value)) + 
  geom_line(aes(y = HLB, linetype = "Control"), size=2) + 
  geom_line(aes(y = DPR, linetype = "Diproline"), size=2)+ labs(list(x = "Time (s)", y = "RMS (µm)")) + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position=c(0.80,0.10),
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=element_blank(), legend.text=element_text(size=20), 
        legend.title=element_blank(),legend.key.width=unit(2,"cm"),legend.key.height=unit(0.8,"cm"),  
        panel.margin=unit (0.5, "lines"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

#colored
ggplot(df, aes(time, y = value)) + 
  geom_line(aes(y = HLB, color = "Control"), size=2) + 
  geom_line(aes(y = DPR, color = "Diproline"), size=2)+ labs(list(x = "Time (s)", y = "RMS (µm)")) + 
  scale_color_manual(values =  c(Control="#f1a340", Diproline="#998ec3"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position=c(0.80,0.10),
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=element_blank(), legend.text=element_text(size=20), 
        legend.title=element_blank(), 
        panel.margin=unit (0.5, "lines"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

#colored without legends
ggplot(df, aes(time, y = value)) + 
  geom_line(aes(y = HLB, color = "Control"), size=2) + 
  geom_line(aes(y = DPR, color = "Diproline"), size=2)+ labs(list(x = "Time (s)", y = "RMS (µm)")) + 
  scale_color_manual(values =  c(Control="#f1a340", Diproline="#998ec3"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=element_blank(), legend.text=element_text(size=20), 
        legend.title=element_blank(), 
        panel.margin=unit (0.5, "lines"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))


