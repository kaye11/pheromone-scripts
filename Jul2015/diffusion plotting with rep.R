library(gdata)
library(data.table)
library(ggplot2)
library(grid)
library(gtable)
library(plyr)

'rep1' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/rep1.csv", sep=";")
'rep2' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/rep2.csv", sep=";")
'rep3' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/rep3.csv", sep=";")

all <- rbind((rep1 [, c (8, 11, 12)]), (rep2 [, c (8, 11, 12)]), 
               (rep3 [, c (8, 11, 12)]))

dev.new(width=6, height=9)
source("resizewin.R")
resize.win(9,6)

grid.newpage()
text <- element_text(size = 18, face="bold") #change the size of the axes
theme_set(theme_bw())

ggplot (data=all, aes(x=rad2, y=nMsq))+geom_line(size=2)+ 
  labs(x="Distance from bead (µm)", 
       y="nM diproline")+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"), axis.text=text, 
        axis.title.y = element_text(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

difall<- summarySE(all, measurevar="nMsq", groupvars=c("rad2"), na.rm=TRUE)

ggplot(data=difall, aes(x=rad2, y=nMsq)) + geom_line(size=0.8)+ 
  geom_ribbon(aes(ymin=nMsq-sd, ymax=nMsq+sd), alpha=0.5, size=3) +  
  labs(x="Distance from bead (µm)", 
       y="nM diproline")+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"),legend.position=c(0.12,0.10), 
        legend.key.width=unit(2,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

