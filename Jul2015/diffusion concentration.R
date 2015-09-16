library(gdata)
library(data.table)
library(ggplot2)
library(grid)
library(gtable)

bead2=15*0.0001 #in cm
rad=seq(1, 5000, 1)*0.0001
D= 10^-5
Ci=3.610971141E-09*0.01 #in mol/cm3
i=Ci*4*pi*bead2*D
i #in mol/s

var2=as.data.frame(rad)

var2$mcm=i/(4*pi*D*var2$rad)
var2$mcmsq=i/(4*pi*D*sqrt(var2$rad))
var2$M=var2$mcm*10^3
var2$Msq=var2$mcmsq*10^3
var2$uM=var2$M*10^6
var2$uMsq=var2$Msq*10^6
var2$rad2=var2$rad/0.0001
var2$ts=var2$rad^2/D
var2$ts2=var2$rad2^2/1000
var2$nMsq=var2$Msq*10^9

var=subset(var2, var2$rad2>20, )

grid.newpage()
text <- element_text(size = 20, face="bold") #change the size of the axes
theme_set(theme_bw())

ggplot (data=var2, aes(x=rad2, y=nMsq))+geom_line(size=2)+
  labs(x="Distance from bead (µm)", 
       y="nM diproline")+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"), axis.text=text, 
        axis.title.y = element_text(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

write.table (var2, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/rep3.csv", 
             sep=";", col.names=T, row.names=F)

p <- ggplot(pp(20), aes(x=x,y=y))