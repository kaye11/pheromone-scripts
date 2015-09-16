library(ggplot2)
library(grid)
library(ggthemes)
library(gridExtra)
library(reshape2)

grid.newpage()
text <- element_text(size = 20, face="bold") #change the size of the axes
theme_set(theme_bw()) 

#21-DPR013
DPR01321 <- qplot(X, Y, data = pherdatacomp [pherdatacomp$AR== "21-DPR013", ]) +
  scale_colour_gradient(low="lightgrey", high="black") + 
  geom_path(aes(color=V), size=2, lineend="square", arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.1, "inches"))) + 
  annotate("path",x = 278 + 16*cos(seq(0,2*pi,length.out=100)), y=270 + 16*sin(seq(0,2*pi,length.out=100)), size=2)+
  xlim(200, 520) + ylim (200, 380)+
  labs(list(X="X", Y = "Y", color="Velocity (µm/s)")) + scale_y_reverse() +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold"), 
        plot.title = element_text(size =35, face="bold"),
        legend.key.width=unit(1,"cm"),legend.key.height=unit(1,"cm"), legend.position=c(0.75, 0.92), 
        legend.direction="horizontal", 
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=element_text(size=15), 
        panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  


#0-DPR018
DPR0180 <- qplot(X, Y, data = pherdatacomp [pherdatacomp$AR== "0-DPR018", ]) +
  scale_colour_gradient(low="lightgrey", high="black") + 
  geom_path(aes(color=V), size=2, lineend="square", arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.1, "inches"))) + 
  annotate("path",x = 277 + 17*cos(seq(0,2*pi,length.out=100)), y=281 + 17*sin(seq(0,2*pi,length.out=100)), size=2)+
  xlim(0, 400) + ylim (200, 380)+
  labs(list(X="X", Y = "Y", color="Velocity (µm/s)")) + scale_y_reverse() +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold"), 
        plot.title = element_text(size =35, face="bold"),
        legend.key.width=unit(1,"cm"),legend.key.height=unit(1,"cm"), legend.position=c(0.75, 0.92), 
        legend.direction="horizontal", 
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=element_text(size=15), 
        panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#0-DPR022
DPR0220 <- qplot(X, Y, data = pherdatacomp [pherdatacomp$AR== "0-DPR022", ]) +
  scale_colour_gradient(low="lightgrey", high="black") + 
  geom_path(aes(color=V), size=2, lineend="square", arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.1, "inches"))) + 
  annotate("path",x = 248 + 12.6*cos(seq(0,2*pi,length.out=100)), y=260 + 12.6*sin(seq(0,2*pi,length.out=100)), size=2)+
  xlim(200, 500) + ylim (180, 420)+
  labs(list(X="X", Y = "Y", color="Velocity (µm/s)")) + scale_y_reverse() +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold"), 
        plot.title = element_text(size =35, face="bold"),
        legend.key.width=unit(1,"cm"),legend.key.height=unit(1,"cm"), legend.position=c(0.75, 0.92), 
        legend.direction="horizontal", 
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=element_text(size=15), 
        panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#229-DPR018
DPR018229 <- qplot(X, Y, data = pherdatacomp [pherdatacomp$AR== "229-DPR018", ]) +
  scale_colour_gradient(low="lightgrey", high="black") + 
  geom_path(aes(color=V), size=2, lineend="square", arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.1, "inches"))) + 
  annotate("path",x = 277 + 17*cos(seq(0,2*pi,length.out=100)), y=281 + 17*sin(seq(0,2*pi,length.out=100)), size=2)+
  xlim(150, 420) + ylim (200, 600)+
  labs(list(X="X", Y = "Y", color="Velocity (µm/s)")) + scale_y_reverse() +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold"), 
        plot.title = element_text(size =35, face="bold"),
        legend.key.width=unit(1,"cm"),legend.key.height=unit(1,"cm"), legend.position=c(0.75, 0.92), 
        legend.direction="horizontal", 
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=element_text(size=15), 
        panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#175-HLB012
HLB012175 <- qplot(X, Y, data = pherdatacomp [pherdatacomp$ARF== "72-HLB014", ]) +
  scale_colour_gradient(low="lightgrey", high="black") + 
  geom_path(aes(color=V), size=2, lineend="square", arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.1, "inches"))) + 
  annotate("path",x = 257 + 14.35*cos(seq(0,2*pi,length.out=100)), y=258 + 14.35*sin(seq(0,2*pi,length.out=100)), size=2)+
  xlim(0, 600) + ylim (0, 600)+
  labs(list(X="X", Y = "Y", color="Velocity (µm/s)")) + scale_y_reverse() +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=35,face="bold"), 
        plot.title = element_text(size =35, face="bold"),
        legend.key.width=unit(1,"cm"),legend.key.height=unit(1,"cm"), legend.position=c(0.75, 0.92), 
        legend.direction="horizontal", 
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=element_text(size=15), 
        panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




qplot(X, Y, data = pherdatacomp [pherdatacomp$rep=="HLB014", ], color = factor(ARF), group = factor(ARF))+
  guides(col = guide_legend(nrow = 25)) + geom_path(aes(group=factor(ARF))) 
