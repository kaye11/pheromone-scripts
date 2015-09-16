library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)

pher360 <- subset (pherspeed, pherspeed$time<361)

#GD
NT = data.table(pherspeed360, key="ARF")

GDtable <- NT[, list(GD=max(GD), T=max(T)), by=c("ARF", "cond")] 

qplot(cond,GD, data = GDtable,  geom = "boxplot") + labs(list(x = "Experimental Condition", y = "Gross distance traveled (µm)"))

shapiro.test(GDtable$GD) #normal
bartlett.test(GD ~ cond, data=GDtable) #homogenous

t.test(GD~cond, data=GDtable) # significant difference in GD


grid.newpage()
text <- element_text(size = 40) #change the size of the axes
theme_set(theme_bw())


ggplot(data=GDtable, aes(x=cond, y=GD, color=cond)) +  geom_boxplot() +
  labs(list(x = "Experimental condition", y = "Gross distance traveled (µm)")) + scale_colour_manual(values = c("#CC1100", "#00008B")) +
    scale_x_discrete(breaks=c("DPR", "HLB"),
                     labels=c("Diproline", "Control")) +
  theme(axis.text.y=element_text(size=30), axis.title.y=element_text(size=40,face="bold", vjust=-0.05), 
        plot.title = element_text(size =40, face="bold"), axis.title.x = element_blank(), 
        axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,2), "cm"), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#ND

## Split the data
dfs <- split(pher360,pher360$AR)

## Find hypotenuse between first and last rows for each A
NDtable2 <- as.data.frame(lapply( dfs , function(x){
  j <- nrow(x)
  str <- x[1,c("X","Y")]
  end <- x[j,c("X","Y")]
  dist <- sqrt( sum( (end - str)^2 ) )
  return( dist )
} ))

library(reshape)
NDtable <- melt(NDtable2)

write.table (NDtable, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/NDtable360.csv", 
             sep=";", col.names=T, row.names=F)

#NDtable saved in processed data folder as NDtable and then imported as ND

shapiro.test(ND$ND) # not normal
library(lawstat)
levene.test(ND$ND, group=ND$cond, location="mean") #homogenous variances

wilcox.test(ND~cond, data=ND, paired=FALSE) #significantly different

grid.newpage()
text <- element_text(size = 30, face="bold") #change the size of the axes
theme_set(theme_bw())


ggplot(data=ND, aes(x=cond, y=ND, color=cond)) +  geom_boxplot() +
  labs(list(x = "Experimental condition", y = "Net distance traveled (µm)")) + scale_colour_manual(values = c("blue", "red"))+
  theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=25, face="bold"),
        axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =20, face="bold"),legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit(1, "lines"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

