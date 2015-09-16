##source this! name the file into s1
library (data.table)
library (ggplot2)
library(grid)

##direct input from trackmate 
NT=data.table(s1, key="TRACK_ID")
NT=NT[, list(A=TRACK_ID, X=POSITION_X, Y=POSITION_Y, T=FRAME), by=c("TRACK_ID")]
s1a= NT[order(A, T)]
s1=s1a[,list(T=T, X=X, Y=Y, V=c(0, sqrt(diff(X)^2+diff(Y)^2))), by=c("A")]

##plotting

t1=s1

qplot(X, Y, data = t1, color = factor(A), group = factor(A))+
  guides(col = guide_legend(nrow = 25)) + scale_y_reverse()+ geom_path(aes(group=factor(A))) 

qplot(X, Y, data=t1, color = factor(A), group = factor(A))+ geom_point(size=1)+
  geom_path(lineend="round", size=1, arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.20, "inches")), 
            mapping=aes(group=factor(A))) + scale_y_reverse()+theme_classic() + 
  labs(list(title="Control Alox Bead"))+
  theme(axis.text=element_text(size=20, face="bold"), axis.title=element_text(size=20,face="bold"), 
        legend.position="none", plot.title = element_text(size =20, face="bold"))+
  annotate("path",x = 337 + 49*cos(seq(0,2*pi,length.out=100)), y=322 + 49*sin(seq(0,2*pi,length.out=100)), 
           color="black", size=2)  

##parameter computations
beadX<-as.numeric(readline("X position of the Bead?"))
beadY<-as.numeric(readline("Y position of the Bead?"))

library(data.table)
library(zoo)

deg<-180/pi


NT<-data.table(t1, key=c("A"))
NT[, V := c(0, V[2:(.N-1)], 0), by = A]

t1=NT[,list(T=T, X=X, Y=Y, V=c(0, sqrt(diff(X)^2+diff(Y)^2)), GD=cumsum(V), ND=sqrt((X-X[1])^2 + (Y-Y[1])^2)), 
      by=c("A")]

t1[,NGDR:=ND/GD]
#t1[,ED:=sqrt((X-X[.N])^2 + (Y-Y[.N])^2), by=A]


t1[,a:=c(NA,(X[2:(.N-1)]),NA),by=A]
t1[,b:=c(NA,(Y[2:(.N-1)]),NA),by=A]
t1[,c:=c(NA,(X[1:(.N-2)]),NA),by=A]
t1[,d:=c(NA,(Y[1:(.N-2)]),NA),by=A]

t1[, scalar:=(a-c)*(a-beadX)+(b-d)*(b-beadY)]
t1[, angle:=acos(((a-c)*(a-beadX)+(b-d)*(b-beadY))/(sqrt((a-c)^2+(b-d)^2)*sqrt((a-beadX)^2+(b-beadY)^2)))*deg]
t1[, angle2:=c(NA, (na.locf(angle [1:(.N-1)])), NA), by=A]
t1[, angs:=sin(angle2)]

t1$a=NULL
t1$b=NULL
t1$c=NULL
t1$d=NULL
t1$scalar=NULL
t1$angle=NULL
t1$angle2=NULL

##saving data
VN<- readline("What data did you analyse?")
Vid<-paste ("d:/Karen's/PhD/R program/Pheromone data/Processed_data/",VN,".csv")
write.table(t1, Vid, sep=";", col.names=T, row.names=F)
