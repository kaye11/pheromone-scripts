HLBspeed <- read.csv("d:/Karen's/PhD/R program/Pheromone data/Processed_data/HLBspeed.csv", sep=";")
DPRspeed <- read.csv("d:/Karen's/PhD/R program/Pheromone data/Processed_data/DPRspeed.csv", sep=";")



#for max number of tracks with the start of each track forced to 0
t1=HLBspeed
NT = data.table(t1)
NT2=NT[, time := seq(from = 1L, by = 1L, length.out = .N), by = ARF]

con=as.data.frame(NT2)
con$A=as.numeric(con$ARF)
con.count2=ddply(con,~time,summarise,Con=(length(ARF)))

t2=DPRspeed
NT = data.table(t2)
NT2=NT[, time := seq(from = 1L, by = 1L, length.out = .N), by = ARF]

DPR=as.data.frame(NT2)
DPR$A=as.numeric(DPR$ARF)
DPR.count2=ddply(DPR,~time,summarise,DPR=(length(A)))

#data binding (DPR=A, Con=Con)
dpr=subset(DPR.count2, DPR.count2$time<601, )
hlb=subset(con.count2, con.count2$time<601, )

CB2=cbind (dpr, hlb)

ggplot(CB2, aes(time)) + geom_line(aes(y = DPR, colour = "DPR"),  DPRze=1) + 
  geom_line(aes(y = Con, colour = "Control"), DPRze=1)+ 
  labs(list(x = "Time (s)", y = "Number of Tracks")) 


ggplot(CB2, aes(time)) + geom_line(aes(y = DPR, colour = "DPR")) + 
  geom_line(aes(y = Con, colour = "Control"))+ 
  labs(list(x = "Time (s)", y = "Number of Tracks")) 
