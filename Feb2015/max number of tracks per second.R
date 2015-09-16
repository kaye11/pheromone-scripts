##for checking max. no of tracks per sec

## for both HLBtrol and treatment
# DPR
DPRspeed$A=as.numeric(DPRspeed$A)
DPR.count=ddply(DPRspeed,~time,summarise,DPR=(length(unique((ARF)))))


#HLB
HLBspeed$A=as.numeric(HLBspeed$A)
HLB.count=ddply(HLBspeed,~time,summarise,HLB=(length(unique((ARF)))))


#data binding (DPR=A, HLB=HLB)
CB=cbind(DPR.count, HLB.count)

ggplot(CB, aes(time)) + geom_line(aes(y = DPR, colour = "DPR"),  DPRze=1) + 
  geom_line(aes(y = HLB, colour = "HLBtrol"), DPRze=1)+ 
  labs(list(x = "Time (s)", y = "Number of Tracks")) 

