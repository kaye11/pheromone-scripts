#preprocessing of data sets, including merging, making new variable factors, binning by time

#merge data sets vertically
DPR013$rep <- "DPR013" 
DPR018$rep <- "DPR018" 
DPR022$rep <- "DPR022" 
DPR <- rbind(DPR013, DPR018, DPR022)
DPR$cond="DPR"

HLB012$rep <- "HLB012" 
HLB014$rep <- "HLB014" 
HLB020$rep <- "HLB020" 
HLB <- rbind (HLB012, HLB014, HLB020)
HLB$cond = "HLB"

pher <- rbind (DPR, HLB)

# Combine A and rep
pher$AR <- paste(pher$A, pher$rep, sep = "-")
pher$ARF <- as.factor (pher$AR) #make AR as factor

#binning by time
pher$time=pher$T #replicate T if something awful happens
tm=seq(0, 600, by = 30)
pher$T <- cut(pher$time, tm, include.lowest=T, labels=paste(tail (tm, -1L)))

#save data table
write.table (pher, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/pherspeed.csv", 
             sep=";", col.names=T, row.names=F)

#save data table for DPR
DPRtable=subset(pher, pher$cond=="DPR", )
write.table (DPRtable, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/DPRspeed.csv", 
             sep=";", col.names=T, row.names=F)

#save data table for HLB
HLBtable=subset(pher, pher$cond=="HLB", )
write.table (HLBtable, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/HLBspeed.csv", 
             sep=";", col.names=T, row.names=F)
