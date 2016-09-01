#read all data files in separate data frames
myDir <- "d:/Karen's/PhD/R program/Pheromone data/Processed_data/trackdata/"
filenames <- list.files(myDir) 
filenames <- filenames[grep("[.]csv", filenames)] 

data_names <- gsub("[.]csv", "", filenames) 

for(i in 1:length(filenames)) assign(data_names[i], read.csv(sep=";", file.path(myDir, filenames[i]))) 

ls() #lists all objects in your workspace

#preprocessing of data sets, including merging, making new variable factors, binning by time

#merge data sets vertically
DPR013$rep <- "DPR013" 
DPR018$rep <- "DPR018" 
DPR022$rep <- "DPR022" 
DPR <- rbind(DPR013, DPR018, DPR022)
DPR$cond="Diproline"

Con012 <- (`HLB012 `)
Con014 <- (`HLB014 `)
Con020 <- (`HLB020 `)

Con012$rep <- "HLB012" 
Con014$rep <- "HLB014" 
Con020$rep <- "HLB020" 
Con <- rbind (Con012, Con014, Con020)
Con$cond = "Control"

pher <- rbind (DPR, Con)


#binning by time
pher$time=pher$T #replicate T if something awful happens
tm=seq(0, 600, by = 30)
pher$time <- cut(pher$T, tm, include.lowest=T, labels=paste(tail (tm, -1L)))

#Vlog, ID and correcting other parameters
pher$Vlog  <- log(pher$V+1)
pher$timef <- as.factor(pher$time)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]} #to make factor into numeric

pher$time <- as.numeric.factor(pher$time)
pher$ID <- as.factor (paste(pher$A, pher$rep, sep = "-")) #make ID as factor
pher$cond <- as.factor(pher$cond)

#save data table
write.table (pher, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/pherspeed.csv", 
             sep=";", col.names=T, row.names=F)

#save data table for DPR
DPRtable=subset(pher, pher$cond=="Diproline", )
write.table (DPRtable, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/DPRspeed.csv", 
             sep=";", col.names=T, row.names=F)

#save data table for HLB
HLBtable=subset(pher, pher$cond=="Control", )
write.table (HLBtable, "d:/Karen's/PhD/R program/Pheromone data/Processed_data/HLBspeed.csv", 
             sep=";", col.names=T, row.names=F)
