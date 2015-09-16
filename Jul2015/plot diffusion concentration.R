
'min0' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/0.csv", sep=";")
'min2'<- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/2.csv", sep=";")
'min4' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/4.csv", sep=";")
'min6' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/6.csv", sep=";")
'min8' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/8.csv", sep=";")
'min10' <- read.csv("D:/Karen's/PhD/R program/Pheromone data/Processed_data/diffusion/10.csv", sep=";")

diff <- rbind((min0 [, c (8, 11, 12)]), (min2 [, c (8, 11, 12)]), (min4 [, c (8, 11, 12)]), (min6 [, c (8, 11, 12)]), 
              (min8 [, c (8, 11, 12)]), (min10 [, c (8, 11, 12)]))

diff2 <- rbind((min4 [, c (8, 11, 12)]), (min6 [, c (8, 11, 12)]), 
               (min8 [, c (8, 11, 12)]), (min10 [, c (8, 11, 12)]))


library(plotly)

py <- plotly(username="karen.bondoc", key="7v7sc4kylg")  # open plotly connection

p <- ggplot(diff, aes(x=rad2,y=rad2))

p <- p + geom_tile(aes(fill=nMsq))

py$ggplotly(p)

diff$time2=as.factor(diff$time)
ggplot (data=diff, aes(x=rad2, y=nMsq, linetype=time2))+geom_point(size=3, aes(shape=time2))

diff2$time2=as.factor(diff2$time)
ggplot (data=diff2, aes(x=rad2, y=nMsq, linetype=time2))+geom_point(size=2, aes(shape=time2))

ggplot (data=diff2, aes(x=rad2, y=nMsq, color=time2, shape=time2))+geom_point(size=3, aes(shape=time2))

pmin10 <- ggplot(min10, aes(x=rad2,y=rad2))
pmin10 <- pmin10 + geom_tile(aes(fill=nMsq))
py$ggplotly(pmin10)
