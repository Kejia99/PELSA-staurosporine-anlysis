setwd("set the working directory")

library(dplyr)
library(tidyverse)

##### import Bayes results
pelsa_k562 <- read.csv("export_bayes_PELSA_stau_K562_P1st.csv",header = T)
pelsa_hela <- read.csv("export_bayes_PELSA_stau_HeLa_P1st.csv",header = T)


#### hela
pelsa_hela $ alteration <- ifelse(pelsa_hela$Log2FC <0, "stabilized","destabilized")
pelsa_hela <- arrange(pelsa_hela,desc(alteration))
pelsa_hela$kinase.count <-0
table(pelsa_hela$is.kinase=="yes")
pelsa_hela$kinase.count[pelsa_hela$is.kinase=="yes"]<- c(1:254)
pelsa_hela$rank <-c(1:nrow(pelsa_hela))


#### k562
pelsa_k562 <- read.csv("export_bayes_PELSA_stau_K562_P1st.csv",header = T)
pelsa_k562 $ alteration <- ifelse(pelsa_k562$Log2FC <0, "stabilized","destabilized")
pelsa_k562 <- arrange(pelsa_k562,desc(alteration))
pelsa_k562$kinase.count <-0
table(pelsa_k562$is.kinase=="yes")
pelsa_k562$kinase.count[pelsa_k562$is.kinase=="yes"]<- c(1:254)
pelsa_k562$rank <-c(1:nrow(pelsa_k562))

str(pelsa_hela)
pelsa_hela$rank <- as.numeric(pelsa_hela$rank)
for (i in pelsa_hela$rank) {
  pelsa_hela[i,]$kinase.count <- ifelse(pelsa_hela[i,]$kinase.count==0,pelsa_hela[i-1,]$kinase.count,pelsa_hela[i,]$kinase.count)
  i<-i+1                
}


str(pelsa_k562)
pelsa_k562$rank <- as.numeric(pelsa_k562$rank)
for (i in pelsa_k562$rank) {
  pelsa_k562[i,]$kinase.count <- ifelse(pelsa_k562[i,]$kinase.count==0,pelsa_k562[i-1,]$kinase.count,pelsa_k562[i,]$kinase.count)
  i<-i+1                
}


write.csv(file=paste0("export_target_ranking_K562_stau",".csv"),pelsa_k562,row.names = F)
write.csv(file=paste0("export_target_ranking_HeLa_stau",".csv"),pelsa_hela,row.names = F)