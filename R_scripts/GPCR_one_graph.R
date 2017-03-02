#!/usr/bin/env Rscript

args <- commandArgs(TRUE) 

png(filename="Hbond_network_test.png", width=600, height=600)
data<-read.table(paste(args[1], "network_distribution_Hbond_test.txt", sep="/"), sep="_", fill=TRUE, row.names=NULL)
data<-data[,-1]
line1<-unlist(data[1,])
plot(line1, xlim=c(1, ncol(data)), ylim=c(0,50), xlab = "networks", ylab = "number of residues for each network", col="forestgreen")
for (i in 2:nrow(data)){
	line<-unlist(data[i,])
	points(line, col="forestgreen")
}
dev.off()

png(filename="Hbond_network_test_norm.png", width=600, height=600)
data<-read.table(paste(args[1], "network_distribution_Hbond_test_norm.txt", sep="/"), sep="_", fill=TRUE, row.names=NULL)
data<-data[,-1]
line1<-unlist(data[1,])
plot(line1, xlim=c(1, ncol(data)), ylim=c(0,50), xlab = "networks", ylab = "number of residues for each network", col="forestgreen")
for (i in 2:nrow(data)){
	line<-unlist(data[i,])
	points(line, col="forestgreen")
}
dev.off()


