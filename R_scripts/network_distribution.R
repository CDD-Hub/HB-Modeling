#!/usr/bin/env Rscript

args <- commandArgs(TRUE) 
data <- as.numeric(unlist(strsplit(args[1],"_")))
png(filename=args[2], width=600, height=600)
if (args[3]=="Hbond") {
	barplot(data, ylim=c(0,50), xlab = "networks", ylab = "number of residues in each network (normalized)", col="brown2")
} else {
	barplot(data, ylim=c(0,25), xlab = "networks", ylab = "number of residues in each network (normalized)", col="brown2")
}
dev.off()

#Figure x : Number of residues in each hydrogen bond network, in descending order, normalized by the sum of acceptors and donors
