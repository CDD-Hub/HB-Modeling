#!/usr/bin/env Rscript

args <- commandArgs(TRUE) 

names <- as.numeric(unlist(strsplit(args[1],"_")))
names <- names/100
data <- as.numeric(unlist(strsplit(args[2],"_")))

png(filename=args[3], width=600, height=600)
barplot(data, xlab = "probability value", ylab = "probability occurrence in matrix", col="grey30", names.arg=names)
dev.off()
