#!/usr/bin/env Rscript

args <- commandArgs(TRUE) 
data <- as.numeric(unlist(strsplit(args[1],"_")))
png(filename=args[2], width=600, height=600)
if (args[3]=="molpdf") {
	plot(data, xlab = "cycles", ylab = "molpdf score")
} else {
	if (args[3]=="DOPE") {
		plot(data, xlab = "cycles", ylab = "DOPE score")
	} else {
		if (args[3]=="GA341") {
			plot(data, xlab = "cycles", ylab = "GA341 score")
		} else {
			if (args[3]=="RMSD_template") {
				data2 <- as.numeric(unlist(strsplit(args[4],"_")))
				plot(data, xlab = "cycles", ylab = "RMSD score (VS template)")
				points(data2, col="turquoise3")
			} else {
				if (args[3]=="RMSD_query_structure") {
					data2 <- as.numeric(unlist(strsplit(args[4],"_")))
					plot(data, xlab = "cycles", ylab = "RMSD score (VS crystal structure)")
					points(data2, col="turquoise3")
				} else {
					if (args[3]=="sum_Hbonds") {
						plot(data, ylim = c(0,500), xlab = "cycles", ylab = "number of hydrogen bonds")
					} else {
						if (args[3]=="sum_probas") {
							plot(data, xlab = "cycles", ylab = "sum of probabilities in the matrix")
						} else {
							if (args[3]=="sum_restraints") {
								plot(data, xlab = "cycles", ylab = "sum of drawn restraints")
							}
						}
					}
				}
			}
		}
	}
}
dev.off()

#Figure x : lambda score evolution during the optimisation cycles
#Figure x : criteria value evolution during the optimisation cycles
