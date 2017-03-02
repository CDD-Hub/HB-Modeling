#!/usr/bin/env Rscript

# provide access to a copy of the command line arguments supplied
args <- commandArgs(TRUE)
# split, unlist and convert as numeric the first argument (vector of type "1_3_1_2_2_")
vect <- as.numeric(unlist(strsplit(args[1],"_")))
# create an empty numeric vector with the vector maximum (size of 3 in this example)
data <- vector("numeric", max(vect, na.rm = TRUE)+1)
# fill vector with the number of occurences between 0 and the maximum
for (i in 0:max(vect, na.rm = TRUE)+1) {
	data[i]=sum(vect==i-1, na.rm = TRUE)
}

# names vector contains labels for axis x
names <- seq(0,length(data)-1,1)
# create a barplot recorded in a PNG graph with specific size
png(filename=args[2], width=600, height=600)
barplot(data, ylim=c(0,500), xlab="number of hydrogen bonds", ylab="number of donor", names.arg=names, col="gray10")
dev.off()


