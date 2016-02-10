
library(ggplot2)
library(reshape)

getwd()
setwd("/Users/Nick/Desktop/Boxplot_analysis/")
getwd()



data <- read.csv("NIPS_10samples_KIMFIX.csv", head=TRUE)
data



datam <- melt(data)
datam

ggplot(datam, aes(x = variable, y = value)) + geom_boxplot(width=1, outlier.shape = 21) + 
  ggtitle("Comparison of mutation prevalence seen in Illumina MiSeq sequence data from 10 NIBSC positive control samples\n") + theme_bw() + xlab("\nMutation") +
  ylab("Percent observed in read data")


