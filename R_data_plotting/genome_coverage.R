# Rscript for plotting read depth over genome position in easy to present format.
# Currently tailored for HIV-1, edit xlim as appropriate for organism.
# Nicholas Gleadall - nick.gleadall@googlemail.com - 20/01/2015
# Usage == $ Rscript HIV-1_genome_coverage.R '.depth_file' 

args <- commandArgs(trailingOnly = TRUE)
#print(args)

out_file = args[1]
out_file = sub("(.*).bam.depth", "\\1.png", out_file)
print(out_file)
png(out_file, width=1048,height=440)


depth <- read.delim( args[1], header=F)
#print(depth)

depth_high <- depth[ depth$V3 > 0,]
#print(depth_high)

sample_name <- args[1]
sample_name <- sub(".*/(.*?).bam.*", "\\1", sample_name)
#print(sample_name)

plot(depth_high$V2, depth_high$V3, xlim=c(0,10000), ylim=c(1,10000), log="y", col='darkgreen', xlab='Genome position (bp) ', ylab='Read depth', type='l')

abline(h=1000, col=2, lty=3)	

title(main = paste("Depth plot for sample ", sample_name))



dev.off()	


