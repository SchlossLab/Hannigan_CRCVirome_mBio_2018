# RunReadCountStats.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Libraries #
#################

library("ggplot2")
library("RColorBrewer")
library("optparse")

#################################
# Parse Input from Command Line #
#################################

option_list = list(
	make_option(c("-c", "--counts"), type="character", default=NULL, 
			help="table with read count information", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, 
			help="output file for count summary", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#################
# Read In Files #
#################

COUNTS <- read.delim(file=opt$counts, sep="\t", header=F)

############
# Run Data #
############

Summary <- summary(COUNTS)

ComparePlot <- ggplot(COUNTS, aes(x=V1, y=V3, colour=V2)) + theme_classic() + geom_bar(stat = "identity") + scale_colour_brewer(palette="Set2")

pdf(file=opt$output, width=8, height=6)
	ComparePlot
dev.off()
