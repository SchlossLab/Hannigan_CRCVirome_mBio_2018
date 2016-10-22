# PlotAnnotations.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Libraries #
#################

library("ggplot2")
library("RColorBrewer")
library("optparse")
library("plyr")
library("reshape2")

#################################
# Parse Input from Command Line #
#################################

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
      help = "Input OTU table", metavar = "character"),
  make_option(c("-m", "--map"), type = "character", default = NULL,
      help = "Mapping file with associated metadata", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###################
# Set Subroutines #
###################

PlotRelAbund <- function(relabundtable, mapp) {
	names(relabundtable) <- gsub("_R2", "", names(relabundtable))
	# collapse the rel abund table
	sumtable <- ddply(relabundtable,"Contig_ID",numcolwise(sum))
	# sumtable <- cbind(sumtable[,1], sweep(sumtable[,-1], 2, colSums(sumtable[,-1]), FUN="/"))
	# colnames(sumtable) <- colnames(ddply(relabundtable,"Contig_ID",numcolwise(sum)))
	melttable <- melt(sumtable)
	merged <- merge(mapp, melttable, by.x="V2", by.y="variable")

	colourCount <- length(unique(melttable$Contig_ID))
	getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

	# Get order
	MeanCounts <- ddply(melttable, "Contig_ID", summarize, mean=mean(value))
	MeanOrder <- MeanCounts[order(MeanCounts$mean, decreasing=TRUE)[],]

	TopAbund <- melttable[c(melttable$Contig_ID == MeanOrder[1,1]),]
	TopAbundOrder <- TopAbund[c(order(TopAbund$value, decreasing=TRUE)),]

	merged$V2 <- factor(merged$V2, levels = TopAbundOrder$variable)
	merged$Contig_ID <- factor(merged$Contig_ID, levels = MeanOrder$Contig_ID)

	plot <- ggplot(merged, aes(x=V2, y=value, fill=Contig_ID)) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
		geom_bar(stat="identity") +
		scale_fill_manual(values = getPalette(colourCount)) +
		facet_grid(~V30, scale="free")
	return(plot)
}

############
# Run Data #
############

input <- read.delim(file=opt$input, sep="\t", header=T)
map <- read.delim(file=opt$map, sep="\t", header=F)

pdf(file="../figures/annotatedviruses.pdf", width=10, height=6)
	PlotRelAbund(input, map)
dev.off()

