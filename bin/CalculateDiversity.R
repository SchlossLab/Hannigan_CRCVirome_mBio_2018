# CalculateDiversity.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Libraries #
#################

library("ggplot2")
library("RColorBrewer")
library("optparse")
library("vegan")

#################################
# Parse Input from Command Line #
#################################

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
      help = "Input OTU table", metavar = "character"),
  make_option(c("-m", "--map"), type = "character", default = NULL,
      help = "Mapping file with associated metadata", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
      help = "Output file for count summary", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###################
# Set Subroutines #
###################
AlphaDiversity <- function(intable, mapping) {
	# Fix the names in the headers
	names(intable) <- gsub("_R2", "", names(intable))
	# Get rid of the missing values
	tabledrop <- intable[, -grep("X", colnames(intable))]
	ttable <- as.matrix(t(tabledrop)[-1,])
	class(ttable) <- "numeric"
	divtable <- diversity(ttable, index=c("shannon"))
	divtable <- data.frame(divtable)
	colnames(divtable) <- "shannon"
	divtable$ID <- row.names(divtable)
	# Merge in the mapping file
	mapsub <- mapping[,c(2,30)]
	merged <- merge(divtable, mapsub, by.x="ID", by.y="V2")
	plot <- ggplot(merged, aes(x=V30, y=shannon, colour=V30)) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
		geom_jitter() +
		scale_colour_brewer(palette="Set1") +
		geom_boxplot(outlier.colour=NA, alpha=0) +
		ggtitle("Virome Shannon Diversity")
	return(plot)
}

BetaDiversity <- function(intable, mapping) {
	# Fix the names in the headers
	names(intable) <- gsub("_R2", "", names(intable))
	# Get rid of the missing values
	tabledrop <- intable[, -grep("X", colnames(intable))]
	ttable <- as.matrix(t(tabledrop)[-1,])
	class(ttable) <- "numeric"
	# Create dist matrix
	gothedistance <- vegdist(ttable, method = "bray")
	ordnmds <- metaMDS(gothedistance,k=2)
	ordnmdsfit = data.frame(MDS1 = ordnmds$points[,1], MDS2 = ordnmds$points[,2])
	ordnmdsfit$ID <- rownames(ordnmdsfit)
	# Merge in the mapping file
	mapsub <- mapping[,c(2,30)]
	merged <- merge(ordnmdsfit, mapsub, by.x="ID", by.y="V2")
	plot <- ggplot(merged, aes(x=MDS1, y=MDS2, colour=V30)) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
		geom_point(size=1.5) +
		scale_colour_brewer(palette="Set1") +
		ggtitle("Virome Bray-Curtis")
	return(plot)
}


############
# Run Data #
############

input <- read.delim(file=opt$input, sep="\t", header=T)
map <- read.delim(file=opt$map, sep="\t", header=F)

pdf(file="../figures/alphadiversity.pdf", width=4, height=4)
	AlphaDiversity(input, map)
dev.off()

pdf(file="../figures/betadiversity.pdf", width=5, height=4)
	BetaDiversity(input, map)
dev.off()
