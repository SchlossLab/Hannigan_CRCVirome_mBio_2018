# QualityControlStats.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

###################
# Set Environment #
###################

write("PROGRESS: Calculating QC Stats", stderr())

library("optparse")
library("ggplot2")
library("wesanderson")

option_list <- list(
  make_option(c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw resulting plot.",
    metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

################
# Run Analysis #
################
metadata <- read.delim("./data/metadata/NexteraXT003Map.tsv", header=FALSE, sep="\t")[,c(2,26,27,30)]
head(metadata)

dnaconcplot <- ggplot(metadata, aes(x = V2, y = V26, fill = V30)) +
	theme_classic() +
	theme(
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		legend.title=element_blank()
	) +
	geom_bar(stat="identity") +
	coord_flip() +
	scale_fill_manual(values = wes_palette("Royal1")[c(1,2,4,3)]) +
	ylab("VLP Genomic DNA Yield (ng/uL)") +
	xlab("Prepared Samples")

pdf("./figures/qualitycontrol.pdf", height=5, width=6)
    dnaconcplot
dev.off()
