# PlotSequencingDepth.R
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
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Input count file.",
    metavar = "character"),
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
input <- read.delim("./data/ProjectSeqDepth.tsv", header=FALSE, sep="\t")
metadata <- read.delim("./data/metadata/NexteraXT003Map.tsv", header=FALSE, sep="\t")[,c(2,26,27,30)]
head(input)
head(metadata)

inputmerge <- merge(input, metadata, by.x="V2", by.y="V2")

dnaconcplot <- ggplot(inputmerge, aes(x = V2, y = V1, fill = V30)) +
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
    xlab("Prepared Samples") +
    geom_hline(yintercept = 10000, linetype="dashed")
dnaconcplot
