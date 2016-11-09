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


