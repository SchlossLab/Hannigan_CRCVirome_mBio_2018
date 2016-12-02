# ContigLengthDist.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

library("optparse")
library("ggplot2")
library("cowplot")
library("plyr")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Contig length table.",
    metavar = "character"),
  make_option(c("-t", "--toplength"),
    type = "integer",
    default = 1,
    help = "Amount of top lengths to report.",
    metavar = "character"),
  make_option(c("-c", "--clusters"),
    type = "character",
    default = NULL,
    help = "Contig cluster id table.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output name.",
    metavar = "character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# Virus

input <- read.delim(opt$input, head = FALSE, sep = "\t")
colnames(input) <- c("ContigID", "Length")

clusters <- read.delim(opt$clusters, head = FALSE, sep = ",")
colnames(clusters) <- c("ContigID", "Cluster")

mergeclust <- merge(input, clusters, by = "ContigID")

# Get top lengths
topcontigsbylength <- ddply(mergeclust, "Cluster", function(x) head(x[order(x$Length, decreasing = TRUE) , ], opt$toplength))

write.table(
  x = topcontigsbylength,
  file = opt$out,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
