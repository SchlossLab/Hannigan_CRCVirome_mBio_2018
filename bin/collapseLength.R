# collapseLength.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

write("Summing contig abundance by cluster", stderr())

library("optparse")
library("plyr")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Contig length table to be collapsed.",
    metavar = "character"),
  make_option(c("-c", "--clusters"),
    type = "character",
    default = NULL,
    help = "Contig clustering table",
    metavar = "character"),
  make_option(c("-b", "--bacclusters"),
    type = "character",
    default = NULL,
    help = "Bacteria contig clustering table",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output table name",
    metavar = "character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

input <- read.delim(opt$input, head=FALSE, sep="\t")

clusters <- read.delim(opt$clusters, head=FALSE, sep=",")

bacclusters <- read.delim(opt$bacclusters, head=FALSE, sep=",")

input$V2 <- as.numeric(as.character(input$V2))

mcl <- rbind(clusters, bacclusters)

mmcl <- merge(input, mcl, by = "V1")

agtable <- ddply(input, c("V3"), summarize, scount=sum(V2))

write.table(
  x = agtable,
  file = opt$out,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write("Completed summing contig abundance by cluster", stderr())
