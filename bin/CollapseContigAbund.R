# CollapseCotingAbund.R
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
    help = "Contig abundance table to be collapsed.",
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

input$V3 <- as.numeric(as.character(input$V3))

agtable <- ddply(input, c("V1", "V2"), summarize, sum=(sum(V3) / length(V3)) )

write.table(
  x = agtable,
  file = opt$out,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write("Completed summing contig abundance by cluster", stderr())
