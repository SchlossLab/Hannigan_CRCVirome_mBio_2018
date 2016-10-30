# ReshapeAlignedAbundance.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

write("PROGRESS: Reshaping contig abundance table.", stderr())

set.seed(1234)

library("optparse")
library("reshape2")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Contig abundance table.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Reshaped output table.",
    metavar = "character"),
  make_option(c("-p", "--percentsubset"),
    type = "numeric",
    default = 0.1,
    help = "Percent to subset samples before printing out table.",
    metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

input <- read.delim(opt$input, head = FALSE, sep = "\t")

collapsed <- dcast(input, V1 ~ V3, value.var="V2")
colnames(collapsed)[1] <- "contig"

collapsed[is.na(collapsed)] <- 0

perc <- round(ncol(collapsed) * opt$percentsubset)

write(paste("PROGRESS: Final sample number is ", perc, sep = ""), stderr())

collapsed <- collapsed[ , c(1, sample(2:ncol(collapsed), perc, replace = FALSE))]

write.table(
  x = collapsed,
  file = opt$out,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

write("PROGRESS: Completed reshaping contig abundance table.", stderr())
