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
library("cowplot")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Input count file."),
  make_option(c("-c", "--contamination"),
    type = "character",
    default = NULL,
    help = "Input count file for contaminants."),
  make_option(c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-t", "--metadatat"),
    type = "character",
    default = NULL,
    help = "Metadata table with disease states for bacteria.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw resulting plot.",
    metavar = "character"),
  make_option(c("-a", "--contout"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw resulting contamination plot.",
    metavar = "character"),
  make_option(c("-s", "--sdepth"),
    type = "integer",
    default = 10000,
    help = "Subsampling depth.",
    metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

################
# Run Analysis #
################
input <- read.delim(opt$input, header=TRUE, sep="\t")
hc <- read.delim(opt$contamination, header=TRUE, sep="\t")
metadata <- read.delim(opt$metadata, header=FALSE, sep="\t")[,c(2,26,27,30)]
metabac <- read.delim(opt$metadatat, header=FALSE, sep="\t")[,c(2,26,27,30)]

metadata$class <- "virus"
metabac$class <- "bacteria"
metabac$V2 <- gsub(".cont", "", metabac$V2)
head(metadata)
head(metabac)
metatot <- rbind(metadata, metabac)

# # DNA concentration
# dnaconcplot <- ggplot(metadata, aes(x = V2, y = V26, fill = V30)) +
# 	theme_classic() +
# 	theme(
# 		axis.text.y=element_blank(),
# 		axis.ticks.y=element_blank(),
# 		legend.title=element_blank(),
# 		legend.position="none"
# 	) +
# 	geom_bar(stat="identity") +
# 	coord_flip() +
# 	scale_fill_manual(values = c(wes_palette("Royal1")[c(1,2,4)], wes_palette("Darjeeling")[5]), name = "Disease") +
# 	ylab("VLP Genomic DNA Yield (ng/uL)") +
# 	xlab("Prepared Samples")

# # Sampling Depth
# inputmerge <- merge(input, metadata, by.x="V2", by.y="V2")

# depthplot <- ggplot(inputmerge, aes(x = V2, y = total, fill = V30)) +
#     theme_classic() +
#     theme(
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.title=element_blank()
#     ) +
#     geom_bar(stat="identity") +
#     coord_flip() +
#     scale_fill_manual(values = c(wes_palette("Royal1")[c(1,2,4)], wes_palette("Darjeeling")[5]), name = "Disease") +
#     ylab("VLP Sequence Count") +
#     xlab("") +
#     geom_hline(yintercept = opt$sdepth, linetype = "dashed")

# gridplot <- plot_grid(dnaconcplot, depthplot, labels = c("A", "B"))

# pdf(opt$out, height=4, width=10)
#     gridplot
# dev.off()

# Human DNA contmaination removal
# Merge all of the data points

colnames(input) <- c("cleancount", "V2")
colnames(hc) <- c("contcount", "V2")
hc$V2 <- gsub(".cont", "", hc$V2)
head(input)
head(hc)
inmerge <- merge(input, hc, by = "V2")
inmerge$totalcount <- inmerge$cleancount + inmerge$contcount
inmerge$percentcont <- 100 * inmerge$contcount / inmerge$totalcount
head(inmerge)

tmm <- merge(inmerge, metatot, by.x="V2", by.y="V2")
head(tmm)

contplot <- ggplot(tmm, aes(x = class, y = percentcont, fill = class)) +
    theme_classic() +
    geom_boxplot()

pdf(opt$contout, height=4, width=6)
    contplot
dev.off()
