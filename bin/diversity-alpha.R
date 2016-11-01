# diversity-alpha.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

###################
# Set Environment #
###################

write("Calculating alpha diversity", stderr())

library("optparse")
library("ggplot2")
library("plotROC")
library("reshape2")
library("wesanderson")
library("vegan")
library("plyr")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "OTU abundance table straight from Mothur.",
    metavar = "character"),
  make_option(c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-s", "--subsample"),
    type = "integer",
    default = 10000,
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw results.",
    metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

################
# Run Analysis #
################
input <- read.delim(opt$input, header=TRUE, sep="\t")
head(input)
datadisease <- read.delim(opt$metadata, header=FALSE, sep="\t")[,c(2,30)]
head(datadisease)

# Format for vegan
inputcast <- dcast(input, V2 ~ V1, value.var = "sum")
inputcast[is.na(inputcast)] <- 0
row.names(inputcast) <- inputcast$V2
inputcast <- inputcast[,-1]

castsum <- rowSums(inputcast)

inputcast <- rrarefy(inputcast, sample=opt$subsample)
inputcast <- inputcast[c(rowSums(inputcast) %in% opt$subsample),]

inputshannon <- data.frame(rownames(inputcast) ,diversity(inputcast, index="shannon"))
colnames(inputshannon) <- c("SampleID", "ShannonDiv")
mergedshannon <- merge(inputshannon, datadisease, by.x="SampleID", by.y="V2")
mergedshannon <- mergedshannon[!c(mergedshannon$V30 %in% "Negative"),]
shannonplot <- ggplot(mergedshannon, aes(x=V30, y=ShannonDiv, fill=V30)) +
    theme_classic() +
    theme(legend.position="none") +
    geom_boxplot(notch=TRUE) +
    xlab("Health Status") +
    ylab("Shannon Diversity")
pairwise.wilcox.test(x=mergedshannon$ShannonDiv, g=mergedshannon$V30, p.adjust.method="bonferroni")

inputrich <- data.frame(rownames(inputcast), rarefy(inputcast, sample=opt$subsample))
colnames(inputrich) <- c("SampleID", "Rich")
mergedrich <- merge(inputrich, datadisease, by.x="SampleID", by.y="V2")
mergedrich <- mergedrich[!c(mergedrich$V30 %in% "Negative"),]
richplot <- ggplot(mergedrich, aes(x=V30, y=Rich, fill=V30)) +
    theme_classic() +
    theme(legend.position="none") +
    geom_boxplot(notch=TRUE) +
    # geom_segment(x=1, xend=3, y=150, yend=150) + annotate("text", x=2, y=152, label="*", size=10) +
    xlab("Health Status") +
    ylab("Richness")
pairwise.wilcox.test(x=mergedrich$Rich, g=mergedrich$V30, p.adjust.method="bonferroni")

pdf(opt$out, height=4.5, width=6)
    shannonplot
    richplot
dev.off()

write("Done calculating alpha diversity", stderr())
