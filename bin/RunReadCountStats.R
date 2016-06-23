# RunReadCountStats.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Libraries #
#################

library("ggplot2")
library("RColorBrewer")
library("optparse")
library("gridExtra")
library("plyr")

#################################
# Parse Input from Command Line #
#################################

option_list <- list(
  make_option(c("-c", "--counts"), type = "character", default = NULL,
      help = "Table with read count information", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
      help = "Output file for count summary", metavar = "character"),
  make_option(c("-p", "--png"), type = "character", default = NULL,
      help = "PNG output file for count summary", metavar = "character"),
  make_option(c("-t", "--title"), type = "character", default = NULL,
      help = "Title for the resulting plot", metavar = "character"),
  make_option(c("-l", "--log"), action = "store_true", default = FALSE,
      help = "Should we use a log scale? [default %default]"),
  make_option(c("-r", "--remove"), action = "store_true", default = FALSE,
      help = "Remove negative controls from view [default %default]"),
  make_option(c("-y", "--ylabel"), type = "character", default = NULL,
      help = "Label for y axis", metavar = "character"),
  make_option(c("-m", "--mean"), type = "store_true", default=FALSE,
      help = "Specify whether mean and std dev should be calculated", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#################
# Read In Files #
#################

COUNTS <- read.delim(file=opt$counts, sep="\t", header=F)

if (opt$remove) {
  COUNTS <- COUNTS[c(!grepl("MG100008", COUNTS$V1)),]
  COUNTS <- COUNTS[c(!grepl("MG100013", COUNTS$V1)),]
}

############
# Run Data #
############

Summary <- summary(COUNTS)

ComparePlot <- ggplot(COUNTS, aes(x=V1, y=V3, fill=V2)) +
  theme_classic() +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab(opt$title) +
  xlab("Samples") +
  coord_flip() +
  ggtitle(opt$title)

if (opt$log) {
   ComparePlot <- ComparePlot + scale_y_log10()
}

if (opt$mean) {
  whyyougottabesomean <- ddply(
      COUNTS,
      "V2",
      summarize,
      mean = mean(V3),
      sem <- sd(V3)/length(V3)
    )
  ComparePlot <- ggplot(whyyougottabesomean, aes(x="Average", y=mean)) +
    theme_classic() +
    geom_bar(position="dodge") +
    geom_errorbar(sem, position="dodge", width=0.25)
}

pdf(file=opt$output, width=8, height=6)
  ComparePlot
dev.off()

png(file=opt$png, width=8, height=6, units="in", res=300)
  ComparePlot
dev.off()
