# diversity-beta.R
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
library("cowplot")

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
  make_option(c("-d", "--divmetric"),
    type = "integer",
    default = "bray",
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
datadisease <- read.delim(opt$metadata, header=FALSE, sep="\t")[,c(2,30)]

# Format for vegan
inputcast <- dcast(input, V2 ~ V1, value.var = "sum")
inputcast[is.na(inputcast)] <- 0
row.names(inputcast) <- inputcast$V2
inputcast <- inputcast[,-1]

castsum <- rowSums(inputcast)

inputcast <- rrarefy(inputcast, sample=opt$subsample)
inputcast <- inputcast[c(rowSums(inputcast) %in% opt$subsample),]

# Make the distance matrix
castdist <- as.data.frame(as.matrix(vegdist(inputcast, method=opt$divmetric)))
ORD_NMDS <- metaMDS(castdist,k=2)
ORD_FIT = data.frame(MDS1 = ORD_NMDS$points[,1], MDS2 = ORD_NMDS$points[,2])
ORD_FIT$SampleID <- rownames(ORD_FIT)
NMDS_AND_MAP <- merge(ORD_FIT, datadisease, by.x="SampleID", by.y="V2")

# Plot everything
plotnmds <- ggplot(NMDS_AND_MAP, aes(x=MDS1, y=MDS2, colour=V30)) +
    theme_classic() +
    geom_point() +
    scale_fill_manual(values = wes_palette("Royal1"))

# Trying dispersion
catdistnames <- castdist
catdistnames$names <- row.names(catdistnames)
mergedist <- merge(catdistnames, datadisease, by.x="names", by.y="V2")
dist <- vegdist(inputcast, method=opt$divmetric)
mod <- betadisper(dist, mergedist[,length(mergedist)])
anova(mod)
permutest(mod, pairwise = TRUE)
mod.HSD <- TukeyHSD(mod)

moddf <- as.data.frame(mod.HSD$group)
moddf$comparison <- row.names(moddf)
limits <- aes(ymax = upr, ymin=lwr)
plotdiffs <- ggplot(moddf, aes(y=diff, x=comparison)) +
    theme_classic() +
    geom_pointrange(limits) +
    geom_hline(yintercept=0, linetype = "dashed") +
    coord_flip() +
    ylab("Differences in Mean Levels of Group") +
    xlab("")

gridplot <- plot_grid(plotnmds, plotdiffs, labels = c("A", "B"))

pdf(file=opt$out, width=10, height=4)
    gridplot
dev.off()
