# diversity-beta.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

###################
# Set Environment #
###################

write("Calculating beta diversity", stderr())

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
    help = "OTU abundance table.",
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
    type = "character",
    default = "bray",
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw results.",
    metavar = "character"),
  make_option(c("-n", "--negout"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw results with neg controls (QC).",
    metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

################
# Run Analysis #
################

##### Beta Diversity #####

# Import otu abundance table
input <- read.delim(opt$input, header=TRUE, sep="\t")
# Get the metadata table
datadisease <- read.delim(opt$metadata, header=FALSE, sep="\t")[,c(2,30)]

# Format for vegan
inputcast <- dcast(input, V2 ~ V1, value.var = "sum")
inputcast[is.na(inputcast)] <- 0
row.names(inputcast) <- inputcast$V2
inputcast <- inputcast[,-1]

castsum <- rowSums(inputcast)

inputcast <- rrarefy(inputcast, sample=opt$subsample)
inputcast <- inputcast[c(rowSums(inputcast) %in% opt$subsample),]

# Make one without negative controls as well
negcastvector <- as.character(datadisease[c(datadisease$V30 %in% "Negative"),1])
negcast <- inputcast[!c(rownames(inputcast) %in% negcastvector),]

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

# Make the distance matrix without negative control data
castdistnoneg <- as.data.frame(as.matrix(vegdist(negcast, method=opt$divmetric)))
ORD_NMDS <- metaMDS(castdistnoneg,k=2)
ORD_FIT = data.frame(MDS1 = ORD_NMDS$points[,1], MDS2 = ORD_NMDS$points[,2])
ORD_FIT$SampleID <- rownames(ORD_FIT)
NMDS_AND_MAP <- merge(ORD_FIT, datadisease, by.x="SampleID", by.y="V2")

# Plot everything
nonegplotnmds <- ggplot(NMDS_AND_MAP, aes(x=MDS1, y=MDS2, colour=V30)) +
    theme_classic() +
    geom_point() +
    scale_fill_manual(values = wes_palette("Royal1"))

write("Dispersion with negative controls", stderr())
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

write("Dispersion without negative controls", stderr())
# Dispersion without negatives as well
catdistnames <- castdistnoneg
catdistnames$names <- row.names(catdistnames)
mergedist <- merge(catdistnames, datadisease, by.x="names", by.y="V2")
dist <- vegdist(negcast, method=opt$divmetric)
mod <- betadisper(dist, mergedist[,length(mergedist)])
anova(mod)
permutest(mod, pairwise = TRUE)
mod.HSD <- TukeyHSD(mod)

moddf <- as.data.frame(mod.HSD$group)
moddf$comparison <- row.names(moddf)
limits <- aes(ymax = upr, ymin=lwr)
nonegplotdiffs <- ggplot(moddf, aes(y=diff, x=comparison)) +
    theme_classic() +
    geom_pointrange(limits) +
    geom_hline(yintercept=0, linetype = "dashed") +
    coord_flip() +
    ylab("Differences in Mean Levels of Group") +
    xlab("")

##### Alpha Diversity #####

inputshannon <- data.frame(rownames(negcast) ,diversity(negcast, index="shannon"))
colnames(inputshannon) <- c("SampleID", "ShannonDiv")
mergedshannon <- merge(inputshannon, datadisease, by.x="SampleID", by.y="V2")
mergedshannon <- mergedshannon[!c(mergedshannon$V30 %in% "Negative"),]
shannonplot <- ggplot(mergedshannon, aes(x=V30, y=ShannonDiv, fill=V30)) +
    theme_classic() +
    theme(legend.position="none") +
    geom_boxplot(notch=FALSE) +
    xlab("Health Status") +
    ylab("Shannon Diversity") +
    scale_fill_manual(values = alpha(wes_palette("Royal1"), 0.75))
pairwise.wilcox.test(x=mergedshannon$ShannonDiv, g=mergedshannon$V30, p.adjust.method="bonferroni")

inputrich <- data.frame(rownames(negcast), rarefy(negcast, sample=opt$subsample))
colnames(inputrich) <- c("SampleID", "Rich")
mergedrich <- merge(inputrich, datadisease, by.x="SampleID", by.y="V2")
mergedrich <- mergedrich[!c(mergedrich$V30 %in% "Negative"),]
richplot <- ggplot(mergedrich, aes(x=V30, y=Rich, fill=V30)) +
    theme_classic() +
    theme(legend.position="none") +
    geom_boxplot(notch=FALSE) +
    xlab("Health Status") +
    ylab("Richness") +
    scale_fill_manual(values = alpha(wes_palette("Royal1"), 0.75))
pairwise.wilcox.test(x=mergedrich$Rich, g=mergedrich$V30, p.adjust.method="bonferroni")


alphaplot <- plot_grid(shannonplot, richplot, labels = c("C", "D"))
gridplot <- plot_grid(plotnmds, plotdiffs, labels = c("A", "B"))
gridplotnoneg <- plot_grid(nonegplotnmds, nonegplotdiffs, alphaplot, labels = c("A", "B"), ncol = 3)

pdf(file=opt$out, width=10, height=4)
    gridplotnoneg
dev.off()

pdf(file=opt$negout, width=10, height=4)
    gridplot
dev.off()
