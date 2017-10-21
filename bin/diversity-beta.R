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
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
# Get the metadata table
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

# Format for vegan
inputcast <- dcast(input, V2 ~ V1, value.var = "sum")
inputcast[is.na(inputcast)] <- 0
row.names(inputcast) <- inputcast$V2
inputcast <- inputcast[,-1]
inputcast <- round(inputcast)
negcastvector <- as.character(datadisease[c(datadisease$V30 %in% "Negative"),1])
negcast <- inputcast[!c(rownames(inputcast) %in% negcastvector),]

castsum <- rowSums(inputcast)

distlist <- lapply(c(1:25), function(i) {
    print(i)
    inputcast <- rrarefy(inputcast, sample=100000)
    inputcast <- inputcast[c(rowSums(inputcast) %in% 100000),]
    outdist <- vegdist(inputcast, method="bray", diag = TRUE, upper = TRUE)
    return(outdist)
})

castdist <- Reduce("+", distlist) / length(distlist)

distlist2 <- lapply(c(1:25), function(i) {
    print(i)
    inputcast <- rrarefy(inputcast, sample=100000)
    inputcast <- inputcast[c(rowSums(inputcast) %in% 100000),]
    # Make one without negative controls as well
    negcastvectorout <- as.character(datadisease[c(datadisease$V30 %in% "Negative"),1])
    negcastout <- inputcast[!c(rownames(inputcast) %in% negcastvectorout),]
    outdist <- vegdist(negcastout, method="bray", diag = TRUE, upper = TRUE)
    return(outdist)
})

castdistnoneg <- Reduce("+", distlist2) / length(distlist2)

# Make the distance matrix

ORD_NMDS <- metaMDS(castdist,k=2, trymax = 50)
ORD_FIT = data.frame(MDS1 = ORD_NMDS$points[,1], MDS2 = ORD_NMDS$points[,2])
ORD_FIT$SampleID <- rownames(ORD_FIT)
NMDS_AND_MAP <- merge(ORD_FIT, datadisease, by.x="SampleID", by.y="V2")

# Plot everything
plotnmds <- ggplot(NMDS_AND_MAP, aes(x=MDS1, y=MDS2, colour=V30)) +
    theme_classic() +
    geom_point() +
    scale_colour_manual(values = wes_palette("Darjeeling"), name = "Disease")

# Make the distance matrix without negative control data
ORD_NMDS <- metaMDS(castdistnoneg,k=2)
ORD_FIT = data.frame(MDS1 = ORD_NMDS$points[,1], MDS2 = ORD_NMDS$points[,2])
ORD_FIT$SampleID <- rownames(ORD_FIT)
NMDS_AND_MAP <- merge(ORD_FIT, datadisease, by.x="SampleID", by.y="V2")

# Plot everything
nonegplotnmds <- ggplot(NMDS_AND_MAP, aes(x=MDS1, y=MDS2, colour=V30)) +
    theme_classic() +
    geom_point() +
    scale_colour_manual(values = wes_palette("Darjeeling"), name = "Disease")

write("Dispersion with negative controls", stderr())
# Trying dispersion
catdistnames <- as.data.frame(as.matrix(castdist))
catdistnames$names <- row.names(catdistnames)
mergedist <- merge(catdistnames, datadisease, by.x="names", by.y="V2")
mod <- betadisper(castdist, mergedist$V30)
anovawithneg <- anova(mod)
anovawithnegsig <- c(anovawithneg[,"Pr(>F)"])[1]
mod.HSD <- TukeyHSD(mod)

anovastatwithneg <- c("anovawithneg", anovawithnegsig)

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
catdistnames <- as.data.frame(as.matrix(castdistnoneg))
catdistnames$names <- row.names(catdistnames)
mergedist <- merge(catdistnames, datadisease, by.x="names", by.y="V2")
mod <- betadisper(castdistnoneg, mergedist$V30)
anovawitouthneg <- anova(mod)
anovawitouthnegsig <- c(anovawitouthneg[,"Pr(>F)"])[1]
permutest(mod, pairwise = TRUE)
mod.HSD <- TukeyHSD(mod)

anovastatwithoutneg <- c("anovawithoutneg", anovawitouthnegsig)

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

inputshannon <- data.frame(rownames(negcast) , diversity(negcast, index="shannon"), sample=100000)
colnames(inputshannon) <- c("SampleID", "ShannonDiv")
mergedshannon <- merge(inputshannon, datadisease, by.x="SampleID", by.y="V2")
mergedshannon <- mergedshannon[!c(mergedshannon$V30 %in% "Negative"),]
shannonplot <- ggplot(mergedshannon, aes(x=V30, y=ShannonDiv, fill=V30)) +
    theme_classic() +
    theme(legend.position="none") +
    geom_boxplot(notch=FALSE) +
    xlab("Health Status") +
    ylab("Shannon Diversity") +
    scale_fill_manual(values = alpha(wes_palette("Darjeeling"), 0.75))
pairwise.wilcox.test(x=mergedshannon$ShannonDiv, g=mergedshannon$V30, p.adjust.method="bonferroni")

inputrich <- data.frame(rownames(negcast), rarefy(negcast, sample=100000))
colnames(inputrich) <- c("SampleID", "Rich")
mergedrich <- merge(inputrich, datadisease, by.x="SampleID", by.y="V2")
mergedrich <- mergedrich[!c(mergedrich$V30 %in% "Negative"),]
richplot <- ggplot(mergedrich, aes(x=V30, y=Rich, fill=V30)) +
    theme_classic() +
    theme(legend.position="none") +
    geom_boxplot(notch=FALSE) +
    xlab("Health Status") +
    ylab("Richness") +
    scale_fill_manual(values = alpha(wes_palette("Darjeeling"), 0.75))
pairwise.wilcox.test(x=mergedrich$Rich, g=mergedrich$V30, p.adjust.method="bonferroni")


alphaplot <- plot_grid(shannonplot, richplot, labels = c("C", "D"))
gridplot <- plot_grid(plotnmds, plotdiffs, labels = c("A", "B"))
gridplotnoneg <- plot_grid(nonegplotnmds, nonegplotdiffs, alphaplot, labels = c("A", "B"), ncol = 3)

pdf(file="./figures/diversity-beta-ogu.pdf", width=14, height=4)
    gridplotnoneg
dev.off()

pdf(file="./figures/diversity-beta-ogu-negative.pdf", width=6, height=4)
    plotdiffs
dev.off()

fortable <- data.frame(rbind(anovastatwithneg, anovastatwithoutneg))
colnames(fortable) <- c("Title", "Stat")
write.table(fortable, file = "./rtables/diversity.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
