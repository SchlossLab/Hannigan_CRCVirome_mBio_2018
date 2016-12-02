# ContigLengthDist.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

library("optparse")
library("ggplot2")
library("cowplot")
library("plyr")
library("hexbin")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Contig length table.",
    metavar = "character"),
  make_option(c("-c", "--clusters"),
    type = "character",
    default = NULL,
    help = "Contig cluster table.",
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

input <- read.delim("./data/VirusContigLength.tsv", head = FALSE, sep = "\t")
colnames(input) <- c("ContigID", "Length")

clusters <- read.delim("./data/ContigClustersVirus/clustering_gt1000.csv", head = FALSE, sep = ",")
colnames(clusters) <- c("ContigID", "Cluster")

abundance <- read.delim("./data/ContigRelAbundForGraphVirus.tsv", head = FALSE, sep = "\t")
abundancesum <- ddply(abundance, "V1", summarize, sum = sum(V2))

mergedabund <- merge(input, abundancesum, by.x = "ContigID", by.y = "V1")

viruscontigplot <- ggplot(mergedabund, aes(x=Length, y=sum)) +
    theme_classic() +
    stat_binhex(bins=100) +
    scale_fill_gradientn(colours=c("steelblue4", "steelblue3", "steelblue2","steelblue1","slategray1")) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("Contig Length (log bp)") +
    ylab("Sequences Mapping to Contigs (log)")

mergedclusters <- merge(mergedabund, clusters, by = "ContigID")

clusterstats <- ddply(mergedclusters, "Cluster", summarize, sumcount = sum(sum), avglength = mean(Length))

virusclusterplot <- ggplot(clusterstats, aes(x=avglength, y=sumcount)) +
    theme_classic() +
    geom_point() +
    scale_fill_gradientn(colours=c("steelblue4", "steelblue3", "steelblue2","steelblue1","slategray1")) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("Average Cluster Length (log bp)") +
    ylab("Sequences Mapping to Clusters (log)")

# Bacteria

inputbac <- read.delim("./data/BacteriaContigLength.tsv", head = FALSE, sep = "\t")
colnames(inputbac) <- c("ContigID", "Length")

clustersbac <- read.delim("./data/ContigClustersBacteria/clustering_gt2000.csv", head = FALSE, sep = ",")
colnames(clustersbac) <- c("ContigID", "Cluster")

abundancebac <- read.delim("./data/ContigRelAbundForGraphBacteria.tsv", head = FALSE, sep = "\t")
abundancebacsum <- ddply(abundancebac, "V1", summarize, sum = sum(V2))

mergedabundbac <- merge(inputbac, abundancebacsum, by.x = "ContigID", by.y = "V1")

bacteriacontigplot <- ggplot(mergedabundbac, aes(x=Length, y=sum)) +
    theme_classic() +
    stat_binhex(bins=100) +
    scale_fill_gradientn(colours=c("steelblue4", "steelblue3", "steelblue2","steelblue1","slategray1")) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("Contig Length (log bp)") +
    ylab("Sequences Mapping to Contigs (log)")

mergedclustersbac <- merge(mergedabundbac, clustersbac, by = "ContigID")

clustersbactats <- ddply(mergedclustersbac, "Cluster", summarize, sumcount = sum(sum), avglength = mean(Length))

bacteriaclusterplot <- ggplot(clustersbactats, aes(x=avglength, y=sumcount)) +
    theme_classic() +
    geom_point() +
    scale_fill_gradientn(colours=c("steelblue4", "steelblue3", "steelblue2","steelblue1","slategray1")) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("Average Cluster Length (log bp)") +
    ylab("Sequences Mapping to Clusters (log)")

totalplot <- plot_grid(viruscontigplot, virusclusterplot, bacteriacontigplot, bacteriaclusterplot, ncol = 2, labels = c("A", "B", "C", "D"))

pdf("./figures/qc-contiglengthcoverage.pdf", height = 8, width = 8)
    totalplot
dev.off()
