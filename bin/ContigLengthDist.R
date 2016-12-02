# ContigLengthDist.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

library("optparse")
library("ggplot2")
library("cowplot")
library("plyr")

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

input <- read.delim("./data/VirusContigLength.tsv", head=FALSE, sep="\t")
colnames(input) <- c("ContigID", "Length")

clusters <- read.delim("./data/ContigClustersVirus/clustering_gt1000.csv", head = FALSE, sep = ",")
colnames(clusters) <- c("ContigID", "Cluster")

mergeddf <- merge(input, clusters, by = "ContigID")

lengthscatter <- ggplot(mergeddf, aes(x = Cluster, y = Length)) +
    theme_classic() +
    geom_point() +
    scale_y_log10() +
    xlab("Viral Cluster ID") +
    coord_flip()

contigspercluster <- count(mergeddf$Cluster)

clusterdens <- ggplot(contigspercluster, aes(x = freq)) +
    theme_classic() +
    geom_density(fill = "gray", alpha=0.5) +
    xlab("Viral Contig Count Frequency per Cluster") +
    scale_x_log10()

clustergrid <- plot_grid(lengthscatter, clusterdens, labels = c("A", "B"))

# Bacteria

inputbac <- read.delim("./data/BacteriaContigLength.tsv", head=FALSE, sep="\t")
colnames(inputbac) <- c("ContigID", "Length")

clustersbac <- read.delim("./data/ContigClustersBacteria/clustering_gt2000.csv", head = FALSE, sep = ",")
colnames(clustersbac) <- c("ContigID", "Cluster")

mergeddfbac <- merge(inputbac, clustersbac, by = "ContigID")

lengthscatterbac <- ggplot(mergeddfbac, aes(x = Cluster, y = Length)) +
    theme_classic() +
    geom_point() +
    scale_y_log10() +
    xlab("Bacteria Cluster ID") +
    coord_flip()

contigsperclusterbac <- count(mergeddfbac$Cluster)

clusterdensbac <- ggplot(contigsperclusterbac, aes(x = freq)) +
    theme_classic() +
    geom_density(fill = "gray", alpha=0.5) +
    xlab("Bacterial Contig Count Frequency per Cluster") +
    scale_x_log10()

clustergridbacteria <- plot_grid(lengthscatterbac, clusterdensbac, labels = c("C", "D"))

totalplot <- plot_grid(clustergrid, clustergridbacteria, ncol = 1)

pdf(file = "./figures/contiglengthbycluster.pdf", width = 10, height = 10)
    totalplot
dev.off()
