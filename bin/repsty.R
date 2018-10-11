gcinfo(FALSE)

library("optparse")
library("ggplot2")
library("plotROC")
library("reshape2")
library("wesanderson")
library("vegan")
library("plyr")
library("dplyr")
library("caret")
library("cowplot")
library("matrixStats")
library("parallel")

minavgrelabund <- -1
minsamps <- 0

input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
tempteratelist <- read.delim("./data/repcycle/repclust.tsv", header = FALSE, sep = "\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

# Remove the bacteria-not-virus clusters
removaltable <- read.delim("./data/contigclustersidentity/BacteriaNotVirus.tsv", header = FALSE, sep = "\t")
removaltable$V1 <- gsub("^", "Cluster_", removaltable$V1, perl = TRUE)
# Clean input
input <- input[!c(input$V1 %in% removaltable$V1),]


# Rarefy input table
minimumsubsample <- 1000000
inputcast <- dcast(input, V1 ~ V2)
inputcast[is.na(inputcast)] <- 0
inputcast[,-1] <-round(inputcast[,-1],0)
row.names(inputcast) <- inputcast[,1]
inputcast <- inputcast[,-1]
inputcast <- as.data.frame(t(inputcast))
counter <- 1
rarefunction <- function(y) {
	if(sum(y) >= 1e9) {
		return(rarefybig(y, minimumsubsample))
	} else {
		return(rrarefy(y, minimumsubsample))
	}
}
rareoutput <- lapply(c(1:length(inputcast[,1])), function(i) {
	write(counter, stderr())
	counter <<- counter + 1; 
	subsetdf <- inputcast[i,]
	if(sum(subsetdf) >= minimumsubsample) {
		return(rarefunction(subsetdf))
	}
})
rareoutputbind <- as.data.frame(do.call(rbind, rareoutput))
length(rareoutputbind[,1])

inputmelt <- melt(as.matrix(rareoutputbind))
colnames(inputmelt) <- c("V2", "V1", "sum")
inputmelt$sum <- ifelse(inputmelt$sum >= (minimumsubsample)^-1, inputmelt$sum, 0)

inputrelabund <- data.frame(inputmelt %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
relabundcast <- dcast(inputrelabund, V2 ~ V1, value.var = "RelAbund")
relabundcast[is.na(relabundcast)] <- 0
row.names(relabundcast) <- relabundcast$V2
relabundcast <- relabundcast[,c(TRUE, colMedians(as.matrix(relabundcast[,-1])) > minavgrelabund)]

castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
castmerge <- castmerge[,-1]

# Filter by presence/absence
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
# abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
abmelt <- melt(abssubset)

abmelt$repsty <- ifelse(abmelt$variable %in% tempteratelist$V1, "Lysogenic", "Lytic")
abply <- ddply(abmelt, c("V30", "V22", "repsty"), summarize, sum = sum(value))
abply$sample <- gsub("\\D+", "", abply$V22, perl = TRUE)

lysoplot <- ggplot(abply[c(abply$repsty %in% "Lysogenic"),], aes(x = V30, y = sum)) +
	theme_classic() +
	geom_boxplot(notch = TRUE, fill = "grey") +
	ylab("Lysogenic Relative Abundance") +
	xlab("")

pairwise.wilcox.test(abply$sum, abply$V30, p.adjust = "bonferroni")

pdf("./figures/lysogenic-relabund.pdf", height = 4, width = 6)
	lysoplot
dev.off()
