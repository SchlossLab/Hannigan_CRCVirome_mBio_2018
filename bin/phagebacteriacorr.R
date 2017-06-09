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
library("cowplot")

minavgrelabund <- 0
minsamps <- 30

caretmodel <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
	repeats = 5,
	number=5,
	classProbs = TRUE,
	summaryFunction = twoClassSummary,
	savePredictions = TRUE)
  model <- train(V30~., data=x, trControl=fitControl, method="rf", metric="ROC", tuneLength=5)
  return(model)
}

impcalcs <- function(x, iterationcount = 10) {
	avgimportance <- lapply(c(1:iterationcount), function(i) GetAverageImportance(x, i))
	avgimportancedf <- ldply(avgimportance, data.frame)
	import <- ddply(avgimportancedf, c("categories"), summarize, mean = mean(Overall))
	importaverage <- merge(avgimportancedf, import, by = "categories")

	virustax <- virustax[,c(1,3,7)]

	importaverage <- merge(importaverage, virustax, by.x = "categories", by.y = "V1", all = TRUE)
	importaverage$V7 <- as.character(importaverage$V7)
	importaverage <- importaverage[!c(importaverage$Overall %in% NA),]
	importaverage[is.na(importaverage)] <- "Unknown"
	importaverage <- importaverage[order(importaverage$mean, decreasing = TRUE),]
	importaverage$categories <- factor(importaverage$categories, levels = importaverage$categories)
	importaverage$V7 <- factor(importaverage$V7, levels = importaverage$V7)
	
	return(importaverage)
}

GetAverageImportance <- function(x, y) {
	write(y, stderr())
	functionmodel <- caretmodel(x)
	resultvardf <- data.frame(varImp(functionmodel$finalModel))
	resultvardf$categories <- rownames(resultvardf)
	resultvardf$categories <- factor(resultvardf$categories, levels = resultvardf$categories)
	resultvardf$iteration <- y
	return(resultvardf)
}

##########
# Virome #
##########

# Run analysis
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]
taxonomy <- read.delim("./data/mothur16S/final.taxonomy", header = TRUE, sep = "\t")
virustax <- read.delim("./data/contigclustersidentity/clustax.tsv", header = FALSE, sep = "\t")
# Format taxonomy table
taxonomy$Taxonomy <- sub(".+\\;(.+)\\(\\d+\\)\\;$", "\\1", taxonomy$Taxonomy, perl=TRUE)
virustax$V1 <- sub("^", "Cluster_", virustax$V1, perl = TRUE)

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

castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
castmerge <- castmerge[,-1]

# Filter by presence/absence
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

############
# Bacteria #
############
inputbacteria <- read.delim("./data/mothur16S/final.shared", header=TRUE, sep="\t")

counter <- 1
rareoutput <- lapply(c(1:length(inputbacteria[,1])), function(i) {
	write(counter, stderr())
	counter <<- counter + 1; 
	subsetdf <- inputbacteria[i,-c(1:3)]
	names <- inputbacteria[i,c(1:3)]
	if(sum(subsetdf) >= 10000) {
		rareoutput <- rrarefy(subsetdf, 10000)
		y <- cbind(names, rareoutput)
	}
	return(y)
})
rareoutputbacteria <- do.call(rbind, rareoutput)

inputbacteria$Group <- sub("^(\\D)(\\d)$","\\10\\2", sub("(.)\\D+(\\d+)", "\\1\\2", inputbacteria$Group, perl=TRUE))
inputbacteria$Group <- as.factor(inputbacteria$Group)
# Calculate as relative abundance
inputbacteriamelt <- melt(inputbacteria[-c(1,3)])
# Get relative abundance
inputbacteriarelabund <- data.frame(inputbacteriamelt %>% group_by(Group) %>% mutate(RelAbund = 100 * value / sum(value)))
relabundcast <- dcast(inputbacteriarelabund, Group ~ variable, value.var = "RelAbund")
# Fix the metadata for this case without duplicate MG IDs
datadiseasesub <- datadisease[,2:3]
datadiseasesub <- datadiseasesub[!duplicated(datadiseasesub),]

# Add the disease classes
relabundclasses <- merge(relabundcast, datadiseasesub, by.x="Group", by.y="V22")

# Compare healthy to cancer
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Adenoma"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]

################
# Correlations #
################
virus <- abssubset[,-(length(abssubset)-1)]
bacteria <- pasubset[,-(length(pasubset))]
virus <- virus[c(virus$V22 %in% bacteria$Group),]
bacteria <- bacteria[c(bacteria$Group %in% virus$V22),]
virus <- virus[,-length(virus)]
bacteria <- bacteria[,-1]

cordf <- melt(cor(virus, bacteria))

# Iterate to get average importance plot
importancevalsvirus <- impcalcs(absmissingid)
importancevalsbacteria <- impcalcs(pasubset)

virusclean <- unique(importancevalsvirus[, c(1,4)])
colnames(virusclean) <- c("categories", "virusmean")
bacteriaclean <- unique(importancevalsbacteria[, c(1,4)])
colnames(bacteriaclean) <- c("categories", "bacteriamean")
cordmerge <- merge(cordf, virusclean, by.x = "Var1", by.y = "categories")
cordmerge <- merge(cordmerge, bacteriaclean, by.x = "Var2", by.y = "categories")
ordermerge <- cordmerge[order(cordmerge[4], decreasing = FALSE),]
ordermerge <- ordermerge[order(ordermerge[5], decreasing = TRUE),]
ordermerge$Var1 <- factor(ordermerge$Var1, levels = ordermerge$Var1)
ordermerge$Var2 <- factor(ordermerge$Var2, levels = ordermerge$Var2)

# Top left has the most important OGUs and OTUs
allplot <- ggplot(ordermerge, aes(x = Var2, y = Var1, fill = value)) +
	theme_bw() +
	geom_tile(color = "white") +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.title.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.y=element_blank()
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1), name = "Correlation") +
	ylab("Viral OGUs") +
	xlab("Bacterial OTUs") +
	geom_rect(xmax = 10, ymax = length(unique(ordermerge$Var1)), xmin = 1, ymin = length(unique(ordermerge$Var1)) - 10, color = "black", alpha = 0)

topdf <- ordermerge[c(ordermerge$Var2 %in% bacteriaclean[1:10,"categories"] & ordermerge$Var1 %in% virusclean[1:10,"categories"]),]
fcordf <- ordermerge[c(ordermerge$value > 0.50 | ordermerge$value < -0.50),]
topvirus <- ordermerge[c(ordermerge$Var1 %in% virusclean[1:5,"categories"] & ordermerge$value > 0.5),]
topbac <- ordermerge[c(ordermerge$Var2 %in% bacteriaclean[1:5,"categories"] & ordermerge$value > 0.5),]

# Try adding a bar plot
vimp <- ggplot(ordermerge, aes(x = Var1, y = virusmean)) +
	theme_classic() +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank()
	) +
	geom_bar(stat = "identity", position = "dodge", width=1) +
	coord_flip() +
	scale_y_reverse() +
	xlab("Viral OVUs")

bimp <- ggplot(ordermerge, aes(x = Var2, y = bacteriamean)) +
	theme_classic() +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank()
	) +
	geom_bar(stat = "identity", position = "dodge", width=1) +
	xlab("Bacterial OTUs") +
	scale_y_reverse()

allleg <- get_legend(allplot)
anl <- allplot + theme(legend.position = "none")
ap1 <- plot_grid(vimp, anl, NULL, bimp, nrow = 2, align = "hv", rel_widths = c(2, 10), rel_heights = c(10, 2))
allplotf <- plot_grid(ap1, allleg, nrow = 1, rel_widths = c(4, 1))

hcplot <- ggplot(fcordf, aes(x = Var2, y = Var1, fill = value)) +
	theme_classic() +
	geom_tile(color = "white") +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 1)
		legend.position = "none"
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	ylab("Viral OVUs") +
	xlab("Bacterial OTUs")

topplot <- ggplot(topdf, aes(x = gsub("Otu0+", "OTU ", Var2), y = Var1, fill = value)) +
	theme_classic() +
	geom_tile(color = "white") +
	theme(
		# axis.text.x=element_blank(),
		# axis.ticks.x=element_blank(),
		# axis.text.y=element_blank(),
		# axis.ticks.y=element_blank()
		axis.text.x = element_text(angle = 90, hjust = 1),
		legend.position = "none"
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	ylab("Viral OVUs") +
	xlab("Bacterial OTUs") +
	scale_x_discrete(labels = gsub("Otu0+", "OTU ", unique(topdf$Var2))) +
	scale_y_discrete(labels = gsub("_", " ", unique(topdf$Var1)))

ggplot(topvirus, aes(x = Var2, y = Var1, fill = value)) +
	theme_classic() +
	geom_tile(color = "white") +
	theme(
		# axis.text.x=element_blank(),
		# axis.ticks.x=element_blank(),
		# axis.text.y=element_blank(),
		# axis.ticks.y=element_blank()
		axis.text.x = element_text(angle = 90, hjust = 1)
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	ylab("Viral OGUs") +
	xlab("Bacterial OTUs")

ggplot(topbac, aes(x = Var1, y = Var2, fill = value)) +
	theme_classic() +
	geom_tile(color = "white") +
	theme(
		# axis.text.x=element_blank(),
		# axis.ticks.x=element_blank(),
		# axis.text.y=element_blank(),
		# axis.ticks.y=element_blank()
		axis.text.x = element_text(angle = 90, hjust = 1)
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	xlab("Viral OGUs") +
	ylab("Bacterial OTUs")

histplot <- ggplot(ordermerge, aes(x = value)) +
	theme_classic() +
	geom_histogram() +
	scale_x_continuous(limits = c(-1, 1)) +
	xlab("Pearson Correlation") +
	ylab("Frequency")

subsets <- plot_grid(topplot, histplot, ncol = 1, labels = c("B", "C"), align = "v")
finalplot <- plot_grid(allplotf, subsets, labels = "A")


pdf("./figures/phage-bacteria-cor.pdf", height = 5, width = 12)
	finalplot
dev.off()


