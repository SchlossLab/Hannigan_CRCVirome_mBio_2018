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

impcalcs <- function(x, iterationcount = 5) {
	avgimportance <- lapply(c(1:iterationcount), function(i) GetAverageImportance(x, i))
	avgimportancedf <- ldply(avgimportance, data.frame)
	import <- ddply(avgimportancedf, c("categories"), summarize, mean = mean(Overall))
	importaverage <- merge(avgimportancedf, import, by = "categories")

	virustax <- virustax[,c(1,3,7)]

	importaverage <- merge(importaverage, virustax, by.x = "categories", by.y = "V1", all = TRUE)
	importaverage$V7 <- as.character(importaverage$V7)
	importaverage <- importaverage[!c(importaverage$Overall %in% NA),]
	importaverage[is.na(importaverage)] <- "Unknown"
	importaverage <- importaverage[order(importaverage$mean, decreasing = FALSE),]
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
fcordf <- cordf[c(cordf$value > 0.75 | cordf$value < -0.75),]

hm <- ggplot(cordf, aes(x = Var2, y = Var1, fill = value)) +
	theme_classic() +
	geom_tile(color = "white") +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank()
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	ylab("Viral OGUs") +
	xlab("Bacterial OTUs")

ggplot(fcordf, aes(x = Var2, y = Var1, fill = value)) +
	theme_classic() +
	geom_tile(color = "white") +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	ylab("Viral OGUs") +
	xlab("Bacterial OTUs") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

outmodel <- caretmodel(absmissingid)
outmodel

# Iterate to get average importance plot
importancevals <- impcalcs(absmissingid)
importancevals <- impcalcs(pasubset)

clusterorder <- data.frame(cluster = rev(unique(importancevals$categories)), order = c(1:length(unique(importancevals$categories))))
cmer <- merge(clusterorder, fcordf, by.x = "cluster", by.y = "Var1")
cordforder <- cmer[order(cmer$order, decreasing = TRUE),]
cordforder$cluster <- factor(cordforder$cluster, levels = cordforder$cluster)
cordforder$Var2 <- factor(cordforder$Var2, levels = cordforder$Var2)
cordforder[c(cordforder$cluster %in% "Cluster_115"),]

ggplot(cordforder, aes(x = Var2, y = cluster, fill = value)) +
	theme_bw() +
	geom_tile(color = "white") +
	theme(
	# 	axis.text.x=element_blank(),
	# 	axis.ticks.x=element_blank(),
	# 	axis.text.y=element_blank(),
	# 	axis.ticks.y=element_blank(),
	axis.text.x = element_text(angle = 90, hjust = 1)
	) +
	scale_fill_gradient2(low = "tomato4", mid = "white", high = "steelblue4", limits=c(-1,1)) +
	ylab("Viral OGUs") +
	xlab("Bacterial OTUs")


pdf("./figures/phage-bacteria-cor.pdf", height = 4, width = 6)
	hm
dev.off()


