library("optparse")
library("ggplot2")
library("plotROC")
library("reshape2")
library("wesanderson")
library("vegan")
library("plyr")
library("dplyr")
library("caret")

processmothurotus <- function(x, removeclass="none") {
	# Clean up disease classes
	x$V30 <- sub("[0-9]+", "", x$Group, perl=TRUE)
	if (removeclass != "none") {
		x <- x[-which(x$V30 %in% removeclass),]
	}
	return(x)
}

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

##########################
# Virus Prediction Model #
##########################

input <- read.delim("./data/ClusteredContigAbund.tsv", header=TRUE, sep="\t")
head(input)
datadisease <- read.delim("./data/metadata/NexteraXT003Map.tsv", header=FALSE, sep="\t")[,c(2,30)]
head(datadisease)

inputrelabund <- data.frame(input %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
relabundcast <- dcast(inputrelabund, V2 ~ V1, value.var = "RelAbund")
relabundcast[is.na(relabundcast)] <- 0
row.names(relabundcast) <- relabundcast$V2
castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
castmerge <- castmerge[,-1]

# # Filter by median relative abundance
# lowabundgone <- castmerge[,c(sapply(castmerge[-length(castmerge)], median) > 0.01,TRUE)]
# lowsubset <- lowabundgone[!c(lowabundgone$V30 %in% "Negative"),]
# lowsubset <- lowsubset[!c(lowsubset$V30 %in% "Adenoma"),]
# lowsubset$V30 <- factor(lowsubset$V30)

# Filter by presence/absence
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset[,c(colSums(abssubset != 0) > 30)]

outmodel <- caretmodel(abssubset)
outmodel

plot(outmodel)

ggplot(outmodel$pred, aes(d = obs, m = Healthy)) +
	geom_roc(n.cuts = 0, color = wes_palette("Royal1")[2]) +
	theme_classic() +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype=2, colour=wes_palette("Royal1")[1]) +
	ylab("Sensitivity") +
	xlab(paste("Inverse Specificity"))

#############################
# Bacteria Prediction Model #
#############################
input <- read.delim("./data/mothur16S/final.shared", header=TRUE, sep="\t")
# Calculate as relative abundance
inputmelt <- melt(input[-c(1,3)])
# Get relative abundance
inputrelabund <- data.frame(inputmelt %>% group_by(Group) %>% mutate(RelAbund = 100 * value / sum(value)))
relabundcast <- dcast(inputrelabund, Group ~ variable, value.var = "RelAbund")

# Compare healthy to cancer
inputnoademona <- processmothurotus(relabundcast, removeclass="Adenoma")

# Filter by presence/absence
pasubset <- inputnoademona[,c(colSums(inputnoademona != 0) > 30)]
subsetmodel <- caretmodel(pasubset)
subsetmodel
plot(subsetmodel)

# Plot the ROC curve
ggplot(subsetmodel$pred, aes(d = obs, m = Healthy)) +
	geom_roc(n.cuts = 0, color = wes_palette("Royal1")[2]) +
	theme_classic() +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype=2, colour=wes_palette("Royal1")[1]) +
	ylab("Sensitivity") +
	xlab(paste("Inverse Specificity"))

###############################
# Compare Bacteria and Virus  #
###############################
subsetmodel$pred$class <- "Bacteria"
outmodel$pred$class <- "Virus"
boundmodel <- rbind(subsetmodel$pred, outmodel$pred)

# Plot the ROC curve
boundplot <- ggplot(boundmodel, aes(d = obs, m = Healthy, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1"))
boundplot

pdf("./figures/predmodel-viromebacteria.pdf", height = 5, width = 6)
	boundplot
dev.off()

