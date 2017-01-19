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

minavgrelabund <- 0
minsamps <- 30

# Check out this reference for some more information on multiclass tuning
# https://github.com/topepo/caret/issues/107

caretmodel <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
	repeats = 5,
	number=5,
	classProbs = TRUE,
	summaryFunction = multiClassSummary,
	savePredictions = TRUE)
  model <- train(V30~., data=x, trControl=fitControl, method="rf", metric="Mean_AUC", tuneLength=5)
  return(model)
}

makemeanroc <- function(x, keep = "Cancer", remove = "Adenoma", one = "Healthy") {
	plotx <- x$pred
	plotx[plotx==remove]<-keep
	plotx$all <- plotx[, keep] + plotx[, remove]
	plotx$one <- plotx[, one]
	plotx <- plotx[ , !(names(plotx) %in% c(remove, keep, one))]
	plotx$obs <- factor(plotx$obs)
	plotx$pred <- factor(plotx$pred)
	plotx$class <- one
	levels(plotx$obs)[match(keep,levels(plotx$obs))] <- "All"
	levels(plotx$obs)[match(one,levels(plotx$obs))] <- "One"
	return(plotx)
}

threeclassmeanroc <- function(x, classone = "Cancer", classtwo = "Adenoma", classthree = "Healthy") {
	HealthyVsAll <- makemeanroc(x, keep = "Cancer", remove = "Adenoma", one = "Healthy")
	CancerVsAll <- makemeanroc(x, keep = "Healthy", remove = "Adenoma", one = "Cancer")
	AdenomaVsAll <- makemeanroc(x, keep = "Healthy", remove = "Cancer", one = "Adenoma")
	boundplot <- rbind(HealthyVsAll, CancerVsAll, AdenomaVsAll)
	return(boundplot)
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

plotimportance <- function(x, iterationcount = 10, topcount = 10) {
	avgimportance <- lapply(c(1:iterationcount), function(i) GetAverageImportance(x, i))
	avgimportancedf <- ldply(avgimportance, data.frame)
	import <- ddply(avgimportancedf, c("categories"), summarize, mean = mean(Overall))
	importaverage <- merge(avgimportancedf, import, by = "categories")
	importaverage <- importaverage[order(importaverage$mean, decreasing = FALSE),]
	importaverage$categories <- factor(importaverage$categories, levels = importaverage$categories)
	
	binlength <- c(1:topcount) + 0.5
	
	importanceplot <- ggplot(importaverage[(length(importaverage[,1])-99):(length(importaverage[,1])),], aes(x=categories, y=Overall)) +
	  theme_classic() +
	  theme(
	    axis.line.x = element_line(colour = "black"),
	    axis.line.y = element_line(colour = "black")
	  ) +
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", binwidth = 0.05, stackdir = "center") +
	  xlab("Categories") +
	  ylab("Mean Decrease in Accuracy") +
	  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5) +
	  geom_vline(xintercept=binlength,color="grey") +
	  coord_flip()

	  return(importanceplot)
}

GetAverageAUC <- function(x, y) {
	write(y, stderr())
	functionmodel <- caretmodel(x)
	highAUC <- functionmodel$results$Mean_AUC[order(functionmodel$results$Mean_AUC, decreasing = TRUE)[1]]
	resultdf <- data.frame(y,highAUC)
	return(resultdf)
}


##########################
# Virus Prediction Model #
##########################
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

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
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

length(absmissingid[,1])

outmodel <- caretmodel(absmissingid)
outmodel

modelmeanroc <- threeclassmeanroc(outmodel)

# Plot the ROC curve
meanrocvirome <- ggplot(modelmeanroc, aes(d = obs, m = one, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1")[c(1,2,4)], name = "Disease")


# Get the variable importance
importanceplotvirus <- plotimportance(absmissingid)

# Get relative abundance of the important OGUs between disease states
vardf <- data.frame(varImp(outmodel$finalModel))
vardf$categories <- rownames(vardf)
vardf <- vardf[order(vardf$Overall, decreasing = FALSE),]
vardf$categories <- factor(vardf$categories, levels = vardf$categories)
varlength <- length(vardf$categories)
pullvalues <- as.character(vardf[c((varlength-5):varlength),c("categories")])
exampletry <- absmissingid[,c(pullvalues,"V30")]
melttry <- melt(exampletry)
binlength <- c(1:5) + 0.5
importanceabundance <- ggplot(melttry, aes(x = factor(variable), y = (value + 1e-06), fill = factor(V30))) +
	theme_classic() +
	geom_dotplot(binaxis = "y", binwidth = 0.2, stackdir = "center", position = "dodge") +
	stat_summary_bin(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5, position = "dodge") +
	scale_y_log10() +
	geom_vline(xintercept=binlength,color="grey") +
	  theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	xlab("Operational Genomic Units") +
	ylab("Relative Abundance") +
	scale_fill_manual(values = wes_palette("Royal1")[c(1,2,4)], name = "Disease")

#############################
# Bacteria Prediction Model #
#############################
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
# relabundcast <- relabundcast[,c(TRUE, colMedians(as.matrix(relabundcast[,-1])) > minavgrelabund)]

# Add the disease classes
relabundclasses <- merge(relabundcast, datadiseasesub, by.x="Group", by.y="V22")
relabundclasses$V30 <- droplevels(relabundclasses$V30)

# Filter by presence/absence
pasubset <- relabundclasses[,c(colSums(relabundclasses != 0) > minsamps)]
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

subsetmodel <- caretmodel(pasubsetmissing)
subsetmodel

modelmeanrocbacteria <- threeclassmeanroc(subsetmodel)

# Plot the ROC curve
meanrocbacteria <- ggplot(modelmeanrocbacteria, aes(d = obs, m = one, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1")[c(1,2,4)], name = "Disease")


# Get the variable importance
# Get the variable importance
importanceplotbac <- plotimportance(pasubsetmissing)



############################
# Quantify AUC Differences #
############################
iterationcount <- 10

viromeauc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
viromeaucdf <- ldply(viromeauc, data.frame)
viromeaucdf$class <- "Virus"

bacteriaauc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(pasubsetmissing, i))
bacteriaaucdf <- ldply(bacteriaauc, data.frame)
bacteriaaucdf$class <- "Bacteria"

megatron <- rbind(viromeaucdf, bacteriaaucdf)

wilcox.test(megatron$highAUC ~ megatron$class)

binlength <- c(1:2) + 0.5

auccompareplot <- ggplot(megatron, aes(x = class, y = highAUC, fill = class)) +
	theme_classic() +
	theme(legend.position="none") +
	geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", binwidth = 0.0025, stackdir = "center") +
	stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5) +
	geom_vline(xintercept=binlength,color="grey") +
	scale_fill_manual(values = wes_palette("Royal1")) +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	xlab("") +
	ylab("Area Under the Curve")


###############################
# Compare Bacteria and Virus  #
###############################
importanceplots <- plot_grid(importanceplotvirus, importanceplotbac, labels = c("D", "E"), ncol = 2)
boundplot <- plot_grid(meanrocvirome, meanrocbacteria, auccompareplot, labels = c("A", "B", "C"), rel_widths = c(2, 2, 1), ncol = 3)
lowerplot <- plot_grid(importanceabundance, labels = c("F"))
topbottomplot <- plot_grid(boundplot, importanceplots, lowerplot, rel_heights = c(3, 2, 3), ncol = 1)

pdf("./figures/predmodel-threewayclassification.pdf", height = 10, width = 12)
	topbottomplot
dev.off()
