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

plotimportance <- function(x, iterationcount = 25, topcount = 10) {
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
	
	binlength <- c(1:topcount) + 0.5

	numberto <- topcount * iterationcount - 1
	dfplot <- importaverage[(length(importaverage[,1])-numberto):(length(importaverage[,1])),]
	
	importanceplot <- ggplot(dfplot, aes(x=categories, y=Overall)) +
	  theme_classic() +
	  theme(
	    axis.line.x = element_line(colour = "black"),
	    axis.line.y = element_line(colour = "black")
	  ) +
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.25, stackratio = 0.5) +
	  xlab("Virus Identity") +
	  ylab("Mean Accuracy Decrease") +
	  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.75) +
	  geom_vline(xintercept=binlength,color="grey") +
	  scale_x_discrete(labels=dfplot[c(0:(topcount - 1))*iterationcount+1,"V7"]) +
	  coord_flip()

	  return(importanceplot)
}

pbi <- function(x, iterationcount = 25, topcount = 10) {
	avgimportance <- lapply(c(1:iterationcount), function(i) GetAverageImportance(x, i))
	avgimportancedf <- ldply(avgimportance, data.frame)
	import <- ddply(avgimportancedf, c("categories"), summarize, mean = mean(Overall))
	importaverage <- merge(avgimportancedf, import, by = "categories")
	importaverage <- merge(importaverage, taxonomy, by.x = "categories", by.y = "OTU")
	importaverage <- importaverage[order(importaverage$mean, decreasing = FALSE),]
	importaverage$categories <- factor(importaverage$categories, levels = importaverage$categories)
	importaverage$Taxonomy <- factor(importaverage$Taxonomy, levels = importaverage$Taxonomy)
	
	binlength <- c(1:topcount) + 0.5

	numberto <- topcount * iterationcount - 1
	dfplot <- importaverage[(length(importaverage[,1])-numberto):(length(importaverage[,1])),]
	
	importanceplot <- ggplot(dfplot, aes(x=categories, y=Overall)) +
	  theme_classic() +
	  theme(
	    axis.line.x = element_line(colour = "black"),
	    axis.line.y = element_line(colour = "black")
	  ) +
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.25, stackratio = 0.5) +
	  xlab("Bacteria Identity") +
	  ylab("Mean Accuracy Decrease") +
	  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.75) +
	  geom_vline(xintercept=binlength,color="grey") +
	  scale_x_discrete(labels=dfplot[c(0:(topcount - 1))*iterationcount+1,"Taxonomy"]) +
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
virustax <- read.delim("./data/contigclustersidentity/clustax.tsv", header = FALSE, sep = "\t")
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
melttry <- merge(melttry, virustax[,c(1,3,7)], by.x = "variable", by.y = "V1", all = TRUE)
melttry$V7 <- as.character(melttry$V7)
melttry <- melttry[!c(melttry$value %in% NA),]
melttry[is.na(melttry)] <- "Unknown"
binlength <- c(1:5) + 0.5
importanceabundance <- ggplot(melttry, aes(x = factor(variable), y = (value + 1e-06), fill = factor(V30))) +
	theme_classic() +
	geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.75, stackratio = 0.5) +
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
	scale_x_discrete(labels=melttry[c(0:(72 - 1))*72+1,"V7"]) +
	scale_fill_manual(values = wes_palette("Royal1")[c(1,2,4)], name = "Disease")

#############################
# Bacteria Prediction Model #
#############################
inputbacteria <- read.delim("./data/mothur16S/final.shared", header=TRUE, sep="\t")
taxonomy <- read.delim("./data/mothur16S/final.taxonomy", header = TRUE, sep = "\t")
taxonomy$Taxonomy <- sub(".+\\;(.+)\\(\\d+\\)\\;$", "\\1", taxonomy$Taxonomy, perl=TRUE)

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
importanceplotbac <- pbi(pasubsetmissing)



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

statsig <- wilcox.test(megatron$highAUC ~ megatron$class)

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
	ylab("Area Under the Curve") +
	geom_segment(x = 1, xend = 2, y = 0.8, yend = 0.8) +
	annotate("text", x = 1.5, y = 0.805, label = paste("p-value = ", signif(statsig$p.value, digits = 2), sep = ""), size = 4) +
	ylim(NA, 0.81)

aucstat <- ddply(megatron, "class", summarize, meanauc = mean(highAUC))


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

write.table(aucstat, file = "./rtables/threeclass-auc.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(statsig$p.value, file = "./rtables/threeclass-sig.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
