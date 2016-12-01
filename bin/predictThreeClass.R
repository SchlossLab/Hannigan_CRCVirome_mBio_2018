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

# Check out this reference for some more information on multiclass tuning
# https://github.com/topepo/caret/issues/107

caretmodel <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
	repeats = 5,
	number=5,
	classProbs = TRUE,
	summaryFunction = multiClassSummary,
	savePredictions = TRUE)
  model <- train(V30~., data=x, trControl=fitControl, method="rf", metric="Mean_ROC", tuneLength=5)
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

##########################
# Virus Prediction Model #
##########################
input <- read.delim("./data/ClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

inputrelabund <- data.frame(input %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
relabundcast <- dcast(inputrelabund, V2 ~ V1, value.var = "RelAbund")
relabundcast[is.na(relabundcast)] <- 0
row.names(relabundcast) <- relabundcast$V2
castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
castmerge <- castmerge[,-1]

# Filter by presence/absence
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > 30)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

outmodel <- caretmodel(absmissingid)
outmodel

makemeanroc(outmodel)

modelmeanroc <- threeclassmeanroc(outmodel)

# Plot the ROC curve
meanrocvirome <- ggplot(modelmeanroc, aes(d = obs, m = one, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1")[c(1,2,4)])
meanrocvirome

# Get the variable importance
vardf <- data.frame(varImp(outmodel$finalModel))
vardf$categories <- rownames(vardf)
vardf <- vardf[order(vardf$Overall, decreasing = FALSE),]
vardf$categories <- factor(vardf$categories, levels = vardf$categories)

importanceplot <- ggplot(vardf[(length(vardf[,1])-10):(length(vardf[,1])),], aes(x=categories, y=Overall)) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  ) +
  geom_bar(stat="identity", fill=wes_palette("Royal1")[2]) +
  xlab("Categories") +
  ylab("Mean Decrease in Accuracy") +
  coord_flip()

importanceplot

#############################
# Bacteria Prediction Model #
#############################
input <- read.delim("./data/mothur16S/final.shared", header=TRUE, sep="\t")
input$Group <- sub("^(\\D)(\\d)$","\\10\\2", sub("(.)\\D+(\\d+)", "\\1\\2", input$Group, perl=TRUE))
input$Group <- as.factor(input$Group)
# Calculate as relative abundance
inputmelt <- melt(input[-c(1,3)])
# Get relative abundance
inputrelabund <- data.frame(inputmelt %>% group_by(Group) %>% mutate(RelAbund = 100 * value / sum(value)))
relabundcast <- dcast(inputrelabund, Group ~ variable, value.var = "RelAbund")
# Fix the metadata for this case without duplicate MG IDs
datadiseasesub <- datadisease[,2:3]
datadiseasesub <- datadiseasesub[!duplicated(datadiseasesub),]

# Add the disease classes
relabundclasses <- merge(relabundcast, datadiseasesub, by.x="Group", by.y="V22")

# Compare healthy to cancer
relabundclasses$V30 <- droplevels(relabundclasses$V30)

# Filter by presence/absence
pasubset <- relabundclasses[,c(colSums(relabundclasses != 0) > 30)]
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

subsetmodel <- caretmodel(pasubsetmissing)
subsetmodel

modelmeanrocbacteria <- threeclassmeanroc(subsetmodel)

# Plot the ROC curve
meanrocbacteria <- ggplot(modelmeanrocbacteria, aes(d = obs, m = one, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1")[c(1,2,4)])
meanrocbacteria

# Get the variable importance
vardfbac <- data.frame(varImp(subsetmodel$finalModel))
vardfbac$categories <- rownames(vardfbac)
vardfbac <- vardfbac[order(vardfbac$Overall, decreasing = FALSE),]
vardfbac$categories <- factor(vardfbac$categories, levels = vardfbac$categories)

importanceplotbac <- ggplot(vardfbac[(length(vardfbac[,1])-10):(length(vardfbac[,1])),], aes(x=categories, y=Overall)) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  ) +
  geom_bar(stat="identity", fill=wes_palette("Royal1")[2]) +
  xlab("Categories") +
  ylab("Mean Decrease in Accuracy") +
  coord_flip()

importanceplotbac

############################
# Quantify AUC Differences #
############################
GetAverageAUC <- function(x, y) {
	functionmodel <- caretmodel(x)
	highAUC <- functionmodel$results$Mean_ROC[order(functionmodel$results$Mean_ROC, decreasing = TRUE)[1]]
	resultdf <- data.frame(y,highAUC)
	return(resultdf)
}

iterationcount <- 10

viromeauc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
viromeaucdf <- ldply(viromeauc, data.frame)
viromeaucdf$class <- "Virus"

bacteriaauc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(pasubsetmissing, i))
bacteriaaucdf <- ldply(bacteriaauc, data.frame)
bacteriaaucdf$class <- "Bacteria"

megatron <- rbind(viromeaucdf, bacteriaaucdf)

wilcox.test(megatron$highAUC ~ megatron$class)

auccompareplot <- ggplot(megatron, aes(x = class, y = highAUC, fill = class)) +
	theme_classic() +
	theme(legend.position="none") +
	geom_boxplot(notch = FALSE) +
	scale_fill_manual(values = wes_palette("Royal1"))
auccompareplot

###############################
# Compare Bacteria and Virus  #
###############################
importanceplots <- plot_grid(importanceplot, importanceplotbac, labels = c("D", "E"), ncol = 2)
boundplot <- plot_grid(meanrocvirome, meanrocbacteria, auccompareplot, labels = c("A", "B", "C"), rel_widths = c(2, 2, 1), ncol = 3)
topbottomplot <- plot_grid(boundplot, importanceplots, rel_heights = c(3, 2), ncol = 1)

pdf("./figures/predmodel-threewayclassification.pdf", height = 5, width = 12)
	topbottomplot
dev.off()
