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
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

inputrelabund <- data.frame(input %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
relabundcast <- dcast(inputrelabund, V2 ~ V1, value.var = "RelAbund")
relabundcast[is.na(relabundcast)] <- 0
row.names(relabundcast) <- relabundcast$V2
castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
castmerge <- castmerge[,-1]

# Filter by presence/absence
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > 30)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

outmodel <- caretmodel(absmissingid)
outmodel

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
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Adenoma"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > 30)]
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

subsetmodel <- caretmodel(pasubsetmissing)
subsetmodel

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

#################################
# Merge Bacteria and Viral Data #
#################################
virusbacteria <- merge(x = abssubset, y = pasubset, by.x = "V22", by.y = "Group")
virusbacteria <- virusbacteria[,-which(names(virusbacteria) %in% c("V30.x", "Group", "V22"))]
colnames(virusbacteria)[colnames(virusbacteria)=="V30.y"] <- "V30"

combomodel <- caretmodel(virusbacteria)
combomodel

# Plot the ROC curve
ggplot(combomodel$pred, aes(d = obs, m = Healthy)) +
	geom_roc(n.cuts = 0, color = wes_palette("Royal1")[2]) +
	theme_classic() +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype=2, colour=wes_palette("Royal1")[1]) +
	ylab("Sensitivity") +
	xlab(paste("Inverse Specificity"))

# Get the variable importance
combovar <- data.frame(varImp(combomodel$finalModel))
combovar$categories <- rownames(combovar)
combovar <- combovar[order(combovar$Overall, decreasing = FALSE),]
combovar$categories <- factor(combovar$categories, levels = combovar$categories)

importanceplotcombo <- ggplot(combovar[(length(combovar[,1])-10):(length(combovar[,1])),], aes(x=categories, y=Overall)) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  ) +
  geom_bar(stat="identity", fill=wes_palette("Royal1")[2]) +
  xlab("Categories") +
  ylab("Mean Decrease in Accuracy") +
  coord_flip()

importanceplotcombo

###############################
# Compare Bacteria and Virus  #
###############################
subsetmodel$pred$class <- "Bacteria"
outmodel$pred$class <- "Virus"
combomodel$pred$class <- "Combined"
boundmodel <- rbind(subsetmodel$pred, outmodel$pred, combomodel$pred)

# Plot the ROC curve
boundplot <- ggplot(boundmodel, aes(d = obs, m = Healthy, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1"))
boundplot

# Compare importance
vardf$class <- "Virus"
vardfbac$class <- "Bacteria"
varimpbind <- rbind(vardf, vardfbac)
ggplot(varimpbind, aes(y = Overall, group = ))
plot_grid(importanceplot, importanceplotbac, labels = c("A", "B"))

pdf("./figures/predmodel-viromebacteria.pdf", height = 5, width = 6)
	boundplot
dev.off()
