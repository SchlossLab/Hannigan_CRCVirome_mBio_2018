library("optparse")
library("ggplot2")
library("plotROC")
library("reshape2")
library("wesanderson")
library("vegan")
library("plyr")
library("dplyr")
library("caret")

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

lowabundgone <- castmerge[,c(sapply(castmerge[-length(castmerge)], median) > 0.01,TRUE)]
lowsubset <- lowabundgone[!c(lowabundgone$V30 %in% "Negative"),]
lowsubset <- lowsubset[!c(lowsubset$V30 %in% "Adenoma"),]
lowsubset$V30 <- factor(lowsubset$V30)

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

outmodel <- caretmodel(lowsubset)
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

ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE)
trainedlogit <- train(V30~., data = lowsubset, method="glm", family="binomial", trControl = ctrl)

ggplot(trainedlogit$pred, aes(d = obs, m = "Healthy")) +
    geom_roc(n.cuts = 0, color = wes_palette("Royal1")[2]) +
    theme_classic() +
    theme(
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    ) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype=2, colour=wes_palette("Royal1")[1]) +
    ylab("Sensitivity") +
    xlab(paste("Inverse Specificity"))