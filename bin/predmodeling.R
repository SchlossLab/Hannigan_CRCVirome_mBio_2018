# predmodeling.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

###################
# Set Environment #
###################

write("Collapsing Contig Counts", stderr())
set.seed(1234)

library("optparse")
library("plyr")
library("dplyr")
library("ggplot2")
# library("caret")
# library("plotROC")
library("reshape2")
library("wesanderson")
library("AUCRF")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "OTU abundance table straight from Mothur.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw results.",
    metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

###################
# Set Subroutines #
###################
processmothurotus <- function(x, removeclass="none") {
    # Clean up disease classes
    x$Group <- sub("[0-9]+", "", x$Group, perl=TRUE)
    if (removeclass != "none") {
        x <- x[-which(x$Group %in% removeclass),]
    }
    return(x)
}

caretmodel <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
    number = 5,
    repeats = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = TRUE)
  model <- train(Group~., data=x, trControl=fitControl, method="rf", metric="ROC")
  return(model)
}

plotscree <- function(x) {
    # Calculate the percent variance accounted for by each component
    pcaPercentVar <- x$sd^2/sum(x$sd^2)*100
    PcaScree <- melt(pcaPercentVar)
    PcaScree$name <- sequence(length(row.names(PcaScree)))
    ggplot(PcaScree, aes(x=name, y=value)) +
        theme_classic() +
        geom_point() +
        geom_path() +
        xlab("PCA Component") +
        ylab("Percent Variance Explained")
}

#######
# PCA #
#######
input <- read.delim("./final.shared", header=TRUE, sep="\t")
# Calculate as relative abundance
inputmelt <- melt(input[-c(1,3)])
# Get relative abundance
inputrelabund <- data.frame(inputmelt %>% group_by(Group) %>% mutate(RelAbund = 100 * value / sum(value)))
relabundcast <- dcast(inputrelabund, Group ~ variable, value.var = "RelAbund")

# Compare healthy to cancer
inputnoademona <- processmothurotus(relabundcast, removeclass="Adenoma")
# pca to reduce dimentions
inputtrans <- t(inputnoademona)
colnames(inputtrans) <- inputtrans[1,]
inputtrans <- data.frame(inputtrans[-1,], stringsAsFactors=FALSE)
inputtrans <- sapply( inputtrans, as.numeric )
pcainput <- prcomp(inputtrans)
screeplot <- plotscree(pcainput)
# Get the categories as the first column
pcadf <- data.frame(pcainput$rotation)
pcadf$Group <- sub(".[0-9]+", "", rownames(pcadf), perl=TRUE)
outmodel <- caretmodel(pcadf)
plot(outmodel)

# Plot the ROC curve
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



#######################
# Abundance Reduction #
#######################
inputnoademona <- processmothurotus(relabundcast, removeclass="Adenoma")
lowabundgone = inputnoademona[,c(TRUE, sapply(inputnoademona[-1], median) > 0.01)]

correlationMatrix <- cor(lowabundgone[,-1])
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75) + 1
lowabundcor <- lowabundgone[,-c(highlyCorrelated)]
outmodellowabund <- caretmodel(lowabundcor)
plot(outmodellowabund)

ggplot(outmodellowabund$pred, aes(d = obs, m = Healthy)) +
    geom_roc(n.cuts = 0, color = wes_palette("Royal1")[2]) +
    theme_classic() +
    theme(
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    ) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype=2, colour=wes_palette("Royal1")[1]) +
    ylab("Sensitivity") +
    xlab(paste("Inverse Specificity"))



#################################
# Recursive Feature Elimination #
#################################
inputnoademona <- processmothurotus(relabundcast, removeclass="Adenoma")
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results <- rfe(lowabundgone[,-1], lowabundgone[,1], rfeControl=control)
plot(results, type=c("g", "o"))

inputnoademona <- processmothurotus(relabundcast, removeclass="Adenoma")
lowabundgone = inputnoademona[,c(TRUE, sapply(inputnoademona[-1], median) > 0.1)]
lowabundgone$Group <- factor(lowabundgone$Group)
fit_les_model <- AUCRF(Group~., data=lowabundgone, pdel=0.05, ranking='MDA')



