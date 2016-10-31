# predmodeling.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

###################
# Set Environment #
###################

write("Collapsing Contig Counts", stderr())

library("optparse")
library("plyr")
library("ggplot2")
library("caret")
library("plotROC")
library("reshape2")
library("wesanderson")

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
    # Remove excess columns
    removecols <- x[,-c(1,3)]
    # Clean up disease classes
    removecols$Group <- sub("[0-9]+", "", removecols$Group, perl=TRUE)
    if (removeclass != "none") {
        removecols <- removecols[-which(removecols$Group %in% removeclass),]
    }
    return(removecols)
}

caretmodel <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
    number = 10,
    repeats = 25,
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

################
# Run Analysis #
################
input <- read.delim(opt$input, head = FALSE, sep = "\t")

input <- read.delim("./final.shared", header=TRUE, sep="\t")

procinput <- processmothurotus(input)

# Compare healthy to cancer
inputnoademona <- processmothurotus(input, removeclass="Adenoma")
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

# outmodel <- caretmodel(inputnoademona)

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

# # Test
# set.seed(2529)
# D.ex <- rbinom(200, size = 1, prob = .5)
# M1 <- rnorm(200, mean = D.ex, sd = .65)
# M2 <- rnorm(200, mean = D.ex, sd = 1.5)

# test <- data.frame(D = D.ex, D.str = c("Healthy", "Ill")[D.ex + 1], 
#                    M1 = M1, M2 = M2, stringsAsFactors = FALSE)

# basicplot <- ggplot(test, aes(d = D, m = M1)) + geom_roc()
# basicplot
