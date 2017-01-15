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
minsamps <- 0

rarefywell <- function(x, sample, average = FALSE)
{
    if (!identical(all.equal(x, round(x)), TRUE)) 
        stop("function is meaningful only for integers (counts)")
    x <- as.matrix(x)
    if (ncol(x) == 1)
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x))
        stop(gettextf(
             "length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length=nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    ## warn if something cannot be rarefied
    if (any(rowSums(x) < sample))
        warning("Some row sums < 'sample' and are not rarefied")
    for (i in 1:nrow(x)) {
        if (sum(x[i,]) <= sample[i]) ## nothing to rarefy: take all
            next
        if (average) {
        	x <- x / rowSums(x) * sample
        } else {
        	row <- sample(rep(nm, times=x[i,]), sample[i])
        	row <- table(row)
        	ind <- names(row)
        	x[i,] <- 0
        	x[i,ind] <- row
        }
    }
    x
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

rarefunction <- function(y, minsub) {
	return(rarefywell(y, minsub, average = TRUE))
}


#####################################
# Whole Metagenome Prediction Model #
#####################################
metainput <- read.delim("./data/BacteriaClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

subdepths <- c(1e4, 1e5, 1e6, 1e7)

subrare <- lapply(subdepths, function(a) {
	write(a, stdout())
	metainputcast <- dcast(metainput, V1 ~ V2)
	metainputcast[is.na(metainputcast)] <- 0
	metainputcast[,-1] <-round(metainputcast[,-1],0)
	row.names(metainputcast) <- metainputcast[,1]
	metainputcast <- metainputcast[,-1]
	metainputcast <- as.data.frame(t(metainputcast))

	metarareoutput <- mclapply(c(1:length(metainputcast[,1])), function(i) {
		subsetdf <- metainputcast[i,]
		if(sum(subsetdf) >= a) {
			metarareoutput <- rarefunction(subsetdf, a)
		}
	}, mc.cores = 4)
	metarareoutputbind <- as.data.frame(do.call(rbind, metarareoutput))

	metainputmelt <- melt(as.matrix(metarareoutputbind))
	colnames(metainputmelt) <- c("V2", "V1", "sum")
	
	metainputrelabund <- data.frame(metainputmelt %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
	metarelabundcast <- dcast(metainputrelabund, V2 ~ V1, value.var = "RelAbund")
	metarelabundcast[is.na(metarelabundcast)] <- 0
	row.names(metarelabundcast) <- metarelabundcast$V2
	metacastmerge <- merge(metarelabundcast, datadisease, by.x="V2", by.y="V2")
	metacastmerge <- metacastmerge[,-1]
	
	# Filter by presence/absence
	metaabssubset <- metacastmerge[!c(metacastmerge$V30 %in% "Negative"),]
	metaabssubset <- metaabssubset[!c(metaabssubset$V30 %in% "Adenoma"),]
	metaabssubset$V30 <- factor(metaabssubset$V30)
	metaabssubset <- metaabssubset[,c(colSums(metaabssubset != 0) > 0)]
	# Get rid of the IDs
	metaabsmissingid <- metaabssubset[,-which(names(metaabssubset) %in% c("V22"))]
	
	metaoutmodel <- caretmodel(metaabsmissingid)

	highAUC <- metaoutmodel$results$ROC[order(metaoutmodel$results$ROC, decreasing = TRUE)[1]]

	outresult <- c(a, highAUC, nrow(metarareoutputbind))

	return(outresult)
})
metaresult <- as.data.frame(do.call(rbind, subrare))

ggplot(metaresult, aes(x = log10(V1), y = V2)) +
	theme_classic() +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	)  +
	geom_line()

