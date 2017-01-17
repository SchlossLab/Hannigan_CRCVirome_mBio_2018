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

rarefybig <- function(x, subdepth) {
    print("Creating empty matrix")

    randnorep <- function(count, min, max) {
        if((max - min + 1) < count) {
            write("ERROR: Sample is out of range!", stderr())
            break
        }
        uct <- 0
        rounded <- NULL
        while(uct < count) {
            itrr <- round(runif((count - uct), min = min, max = max))
            rounded <- c(rounded, itrr)
            rounded <- unique(rounded)
            uct <- length(rounded)
        }
        rounded
    }

    m <- as.matrix(x)
    m[] <- 0
    m <- as.data.frame(m)

    df <- data.frame(value = c(0, 0))
    xcount <- 0
    finalcount <- 0
    # Create reference table
    print("Creating reference table")
    for (i in 1:length(x)) {
        cval <- x[,i]
        name <- colnames(x)[i]
        finalcount <- xcount + cval
        df[,name] <- c(xcount, finalcount)
        xcount <- finalcount
    }
    df <- as.data.frame(t(df[,-1]))
    df$names <- rownames(df)

    randdist <- randnorep(subdepth, 1, sum(x))

    print("Annotating random distribution")
    getValue <- function(x, data) {
        tmp <- data %>%
        filter(x <= V2, x > V1)
        return(tmp$names)
    }

    row <- sapply(randdist, getValue, data=df)
    row <- table(row)
    ind <- names(row)
    x[] <- 0
    x[,ind] <- row

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

rarefunction <- function(y, minsub) {
	return(rarefywell(y, minsub, average = TRUE))
}


#####################################
# Whole Metagenome Prediction Model #
#####################################
metainput <- read.delim("./data/BacteriaClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

subdepths <- c(1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7)
numitrs <- 5

subrare <- lapply(subdepths, function(a) {
	write(a, stdout())
	metainputcast <- dcast(metainput, V1 ~ V2)
	metainputcast[is.na(metainputcast)] <- 0
	metainputcast[,-1] <-round(metainputcast[,-1],0)
	row.names(metainputcast) <- metainputcast[,1]
	metainputcast <- metainputcast[,-1]
	metainputcast <- as.data.frame(t(metainputcast))

	metarareoutput <- lapply(c(1:length(metainputcast[,1])), function(i) {
		subsetdf <- metainputcast[i,]
		if(sum(subsetdf) >= a) {
			metarareoutput <- rarefunction(subsetdf, a)
		}
	})
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

	outresult <- c(a, highAUC, nrow(metaabsmissingid))

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

#####################
# Virome Prediction #
#####################
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]
taxonomy <- read.delim("./data/mothur16S/final.taxonomy", header = TRUE, sep = "\t")
# Format taxonomy table
taxonomy$Taxonomy <- sub(".+\\;(.+)\\(\\d+\\)\\;$", "\\1", taxonomy$Taxonomy, perl=TRUE)

viromesubrare <- lapply(subdepths, function(a) {
	write(a, stdout())
	# Rarefy input table
	inputcast <- dcast(input, V1 ~ V2)
	inputcast[is.na(inputcast)] <- 0
	inputcast[,-1] <-round(inputcast[,-1],0)
	row.names(inputcast) <- inputcast[,1]
	inputcast <- inputcast[,-1]
	inputcast <- as.data.frame(t(inputcast))
	counter <- 1
	rarefunction <- function(y) {
		return(rarefywell(y, a, average = TRUE))
	}
	rareoutput <- lapply(c(1:length(inputcast[,1])), function(i) {
		counter <<- counter + 1; 
		write(counter, stderr())
		subsetdf <- inputcast[i,]
		if(sum(subsetdf) >= a) {
			rareoutput <- rarefunction(subsetdf)
		}
	})
	rareoutputbind <- as.data.frame(do.call(rbind, rareoutput))
	# save(rareoutputbind, file = "./data/subsamplevirome.Rdata")
	# length(rareoutput)
	
	# load(file = "./data/subsamplevirome.Rdata")
	
	inputmelt <- melt(as.matrix(rareoutputbind))
	colnames(inputmelt) <- c("V2", "V1", "sum")
	
	inputrelabund <- data.frame(inputmelt %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
	relabundcast <- dcast(inputrelabund, V2 ~ V1, value.var = "RelAbund")
	relabundcast[is.na(relabundcast)] <- 0
	row.names(relabundcast) <- relabundcast$V2
	relabundcast <- relabundcast[,c(TRUE, colMedians(as.matrix(relabundcast[,-1])) > minavgrelabund)]
	
	castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
	castmerge <- castmerge[,-1]
	
	# Filter by presence/absence
	abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
	abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
	abssubset$V30 <- factor(abssubset$V30)
	abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
	# Get rid of the IDs
	absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]
	
	outmodel <- caretmodel(absmissingid)

	highAUC <- outmodel$results$ROC[order(outmodel$results$ROC, decreasing = TRUE)[1]]

	outresult <- c(a, highAUC, nrow(rareoutputbind))

	return(outresult)
})

viromeresult <- as.data.frame(do.call(rbind, viromesubrare))

ggplot(viromeresult, aes(x = log10(V1), y = V2)) +
	theme_classic() +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	)  +
	geom_line()

################
# Bacteria 16S #
################
inputbacteria <- read.delim("./data/mothur16S/final.shared", header=TRUE, sep="\t")

subdepths16S <- c(1e3, 5e3, 1e4, 5e4, 1e5)

bacsubrare <- lapply(subdepths16S, function(a) {
	write(a, stdout())
	counter <- 1
	rareoutput <- lapply(c(1:length(inputbacteria[,1])), function(i) {
		write(counter, stderr())
		counter <<- counter + 1; 
		subsetdf <- inputbacteria[i,-c(1:3)]
		names <- inputbacteria[i,c(1:3)]
		if(sum(subsetdf) >= a) {
			rareoutput <- rarefywell(subsetdf, a, average = TRUE)
			y <- cbind(names, rareoutput)
			y
		}
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
	
	# Compare healthy to cancer
	relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Adenoma"),]
	relabundremove$V30 <- droplevels(relabundremove$V30)
	
	# Filter by presence/absence
	pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]
	pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]
	
	subsetmodel <- caretmodel(pasubsetmissing)
	
	highAUC <- subsetmodel$results$ROC[order(subsetmodel$results$ROC, decreasing = TRUE)[1]]

	outresult <- c(a, highAUC, nrow(rareoutputbacteria))

	return(outresult)
})

bacresult <- as.data.frame(do.call(rbind, bacsubrare))

ggplot(bacresult, aes(x = log10(V1), y = V2)) +
	theme_classic() +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	)  +
	geom_line()
