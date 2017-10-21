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
iterationcount <- 15


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
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 0.5) +
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
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 0.5) +
	  xlab("Bacteria Identity") +
	  ylab("Mean Accuracy Decrease") +
	  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.75) +
	  geom_vline(xintercept=binlength,color="grey") +
	  scale_x_discrete(labels=dfplot[c(0:(topcount - 1))*iterationcount+1,"Taxonomy"]) +
	  coord_flip()

	  return(importanceplot)
}

pcomb <- function(x, iterationcount = 25, topcount = 10) {
	avgimportance <- lapply(c(1:iterationcount), function(i) GetAverageImportance(x, i))
	avgimportancedf <- ldply(avgimportance, data.frame)
	import <- ddply(avgimportancedf, c("categories"), summarize, mean = mean(Overall))
	importaverage <- merge(avgimportancedf, import, by = "categories")
	coltax <- virustax[,c(1, 3, 7)]
	colnames(coltax) <- c("OTU", "Size", "Taxonomy")
	mtax <- rbind(taxonomy, coltax)
	importaverage <- merge(importaverage, mtax, by.x = "categories", by.y = "OTU", all = TRUE)
	importaverage$Taxonomy <- as.character(importaverage$Taxonomy)
	importaverage <- importaverage[!c(importaverage$Overall %in% NA),]
	importaverage[is.na(importaverage)] <- "Unknown"
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
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 0.5) +
	  xlab("Microbe Identity") +
	  ylab("Mean Accuracy Decrease") +
	  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.75) +
	  geom_vline(xintercept=binlength,color="grey") +
	  scale_x_discrete(labels=dfplot[c(0:(topcount - 1))*iterationcount+1,"Taxonomy"]) +
	  coord_flip()

	  return(importanceplot)
}

nestedcv <- function(x, iterations = 5, split = 0.75) {
  outlist <- lapply(1:iterations, function(i) {
    write(i, stdout())

    trainIndex <- createDataPartition(
      x$V30,
      p = split,
      list = FALSE, 
      times = 1
    )

    dftrain <- x[trainIndex,]
    dftest <- x[-trainIndex,]
    
    outmodel <- caretmodel(dftrain)
    
    x_test <- dftest[,-length(dftest)]
    y_test <- dftest[,length(dftest)]
    
    outpred <- predict(outmodel, x_test, type="prob")
    outpred$pred <- predict(outmodel, x_test)
    outpred$obs <- y_test

    # confusionMatrix(outpred, y_test)
    # postResample(pred = outpred$pred, obs = outpred$obs)
    sumout <- twoClassSummary(outpred, lev = levels(outpred$obs))
    sumroc <- sumout[[1]]
    sumsens <- sumout[[2]]
    sumspec <- sumout[[3]]
    return(list(sumroc, outpred, sumsens, sumspec))
  })
  
  # Get the max and min values
  rocpositions <- sapply(outlist,`[`,1)
  maxl <- outlist[[match(max(unlist(rocpositions)), rocpositions)]]
  medl <- outlist[[match(median(unlist(rocpositions)), rocpositions)]]
  minl <- outlist[[match(min(unlist(rocpositions)), rocpositions)]]

  return(c(maxl, medl, minl))
}

GetAverageAUC <- function(x, y) {
	write(y, stderr())
	functionmodel <- nestedcv(x)
	avgAUC <- functionmodel[[5]]
	avgSpec <- functionmodel[[8]]
	avgSens <- functionmodel[[7]]
	resultdf <- data.frame(y, avgAUC, avgSpec, avgSens)
	return(resultdf)
}

######################
# Start the Analysis #
######################
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]
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

###################################
# Bacteria Without Zero Filtering #
###################################

inputbacteria$Group <- sub("^(\\D)(\\d)$","\\10\\2", sub("(.)\\D+(\\d+)", "\\1\\2", inputbacteria$Group, perl=TRUE))
inputbacteria$Group <- as.factor(inputbacteria$Group)
# Calculate as relative abundance
inputbacteriamelt <- melt(inputbacteria[-c(1,3)])
# Get relative abundance
inputbacteriarelabund <- data.frame(inputbacteriamelt %>% group_by(Group) %>% mutate(RelAbund = 100 * value / sum(value)))
relabundcast <- dcast(inputbacteriarelabund, Group ~ variable, value.var = "RelAbund")
# If you change to less than instead of less than OR equal, you lose most performance
relabundcast <- relabundcast[,c(TRUE, colMedians(as.matrix(relabundcast[,-1])) >= minavgrelabund)]

# Fix the metadata for this case without duplicate MG IDs
datadiseasesub <- datadisease[,2:3]
datadiseasesub <- datadiseasesub[!duplicated(datadiseasesub),]

# Add the disease classes
relabundclasses <- merge(relabundcast, datadiseasesub, by.x="Group", by.y="V22")

# Compare healthy to cancer
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Adenoma"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]
pasubsetmissingwozf <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

caretmodel(pasubsetmissingwozf)
nestedsubset <- nestedcv(pasubsetmissingwozf)

withoutzerofilter <- lapply(c(1:iterationcount), function(i) GetAverageAUC(pasubsetmissingwozf, i))
withoutzerofilterdf <- ldply(withoutzerofilter, data.frame)
withoutzerofilterdf$class <- "WithoutZeroFilter"

###################################
# Bacteria Without Zero Filtering #
###################################
inputbacteria$Group <- sub("^(\\D)(\\d)$","\\10\\2", sub("(.)\\D+(\\d+)", "\\1\\2", inputbacteria$Group, perl=TRUE))
inputbacteria$Group <- as.factor(inputbacteria$Group)
# Calculate as relative abundance
inputbacteriamelt <- melt(inputbacteria[-c(1,3)])
# Get relative abundance
inputbacteriarelabund <- data.frame(inputbacteriamelt %>% group_by(Group) %>% mutate(RelAbund = 100 * value / sum(value)))
relabundcast <- dcast(inputbacteriarelabund, Group ~ variable, value.var = "RelAbund")
# If you change to less than instead of less than OR equal, you lose most performance
relabundcast <- relabundcast[,c(TRUE, colMedians(as.matrix(relabundcast[,-1])) > minavgrelabund)]

# Fix the metadata for this case without duplicate MG IDs
datadiseasesub <- datadisease[,2:3]
datadiseasesub <- datadiseasesub[!duplicated(datadiseasesub),]

# Add the disease classes
relabundclasses <- merge(relabundcast, datadiseasesub, by.x="Group", by.y="V22")

# Compare healthy to cancer
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Adenoma"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]
pasubsetmissingwzf <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

caretmodel(pasubsetmissingwzf)

withzerofilter <- lapply(c(1:iterationcount), function(i) GetAverageAUC(pasubsetmissingwzf, i))
withzerofilterdf <- ldply(withzerofilter, data.frame)
withzerofilterdf$class <- "WithZeroFilter"

#########################
# ID Those Filtered Out #
#########################
# ID those that were filtered out
FilteredID <- setdiff(colnames(pasubsetmissingwozf), colnames(pasubsetmissingwzf))
FilteredID

filtra <- pasubsetmissingwozf[,c(colnames(pasubsetmissingwozf) %in% FilteredID | colnames(pasubsetmissingwozf) %in% "V30")]
filtram <- melt(filtra)

binlength <- c(1:6) + 0.5

frabund <- ggplot(filtram, aes(x = factor(variable), y = (value + 1e-06), fill = factor(V30))) +
	theme_classic() +
	geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.75, stackratio = 0.35) +
	stat_summary_bin(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5, position = "dodge") +
	scale_y_log10() +
	geom_vline(xintercept=binlength,color="grey") +
	  theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	xlab("Operational Taxonomic Units") +
	ylab("Relative Abundance") +
	scale_fill_manual(values = wes_palette("Royal1")[c(1,2,4)], name = "Disease") +
	geom_hline(yintercept = 1, linetype = "dashed")

fmod <- caretmodel(filtra)
fmod

##################
# Combo Analysis #
##################
# Combine the data frames with and without filtering
cf <- rbind(withoutzerofilterdf, withzerofilterdf)

statsig <- wilcox.test(cf$avgAUC ~ cf$class)

compplot <- ggplot(cf, aes(x = class, y = avgAUC, fill = class)) +
	theme_classic() +
	theme(legend.position="none") +
	geom_dotplot(fill=wes_palette("Royal1")[4], binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.6) +
	stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5) +
	scale_fill_manual(values = wes_palette("Royal1")) +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	xlab("") +
	ylab("Area Under the Curve") +
	geom_segment(x = 1, xend = 2, y = 0.90, yend = 0.90) +
	annotate("text", x = 1.5, y = 0.92, label = paste("p-value = ", signif(statsig$p.value, digits = 2), sep = ""), size = 4) +
	ylim(0, 1) +
	geom_hline(yintercept = 0.5, linetype = "dashed")

allplot <- plot_grid(compplot, frabund, rel_widths = c(1.5,5), labels = c("A", "B"))

pdf("./figures/filtered16S.pdf", height = 6, width = 16)
	allplot
dev.off()

