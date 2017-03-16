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
iterationcount <- 10


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

plotimportance <- function(x, iterationcount = 10, topcount = 10) {
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
	  # geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 0.5) +
	  xlab("Virus Identity") +
	  ylab("Mean Accuracy Decrease") +
	  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.75) +
	  stat_summary(fun.y = mean, geom = "point", width = 0.75) +
	  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
	  geom_vline(xintercept=binlength,color="grey") +
	  scale_x_discrete(
	  	labels=parse(
	  		text = paste0(
	  			"italic('",
	  			dfplot[c(0:(topcount - 1))*iterationcount+1,"V7"],
	  			"')~",
	  			paste(
	  				" (",
	  				dfplot[c(0:(topcount - 1))*iterationcount+1,"categories"],
	  				")",
	  				sep = ""
	  			)
	  		)
	  	)
	  ) +
	  coord_flip()

	  return(importanceplot)
}

pbi <- function(x, iterationcount = 10, topcount = 10) {
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
	  # geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 0.5) +
	  xlab("Bacteria Identity") +
	  ylab("Mean Accuracy Decrease") +
	  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.75) +
	  stat_summary(fun.y = mean, geom = "point", width = 0.75) +
	  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
	  geom_vline(xintercept=binlength,color="grey") +
	  scale_x_discrete(
	  	labels=parse(
	  		text = paste0(
	  			"italic('",
	  			dfplot[c(0:(topcount - 1))*iterationcount+1,"Taxonomy"],
	  			"')~",
	  			paste(
	  				" (",
	  				dfplot[c(0:(topcount - 1))*iterationcount+1,"categories"],
	  				")",
	  				sep = ""
	  			)
	  		)
	  	)
	  ) +
	  coord_flip()

	  return(importanceplot)
}

GetAverageAUC <- function(x, y) {
	write(y, stderr())
	functionmodel <- caretmodel(x)
	highAUC <- functionmodel$results$ROC[order(functionmodel$results$ROC, decreasing = TRUE)[1]]
	highSpec <- functionmodel$results$Spec[order(functionmodel$results$Spec, decreasing = TRUE)[1]]
	highSens <- functionmodel$results$Sens[order(functionmodel$results$Sens, decreasing = TRUE)[1]]
	resultdf <- data.frame(y,highAUC, highSpec, highSens)
	return(resultdf)
}

##########################
# Virus Prediction Model #
##########################
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]
taxonomy <- read.delim("./data/mothur16S/final.taxonomy", header = TRUE, sep = "\t")
virustax <- read.delim("./data/contigclustersidentity/clustax.tsv", header = FALSE, sep = "\t")
# Format taxonomy table
taxonomy$Taxonomy <- sub(".+\\;(.+)\\(\\d+\\)\\;$", "\\1", taxonomy$Taxonomy, perl=TRUE)
virustax$V1 <- sub("^", "Cluster_", virustax$V1, perl = TRUE)
virustax <- unique(virustax)

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

castmerge <- merge(relabundcast, datadisease, by.x="V2", by.y="V2")
castmerge <- castmerge[,-1]

### Healthy VS Adenoma
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Cancer"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

length(absmissingid[,1])

hva <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
hvadf <- ldply(hva, data.frame)
hvadf$class <- "hva"

# Iterate to get average importance plot
hvai <- plotimportance(absmissingid)

### Adenoma VS Cancer
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Healthy"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

length(absmissingid[,1])

avc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
avcdf <- ldply(avc, data.frame)
avcdf$class <- "avc"

avci <- plotimportance(absmissingid)

### Healthy VS Cancer
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

length(absmissingid[,1])

hvc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
hvcdf <- ldply(hvc, data.frame)
hvcdf$class <- "hvc"

hvci <- plotimportance(absmissingid)

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
# If you change to less than instead of less than OR equal, you lose most performance
relabundcast <- relabundcast[,c(TRUE, colMedians(as.matrix(relabundcast[,-1])) >= minavgrelabund)]

# Fix the metadata for this case without duplicate MG IDs
datadiseasesub <- datadisease[,2:3]
datadiseasesub <- datadiseasesub[!duplicated(datadiseasesub),]

# Add the disease classes
relabundclasses <- merge(relabundcast, datadiseasesub, by.x="Group", by.y="V22")

### Healthy vs Adenoma
# Filter
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Cancer"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

bhva <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
bhvadf <- ldply(bhva, data.frame)
bhvadf$class <- "bhva"

bhvai <- pbi(pasubsetmissing)

### Adenoma vs Cancer
# Filter
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Healthy"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

bavc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
bavcdf <- ldply(bavc, data.frame)
bavcdf$class <- "bavc"

bavci <- pbi(pasubsetmissing)

### Healthy vs Cancer
# Filter
relabundremove <- relabundclasses[!c(relabundclasses$V30 %in% "Adenoma"),]
relabundremove$V30 <- droplevels(relabundremove$V30)

# Filter by presence/absence
pasubset <- relabundremove[,c(colSums(relabundremove != 0) > minsamps)]
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

bhvc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(absmissingid, i))
bhvcdf <- ldply(bhvc, data.frame)
bhvcdf$class <- "bhvc"

bhvci <- pbi(pasubsetmissing)


###############################
# Compare Bacteria and Virus  #
###############################
vauc <- rbind(hvadf, avcdf, hvcdf)
vauc$class <- factor(vauc$class, levels = vauc$class)

bauc <- rbind(bhvadf, bavcdf, bhvcdf)
bauc$class <- factor(bauc$class, levels = bauc$class)

vplot <- ggplot(vauc, aes(x = class, y = highAUC)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle=45, hjust = 1)) +
	geom_point() +
	ylim(0, 1) +
	stat_summary(fun.y=mean, colour="blue", geom="line", aes(group = 1)) +
	ylab("Virus Model AUC") +
	xlab("") +
	scale_x_discrete(labels = c("Healthy vs Adenoma", "Adenoma vs Cancer", "Healthy vs Cancer"))

bplot <- ggplot(bauc, aes(x = class, y = highAUC)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle=45, hjust = 1)) +
	geom_point() +
	ylim(0, 1) +
	stat_summary(fun.y=mean, colour="blue", geom="line", aes(group = 1)) +
	ylab("Bacteria Model AUC") +
	xlab("") +
	scale_x_discrete(labels = c("Healthy vs Adenoma", "Adenoma vs Cancer", "Healthy vs Cancer"))

top <- plot_grid(vplot, bplot, labels = c("A", "B"))
vcol <- plot_grid(
	hvai + ggtitle("Healthy vs Adenoma") + theme(plot.title = element_text(hjust = 0)),
	avci + ggtitle("Adenoma vs Cancer") + theme(plot.title = element_text(hjust = 0)),
	hvci + ggtitle("Healthy vs Cancer") + theme(plot.title = element_text(hjust = 0)),
	ncol = 1,
	labels = c("C", "D", "E"))
bcol <- plot_grid(
	bhvai + ggtitle("Healthy vs Adenoma") + theme(plot.title = element_text(hjust = 0)),
	bavci + ggtitle("Adenoma vs Cancer") + theme(plot.title = element_text(hjust = 0)),
	bhvci + ggtitle("Healthy vs Cancer") + theme(plot.title = element_text(hjust = 0)),
	ncol = 1,
	labels = c("F", "G", "H"))
vb <- plot_grid(vcol, bcol, nrow = 1)
finalp <- plot_grid(top, vb, ncol = 1, rel_heights = c(1,1.75))
finalp

pdf("./figures/temporalmodels.pdf", height = 12, width = 10)
	finalp
dev.off()

write.table(aucstats, file = "./rtables/twoclass-auc.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
colnames(aucpv) <- c("First", "Second", "pval")
write.table(aucpv, file = "./rtables/twoclass-auc-statsig.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

aucpv[c(aucpv$First %in% "Virus" & aucpv$Second %in% "Bacteria"),"pval"]
