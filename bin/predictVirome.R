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
	importaverage <- importaverage[order(importaverage$mean, decreasing = FALSE),]
	importaverage$categories <- factor(importaverage$categories, levels = importaverage$categories)
	
	binlength <- c(1:topcount) + 0.5
	
	importanceplot <- ggplot(importaverage[(length(importaverage[,1])-99):(length(importaverage[,1])),], aes(x=categories, y=Overall)) +
	  theme_classic() +
	  theme(
	    axis.line.x = element_line(colour = "black"),
	    axis.line.y = element_line(colour = "black")
	  ) +
	  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", binwidth = 0.05, stackdir = "center") +
	  xlab("Categories") +
	  ylab("Mean Accuracy Decrease") +
	  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5) +
	  geom_vline(xintercept=binlength,color="grey") +
	  coord_flip()

	  return(importanceplot)
}

##########################
# Virus Prediction Model #
##########################
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]
taxonomy <- read.delim("./data/mothur16S/final.taxonomy", header = TRUE, sep = "\t")
# Format taxonomy table
taxonomy$Taxonomy <- sub(".+\\;(.+)\\(\\d+\\)\\;$", "\\1", taxonomy$Taxonomy, perl=TRUE)

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

# Filter by presence/absence
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Adenoma"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > minsamps)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

length(absmissingid[,1])

outmodel <- caretmodel(absmissingid)
outmodel

# Iterate to get average importance plot
importanceplot <- plotimportance(absmissingid)

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
pasubsetmissing <- pasubset[,-which(names(pasubset) %in% c("Group", "V2"))]

subsetmodel <- caretmodel(pasubsetmissing)
subsetmodel

importanceplotbac <- plotimportance(pasubsetmissing)

#####################################
# Whole Metagenome Prediction Model #
#####################################
metainput <- read.delim("./data/BacteriaClusteredContigAbund.tsv", header=TRUE, sep="\t")

# Rarefy input table
minimumsubsample <- 250000
metainputcast <- dcast(metainput, V1 ~ V2)
metainputcast[is.na(metainputcast)] <- 0
metainputcast[,-1] <-round(metainputcast[,-1],0)
row.names(metainputcast) <- metainputcast[,1]
metainputcast <- metainputcast[,-1]
metainputcast <- as.data.frame(t(metainputcast))
counter <- 1

metarareoutput <- lapply(c(1:length(metainputcast[,1])), function(i) {
	write(counter, stderr())
	counter <<- counter + 1; 
	subsetdf <- metainputcast[i,]
	if(sum(subsetdf) >= minimumsubsample) {
		rareoutput <- rrarefy(subsetdf, minimumsubsample)
		return(rareoutput)
	}
})
metarareoutputbind <- as.data.frame(do.call(rbind, metarareoutput))
length(metarareoutputbind[,1])

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
metaoutmodel

#################################
# Merge Bacteria and Viral Data #
#################################
virusbacteria <- merge(x = abssubset, y = pasubset, by.x = "V22", by.y = "Group")
virusbacteria <- virusbacteria[,-which(names(virusbacteria) %in% c("V30.x", "Group", "V22"))]
colnames(virusbacteria)[colnames(virusbacteria)=="V30.y"] <- "V30"

combomodel <- caretmodel(virusbacteria)
combomodel

# Get the variable importance
importanceplotcombo <- plotimportance(virusbacteria)

############################
# Quantify AUC Differences #
############################
GetAverageAUC <- function(x, y) {
	write(y, stderr())
	functionmodel <- caretmodel(x)
	highAUC <- functionmodel$results$ROC[order(functionmodel$results$ROC, decreasing = TRUE)[1]]
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

comboauc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(virusbacteria, i))
comboaucdf <- ldply(comboauc, data.frame)
comboaucdf$class <- "Combined"

metagenomeauc <- lapply(c(1:iterationcount), function(i) GetAverageAUC(metaabsmissingid, i))
metagenomeaucdf <- ldply(metagenomeauc, data.frame)
metagenomeaucdf$class <- "Metagenomic"

megatron <- rbind(viromeaucdf, bacteriaaucdf, metagenomeaucdf, comboaucdf)

pairwise.wilcox.test(x=megatron$highAUC, g=megatron$class, p.adjust.method="bonferroni")

binlength <- c(1:4) + 0.5

auccompareplot <- ggplot(megatron, aes(x = class, y = highAUC, fill = class)) +
	theme_classic() +
	theme(legend.position="none") +
	geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", binwidth = 0.005, stackdir = "center") +
	stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5) +
	geom_vline(xintercept=binlength,color="grey") +
	scale_fill_manual(values = wes_palette("Royal1")) +
	theme(
	  axis.line.x = element_line(colour = "black"),
	  axis.line.y = element_line(colour = "black")
	) +
	xlab("") +
	ylab("Area Under the Curve")

auccompareplot

###############################
# Compare Bacteria and Virus  #
###############################
subsetmodel$pred$class <- "Bacteria"
outmodel$pred$class <- "Virus"
combomodel$pred$class <- "Combined"
metaoutmodel$pred$class <- "Metagenome"
boundmodel <- rbind(subsetmodel$pred, outmodel$pred, combomodel$pred, metaoutmodel$pred)

# Plot the ROC curve
boundplot <- ggplot(boundmodel, aes(d = obs, m = Healthy, color = class)) +
	geom_roc(n.cuts = 0) +
	style_roc() +
	scale_color_manual(values = wes_palette("Royal1"), name = "Disease")
boundplot

# Compare importance
subimportance <- plot_grid(importanceplot, importanceplotcombo, labels = c("D", "E"), ncol = 2)
bottomgrid <- plot_grid(importanceplotbac, subimportance, labels = c("C"), ncol = 2)
topgrid <- plot_grid(boundplot, auccompareplot, labels = c("A", "B"), ncol = 2, rel_widths = c(6,5))
finalgridplot <- plot_grid(topgrid, bottomgrid, ncol=1, rel_heights = c(2,1))

pdf("./figures/predmodel-viromebacteria.pdf", height = 8, width = 12)
	finalgridplot
dev.off()


