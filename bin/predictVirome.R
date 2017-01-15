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

rarefybig <- function(x, subdepth) {
	print("Creating empty matrix")
	m <- as.matrix(x)
	m[] <- 0
	m <- as.data.frame(m)

	lowbound <- 0
	
	while(lowbound < subdepth) {
		# print("Start loop")
		df <- data.frame(value = c(0, 0))
		xcount <- 0
		finalcount <- 0
		# Create reference table
		# print("Creating reference table")
		for (i in 1:length(x)) {
			cval <- x[,i]
			name <- colnames(x)[i]
			finalcount <- xcount + cval
			df[,name] <- c(xcount, finalcount)
			xcount <- finalcount
		}
		df <- as.data.frame(t(df[,-1]))
		
		# print("Evaluate random number position")
		# See where number fits in table
		numbertest <- runif(1, 0, sum(x))
		pullval <- row.names(df[c(df$V1 < numbertest & df$V2 >= numbertest),])

		# print("Save count information")
		# Save count in data frame
		m[,pullval] <- m[,pullval] + 1
		
		# print("Update reference table")
		# Update the reference table
		x[,pullval] <- x[,pullval] - 1
		
		# print("Remove zeros")
		# Remove columns whose values reached zero
		x <- x[, colSums(x) > 0]
	
		lowbound <- lowbound + 1
		write(lowbound, stderr())
	}
	return(m)
}

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

##########################
# Virus Prediction Model #
##########################
input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]
taxonomy <- read.delim("./data/mothur16S/final.taxonomy", header = TRUE, sep = "\t")
# Format taxonomy table
taxonomy$Taxonomy <- sub(".+\\;(.+)\\(\\d+\\)\\;$", "\\1", taxonomy$Taxonomy, perl=TRUE)

# Rarefy input table
minimumsubsample <- round(min(ddply(input, c("V2"), summarize, sumsamples = sum(sum))$sumsamples), 0)
minimumsubsample <- 25000
inputcast <- dcast(input, V1 ~ V2)
inputcast[is.na(inputcast)] <- 0
inputcast[,-1] <-round(inputcast[,-1],0)
row.names(inputcast) <- inputcast[,1]
inputcast <- inputcast[,-1]
inputcast <- as.data.frame(t(inputcast))
counter <- 1
rarefunction <- function(y) {
	if (sum(y) > 1e9) {
		write("Using long rarefy", stderr())
		return(rarefywell(y, minimumsubsample))
	} else {
		write("Using short rarefy", stderr())
		return(rrarefy(y, minimumsubsample))
	}
}
rareoutput <- mclapply(c(1:length(inputcast[,1])), function(i) {
	counter <<- counter + 1; 
	write(counter, stderr())
	subsetdf <- inputcast[i,]
	if(sum(subsetdf) >= minimumsubsample) {
		rareoutput <- rarefunction(subsetdf)
	}
}, mc.cores = 4)
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
outmodel

ggplot(outmodel$pred, aes(d = obs, m = Healthy)) +
	plotROC::geom_roc(n.cuts = 0, color = wes_palette("Royal1")[2]) +
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

# Iterate to get average importance plot
GetAverageImportance <- function(x, y) {
	write(y, stderr())
	functionmodel <- caretmodel(x)
	resultvardf <- data.frame(varImp(functionmodel$finalModel))
	resultvardf$categories <- rownames(resultvardf)
	resultvardf$categories <- factor(resultvardf$categories, levels = resultvardf$categories)
	resultvardf$iteration <- y
	return(resultvardf)
}

iterationcount <- 10

viromeavgimportance <- lapply(c(1:iterationcount), function(i) GetAverageImportance(absmissingid, i))
viromeavgimportancedf <- ldply(viromeavgimportance, data.frame)
virimport <- ddply(viromeavgimportancedf, c("categories"), summarize, median = median(Overall))
virimportaverage <- merge(viromeavgimportancedf, virimport, by = "categories")
virimportaverage <- virimportaverage[order(virimportaverage$median, decreasing = FALSE),]
virimportaverage$categories <- factor(virimportaverage$categories, levels = virimportaverage$categories)

importanceplot <- ggplot(virimportaverage[(length(virimportaverage[,1])-99):(length(virimportaverage[,1])),], aes(x=categories, y=Overall)) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  ) +
  geom_boxplot(fill=wes_palette("Royal1")[2]) +
  xlab("Categories") +
  ylab("Mean Decrease in Accuracy") +
  coord_flip()

importanceplot

#############################
# Bacteria Prediction Model #
#############################
inputbacteria <- read.delim("./data/mothur16S/final.shared", header=TRUE, sep="\t")
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
vardfbacmerge <- merge(vardfbac, taxonomy, by.x = "categories", by.y = "OTU")
vardfbacmerge <- vardfbacmerge[order(vardfbacmerge$Overall, decreasing = FALSE),]
vardfbacmerge <- vardfbacmerge[(length(vardfbacmerge[,1])-10):(length(vardfbacmerge[,1])),]
vardfbacmerge$categories <- factor(vardfbacmerge$categories, levels = vardfbacmerge$categories)


importanceplotbac <- ggplot(vardfbacmerge, aes(x=categories, y=Overall, fill = Taxonomy)) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  ) +
  geom_bar(stat="identity") +
  xlab("Categories") +
  ylab("Mean Decrease in Accuracy") +
  coord_flip()

importanceplotbac

#####################################
# Whole Metagenome Prediction Model #
#####################################
metainput <- read.delim("./data/BacteriaClusteredContigAbund.tsv", header=TRUE, sep="\t")

# Rarefy input table
minimumsubsample <- round(min(ddply(metainput, c("V2"), summarize, sumsamples = sum(sum))$sumsamples), 0)
metainputcast <- dcast(metainput, V1 ~ V2)
metainputcast[is.na(metainputcast)] <- 0
metainputcast[,-1] <-round(metainputcast[,-1],0)
row.names(metainputcast) <- metainputcast[,1]
metainputcast <- metainputcast[,-1]
metainputcast <- as.data.frame(t(metainputcast))
counter <- 1
rarefunction <- function(y) {
	write(counter, stderr())
	counter <<- counter + 1; 
	if (sum(y) > 1e9) {
		write("Using long rarefy", stderr())
		return(rarefywell(y, minimumsubsample))
	} else {
		write("Using short rarefy", stderr())
		return(rrarefy(y, minimumsubsample))
	}
}
metarareoutput <- lapply(c(1:length(metainputcast[,1])), function(i) {
	subsetdf <- metainputcast[i,]
	metarareoutput <- rarefunction(subsetdf)
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
metaoutmodel

ggplot(metaoutmodel$pred, aes(d = obs, m = Healthy)) +
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
vardf <- data.frame(varImp(metaoutmodel$finalModel))
vardf$categories <- rownames(vardf)
vardf <- vardf[order(vardf$Overall, decreasing = FALSE),]
vardf$categories <- factor(vardf$categories, levels = vardf$categories)

metaimportanceplot <- ggplot(vardf[(length(vardf[,1])-10):(length(vardf[,1])),], aes(x=categories, y=Overall)) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  ) +
  geom_bar(stat="identity", fill=wes_palette("Royal1")[2]) +
  xlab("Categories") +
  ylab("Mean Decrease in Accuracy") +
  coord_flip()

metaimportanceplot

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

megatron <- rbind(viromeaucdf, bacteriaaucdf, metagenomeaucdf)

pairwise.wilcox.test(x=megatron$highAUC, g=megatron$class, p.adjust.method="bonferroni")

auccompareplot <- ggplot(megatron, aes(x = class, y = highAUC, fill = class)) +
	theme_classic() +
	theme(legend.position="none") +
	geom_boxplot(notch = FALSE) +
	scale_fill_manual(values = wes_palette("Royal1"))
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
	scale_color_manual(values = wes_palette("Royal1"))
boundplot

# Compare importance
subimportance <- plot_grid(importanceplot, importanceplotcombo, labels = c("D", "E"), ncol = 2)
bottomgrid <- plot_grid(importanceplotbac, subimportance, labels = c("C"), ncol = 2)
topgrid <- plot_grid(boundplot, auccompareplot, labels = c("A", "B"), ncol = 2, rel_widths = c(6,5))
finalgridplot <- plot_grid(topgrid, bottomgrid, ncol=1, rel_heights = c(2,1))

pdf("./figures/predmodel-viromebacteria.pdf", height = 8, width = 12)
	finalgridplot
dev.off()


