# ConnectionImportanceCorr.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##################
# Load Libraries #
##################
packagelist <- c("RNeo4j", "ggplot2", "wesanderson", "igraph", "visNetwork", "scales", "plyr", "vegan")
new.packages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
lapply(packagelist, library, character.only = TRUE)
library("ggraph")
library("caret")
library("dplyr")
library("reshape2")
library("cowplot")
library("vegan")
library("extrafont")
loadfonts()
library("parallel")


###################
# Set Subroutines #
###################

importgraphtodataframe <- function (
graphconnection=graph,
cypherquery=query,
filter=0) {
  write("Retrieving Cypher Query Results", stderr())
  # Use cypher to get the edges
  edges <- cypher(graphconnection, cypherquery)
  # Filter out nodes with fewer edges than specified
  if (filter > 0) {
    # Remove the edges to singleton nodes
    singlenodes <- ddply(edges, c("to"), summarize, length=length(to))
    # # Subset because the it is not visible with all small clusters
    # singlenodesremoved <- singlenodes[c(singlenodes$length > filter),]
    multipleedge <- edges[c(which(edges$to %in% singlenodesremoved$to)),]
  } else {
    multipleedge <- edges
  }
  # Set nodes
  nodes <- data.frame(id=unique(c(multipleedge$from, multipleedge$to)))
  nodes$label <- nodes$id
  return(list(nodes, multipleedge))
}

central <- function (nodeframe=nodeout, edgeframe=edgeout) {
  # Pull out the data for clustering
  ig <- graph_from_data_frame(edgeframe, directed = FALSE)
  ec <- c(alpha_centrality(ig))
  return(ec)
}

caretmodel <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
  repeats = 5,
  number=5,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = TRUE)
  model <- train(V30~., data=x, trControl=fitControl, method="rf", metric="Mean_ROC", tuneLength=5)
  return(model)
}

tcar <- function(x) {
  fitControl <- trainControl(method = "repeatedcv",
  repeats = 5,
  number=5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE)
  model <- train(V30~., data=x, trControl=fitControl, method="rf", metric="ROC", tuneLength=5)
  return(model)
}

plotnetwork <- function (nodeframe=nodeout, edgeframe=edgeout) {
  write("Preparing Data for Plotting", stderr())
  # Pull out the data for clustering
  ig <- graph_from_data_frame(edgeframe, directed=F)
  # Set plot paramters
  V(ig)$label <- ifelse(grepl("^Bacteria", nodeframe$id),
    "Bacteria",
    "Phage")
  # Create the plot
  outputgraph <- ggraph(ig, 'igraph', algorithm = 'kk') + 
        coord_fixed() + 
        geom_edge_link0(edge_alpha = 0.05) +
        geom_node_point(aes(color = label), size = 1.5) +
        ggforce::theme_no_axes() +
        scale_color_manual(values = wes_palette("Royal1")[c(1,2)]) +
        theme_graph(border = FALSE) +
        theme(legend.position = "bottom", legend.title = element_blank())
  return(outputgraph)
}

GetAverageImportance <- function(x, y, twomodel = FALSE) {
  write(y, stderr())
  if(twomodel == TRUE) {
    print("Running two class model.")
    functionmodel <- tcar(x)
  } else {
    print("Running multiclass model.")
    functionmodel <- caretmodel(x)
  }
  resultvardf <- data.frame(varImp(functionmodel$finalModel))
  resultvardf$categories <- rownames(resultvardf)
  resultvardf$categories <- factor(resultvardf$categories, levels = resultvardf$categories)
  resultvardf$iteration <- y
  return(resultvardf)
}

##############################
# Run Analysis & Save Output #
##############################

# Start the connection to the graph
# If you are getting a lack of permission, disable local permission on Neo4J
graph <- startGraph("http://localhost:7474/db/data/", "neo4j", "root")

# Use Cypher query to get a table of the table edges
query <- "
MATCH (n)-[r]->(m)
WHERE r.Prediction = 'Interacts'
RETURN n.Name AS from, m.Species AS to;
"

graphoutputlist <- importgraphtodataframe()
nodeout <- as.data.frame(graphoutputlist[1])
edgeout <- as.data.frame(graphoutputlist[2])
head(nodeout)
head(edgeout)

alphacent <- as.data.frame(central())
colnames(alphacent) <- "Centrality"
alphacent$sampleID <- rownames(alphacent)


input <- read.delim("./data/VirusClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

# Remove the bacteria-not-virus clusters
removaltable <- read.delim("./data/contigclustersidentity/BacteriaNotVirus.tsv", header = FALSE, sep = "\t")
removaltable$V1 <- gsub("^", "Cluster_", removaltable$V1, perl = TRUE)
# Clean input
input <- input[!c(input$V1 %in% removaltable$V1),]
# Remove the nodes from the network
removenodevector <- data.frame(V1 = gsub("Cluster_", "Phage_", removaltable$V1, perl = TRUE))
edgeout <- edgeout[!c(edgeout$from %in% removenodevector$V1),]


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
  return(rrarefy(y, minimumsubsample))
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
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > 30)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

outmodel <- caretmodel(absmissingid)
outmodel

avgimportance <- mclapply(c(1:25), mc.cores = 5, function(i) GetAverageImportance(absmissingid, i))
avgimportance <- ldply(avgimportance, data.frame)
avgimportance$sampleID <- gsub("Cluster", "Phage", avgimportance$categories)
davg <- ddply(avgimportance, c("sampleID"), summarize, mean = mean(Overall))
mergeddf <- merge(alphacent, davg, by = "sampleID")

scor <- cor.test(mergeddf$Centrality, mergeddf$mean, method = "spearman")

scatterplotconnect <- ggplot(mergeddf, aes(x = Centrality, y = log10(mean))) +
  theme_classic() +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  ylab("OVU Importance (log10)") +
  xlab("OVU Alpha Centrality")

scatterplotconnect

networkplot <- plotnetwork()

finalplot <- plot_grid(networkplot, scatterplotconnect, ncol = 2, labels = c("A", "B"))
finalplot

pdf(file="./figures/NetworkAndScatter.pdf",
width=10,
height=5)
  finalplot
dev.off()

scorw <- data.frame("names" = c("pval", "rho"), "scores" = c(scor$p.value, as.numeric(scor$estimate)))

write.table(scorw, file = "./rtables/cor-stats.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##############################################
# Compare Importance Centrality Across Times #
##############################################
### Healthy Vs Adenoma
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Cancer"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > 30)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

caretmodel(absmissingid)

avgimportance <- lapply(c(1:25), function(i) GetAverageImportance(absmissingid, i, twomodel = TRUE))
avgimportance <- ldply(avgimportance, data.frame)
avgimportance$sampleID <- gsub("Cluster", "Phage", avgimportance$categories)
davg <- ddply(avgimportance, c("sampleID"), summarize, mean = mean(Overall))
hva <- merge(alphacent, davg, by = "sampleID")

### Adenoma Vs Cancer
abssubset <- castmerge[!c(castmerge$V30 %in% "Negative"),]
abssubset <- abssubset[!c(abssubset$V30 %in% "Healthy"),]
abssubset$V30 <- factor(abssubset$V30)
abssubset <- abssubset[,c(colSums(abssubset != 0) > 30)]
# Get rid of the IDs
absmissingid <- abssubset[,-which(names(abssubset) %in% c("V22"))]

avgimportance <- lapply(c(1:25), function(i) GetAverageImportance(absmissingid, i, twomodel = TRUE))
avgimportance <- ldply(avgimportance, data.frame)
avgimportance$sampleID <- gsub("Cluster", "Phage", avgimportance$categories)
davg <- ddply(avgimportance, c("sampleID"), summarize, mean = mean(Overall))
avc <- merge(alphacent, davg, by = "sampleID")


hvao <- hva[c(order(hva$mean, decreasing = TRUE)),][1:10,]
hvao$class <- "hva"
avco <- avc[c(order(avc$mean, decreasing = TRUE)),][1:10,]
avco$class <- "avc"

bo <- rbind(hvao, avco)

ggplot(bo, aes(x = class, y = Centrality, group = class)) +
  theme_classic() +
  geom_dotplot(fill=wes_palette("Royal1")[2], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 0.5)

wilcox.test(x = bo$Centrality, g = bo$class)
