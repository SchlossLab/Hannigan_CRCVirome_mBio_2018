# ConnectionImportanceCorr.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##################
# Load Libraries #
##################
packagelist <- c("RNeo4j", "ggplot2", "wesanderson", "igraph", "visNetwork", "scales", "plyr")
new.packages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
lapply(packagelist, library, character.only = TRUE)
library("ggraph")
library("caret")
library("dplyr")
library("reshape2")
library("cowplot")

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

edgecount <- function (nodeframe=nodeout, edgeframe=edgeout) {
  # Pull out the data for clustering
  ig <- graph_from_data_frame(edgeframe, directed = TRUE)
  connectionresult <- gsize(ig)
  meandistance <- mean_distance(ig, directed = FALSE)
  degreedist <- c(degree_distribution(ig, cumulative = TRUE))
  acent <- c(alpha_centrality(ig))
  pc <- c(power_centrality(ig))
  return(pc)
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
        scale_color_manual(values = wes_palette("Royal1")[c(1,2)])
  return(outputgraph)
}

##############################
# Run Analysis & Save Output #
##############################

# Start the connection to the graph
# If you are getting a lack of permission, disable local permission on Neo4J
graph <- startGraph("http://localhost:7474/db/data/", "neo4j", "neo4j")

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

alphacent <- as.data.frame(edgecount())
colnames(alphacent) <- "Centrality"
alphacent$sampleID <- rownames(alphacent)


input <- read.delim("./data/ClusteredContigAbund.tsv", header=TRUE, sep="\t")
datadisease <- read.delim("./data/metadata/MasterMeta.tsv", header=FALSE, sep="\t")[,c(2,30,22)]

inputrelabund <- data.frame(input %>% group_by(V2) %>% mutate(RelAbund = 100 * sum / sum(sum)))
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

vardf <- data.frame(varImp(outmodel$finalModel))
vardf$sampleID <- gsub("Cluster", "Phage", rownames(vardf))

mergeddf <- merge(vardf, alphacent, by = "sampleID")

cor.test(log10(mergeddf$Overall), log10(mergeddf$Centrality), method = c("pearson"))

scatterplotconnect <- ggplot(mergeddf, aes(x = log10(Overall), y = log10(Centrality))) +
  theme_classic() +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  xlab("OGU Importance (log10)") +
  ylab("OGU Alpha Centrality (log10)")

networkplot <- plotnetwork()

finalplot <- plot_grid(networkplot, scatterplotconnect, ncol = 2, labels = c("A", "B"))

pdf(file="./figures/NetworkAndScatter.pdf",
width=10,
height=5)
  finalplot
dev.off()
