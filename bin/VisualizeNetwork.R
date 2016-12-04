# VisualizeNetwork.R
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

graphDiameter <- function (nodeframe=nodeout, edgeframe=edgeout) {
  write("Calculating Graph Diameter", stderr())  
  # Pull out the data for clustering
  ig <- graph_from_data_frame(edgeframe, directed=F)
  connectionresult <- diameter(ig, directed=F)
  return(connectionresult)
}

connectionstrength <- function (nodeframe=nodeout, edgeframe=edgeout) {
  write("Testing Connection Strength", stderr())  
  # Pull out the data for clustering
  ig <- graph_from_data_frame(edgeframe, directed=F)
  connectionresult <- is.connected(ig, mode="strong")
  if (!connectionresult) {
    connectionresult <- is.connected(ig, mode="weak")
      if (connectionresult) {
      result <- "RESULT: Graph is weakly connected."
    } else {
      result <- "RESULT: Graph is not weakly or strongly connected."
    }
  } else {
    result <- "RESULT: Graph is strongly connected."
  }
  return(result)
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

# Test connection strength of the network
write(connectionstrength(), stderr())
write(graphDiameter(), stderr())

# Save as PDF & PNG
pdf(file="./figures/BacteriaPhageNetworkDiagram.pdf",
width=8,
height=8)
  a <- dev.cur()
  png(file="./figures/BacteriaPhageNetworkDiagram.png",
  width=8,
  height=8,
  units="in",
  res=800)
    dev.control("enable")
    plotnetwork()
    dev.copy(which=a)
  dev.off()
dev.off()

# Get number of hosts for each phage as histogram
edgecount <- count(edgeout$from)
edgecount <- edgecount[order(edgecount$freq, decreasing=FALSE),]
edgecount$x <- factor(edgecount$x, levels=edgecount$x)

edgehist <- ggplot(edgecount, aes(x=freq)) +
  theme_classic() +
  theme(axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black")) +
  geom_histogram(fill="tomato3") +
  xlab("Phage Host Count (Bacterial Strains)") +
  ylab("Frequency")

pdf(file="./figures/PhageHostHist.pdf",
width=8,
height=8)
  a <- dev.cur()
  png(file="./figures/PhageHostHist.png",
  width=8,
  height=8,
  units="in",
  res=800)
    dev.control("enable")
    edgehist
    dev.copy(which=a)
  dev.off()
dev.off()

# Get number of phages hitting bacteria
querygenus <- "
START n=node(*) MATCH (n)-[r]->(m) RETURN n.Name AS from, m.Species AS to;
"
graphoutputlist <- importgraphtodataframe(filter=0, cypherquery=querygenus)
nodeout <- as.data.frame(graphoutputlist[1])
edgeout <- as.data.frame(graphoutputlist[2])
head(nodeout)
head(edgeout)

edgecount <- count(edgeout$to)
edgecount <- edgecount[order(edgecount$freq, decreasing=FALSE),]
edgecount$x <- factor(edgecount$x, levels=edgecount$x)

edgeplotgg <- ggplot(edgecount, aes(x=x, y=freq)) +
  theme_classic() +
  theme(axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black")) +
  geom_bar(stat="identity", fill="tomato3") +
  coord_flip() +
  ylab("Unweighted Count of Targeted Phage") +
  xlab("")

pdf(file="./figures/BacteriaEdgeCount.pdf",
width=8,
height=8)
  a <- dev.cur()
  png(file="./figures/BacteriaEdgeCount.png",
  width=8,
  height=8,
  units="in",
  res=800)
    dev.control("enable")
    edgeplotgg
    dev.copy(which=a)
  dev.off()
dev.off()
