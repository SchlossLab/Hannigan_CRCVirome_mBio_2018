#! /usr/local/bin/R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

write("Predicting Phage-Bacteria Interactions", stderr())

##################################
# Install Dependencies if Needed #
##################################
packagelist <- c("RNeo4j", "ggplot2", "optparse", "caret", "wesanderson", "plotROC")
new.packages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packagelist, library, character.only = TRUE)

option_list <- list(
  make_option(c("-m", "--input"),
    type = "character",
    default = NULL,
    help = "Contig count table formatted from bowtie2.",
    metavar = "character"),
  make_option(c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Table of predicted interactions.",
    metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

###################
# Set Subroutines #
###################
getresults <- function(x) {
  x[,3:6] <- as.data.frame(sapply(x[,3:6], as.numeric))
  rownames(x) <- NULL
  x[is.na(x)] <- 0
  return(x)
}

################
# Run Analysis #
################
# Load in R data file that contains predictive model
load(opt$input)

# Start the connection to the graph
# If you are getting a lack of permission, disable local permission on Neo4J
graph <- startGraph("http://localhost:7474/db/data/", "neo4j", "neo4j")

queryresults <- "
MATCH (n)-[r]->(m)
RETURN
m.Name as Bacteria,
n.Name as Phage,
r.CRISPR as CRISPR,
r.BLAST as Blast,
r.BLASTX as Blastx,
r.PFAM as Pfam;
"

# Run the cypher queries
querydata <- cypher(graph, queryresults)

datadef <- getresults(querydata)

comdf <- data.frame(datadef[apply(datadef[, -c(1:2)], MARGIN = 1, function(x) any(x > 0)),])

predoutput <- predict(outmodel, newdata=comdf)

summary(predoutput)

predresulttable <- cbind(comdf[,c(1:2)], predoutput)
colnames(predresulttable) <- c("Bacteria", "Phage", "InteractionScore")

write.table(
  x = predresulttable,
  file = opt$out,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write("Completed Predicting Interactions", stderr())
