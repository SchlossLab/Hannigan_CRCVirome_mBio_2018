#! /bin/bash

module load R/3.2.3

# Start neo4j server locally
/mnt/EXT/Schloss-data/bin/neo4j-enterprise-2.3.0/bin/neo4j start

Rscript ./bin/PredictRelationships.R -m "${1}" -o "${2}"

# Stop local neo4j server
/mnt/EXT/Schloss-data/bin/neo4j-enterprise-2.3.0/bin/neo4j stop
