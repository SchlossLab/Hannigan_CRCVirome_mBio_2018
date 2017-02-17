#! /bin/bash
# DownloadVirusesFromENA.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

export referencelist=$1
export outputfilename=$2

mkdir -p ./tmp-database-download

while read line; do
	wget "http://www.ebi.ac.uk/ena/data/view/${line}&display=fasta" -O ./tmp-database-download/${line}
done < ${referencelist}

cat ./tmp-database-download/* > ${outputfilename}

rm -rf ./tmp-database-download
