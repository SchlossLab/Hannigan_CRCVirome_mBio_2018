#! /bin/bash
# DownloadVirusesFromENA.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

export referencelist=$1
export outputfilename=$2

export downloadvar=$(tr '\n' ',' < ${referencelist} | sed 's/\,$//')

echo $downloadvar

wget "http://www.ebi.ac.uk/ena/data/view/$downloadvar&display=fasta" -O ${outputfilename}
