#! /bin/bash
# IdentifyContigs.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

export fastafile=$1
export referencefile=$2
export idfile=$3
export outputfile=$4
export BlastPath=$5

# Get the contig representative sequences
mkdir -p ./data/tmpid

RunBlast () {
	# 1 = Query Seqs
	# 2 = ReferenceSeqs

	echo Making blast database...
	${BlastPath}makeblastdb \
		-dbtype nucl \
		-in "${2}" \
		-out ./data/tmpid/ReferenceGenomes

	echo Running tblastx...
	${BlastPath}tblastx \
    	-query "${1}" \
    	-out ${outputfile} \
    	-db ./data/tmpid/ReferenceGenomes \
    	-evalue 1e-25 \
    	-num_threads 8 \
    	-max_target_seqs 1 \
    	-outfmt 6
}

export -f RunBlast

cut -f 1 ${idfile} | tail -n +2 > ./data/tmpid/contiglist.tsv
grep -A 1 -f ./data/tmpid/contiglist.tsv ${fastafile} \
	| egrep -v "\-\-" \
	> ./data/tmpid/contigrepset.fa

RunBlast ./data/tmpid/contigrepset.fa ${referencefile}
