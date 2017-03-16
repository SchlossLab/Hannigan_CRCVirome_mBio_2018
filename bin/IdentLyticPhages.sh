#! /bin/bash
# IdentLyticPhages.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Identify contigs as being lysogenic if they are similar to integrase genes, bacterial
# genomes, or ACLAME elements.

export fastafile=$1
export idfile=$2
export BlastPath=$3
export integrase="./data/metadata/phage-integrase.fa"
export aclame="./data/metadata/aclame_proteins_prophages_0.4.fa"
export bacteriadb="./data/metadata/BacteriaReference.fa"

mkdir -p ./data/tmpidlytic
mkdir -p ./data/repcycle

RunBlastx () {
	# 1 = Query Seqs
	# 2 = ReferenceSeqs

	echo Making blast database...
	${BlastPath}makeblastdb \
		-dbtype prot \
		-in "${2}" \
		-out ./data/tmpidlytic/ReferenceGenomes

	echo Running blastx...
	${BlastPath}blastx \
    	-query "${1}" \
    	-out "${3}" \
    	-db ./data/tmpidlytic/ReferenceGenomes \
    	-evalue 1e-5 \
    	-num_threads 8 \
    	-max_target_seqs 1 \
    	-outfmt 6
}

RunBlastn () {
	# 1 = Query Seqs
	# 2 = ReferenceSeqs

	echo Making blast database...
	${BlastPath}makeblastdb \
		-dbtype nucl \
		-in "${2}" \
		-out ./data/tmpidlytic/ReferenceGenomes

	echo Running blastn...
	${BlastPath}blastn \
    	-query "${1}" \
    	-out "${3}" \
    	-db ./data/tmpidlytic/ReferenceGenomes \
    	-evalue 1e-25 \
    	-num_threads 8 \
    	-max_target_seqs 1 \
    	-outfmt 6
}

export -f RunBlastx
export -f RunBlastn

cut -f 1 ${idfile} | tail -n +2 > ./data/tmpidlytic/contiglist.tsv
grep -A 1 -f ./data/tmpidlytic/contiglist.tsv ${fastafile} \
	| egrep -v "\-\-" \
	> ./data/tmpidlytic/contigrepset.fa

# Align contigs to the integrase dataset
RunBlastx \
	./data/tmpidlytic/contigrepset.fa \
	${integrase} \
	./data/repcycle/intblastout.tsv

# Align contigs to ACLAME database
RunBlastx \
	./data/tmpidlytic/contigrepset.fa \
	${aclame} \
	./data/repcycle/aclameblastout.tsv

# Align contigs to bacterial database
RunBlastn \
	./data/tmpidlytic/contigrepset.fa \
	${bacteriadb} \
	./data/repcycle/phagebacteriablastout.tsv

# Combine the blast hits
cat \
	./data/repcycle/intblastout.tsv \
	./data/repcycle/aclameblastout.tsv \
	./data/repcycle/phagebacteriablastout.tsv \
	| cut -f 1 \
	| sort \
	| uniq \
	> ./data/repcycle/phagereplicationcycle.tsv

# Add the contig IDs to the table
awk -F "\t" 'FNR==NR { a[$1] = $3; next } { for( i in a ) if($1 ~ i) {print "Cluster_"a[$1]"\tLysogenic"} next }' \
	./data/contigclustersidentity/longestcontigsvirus.tsv \
	./data/repcycle/phagereplicationcycle.tsv \
	> ./data/repcycle/repclust.tsv



