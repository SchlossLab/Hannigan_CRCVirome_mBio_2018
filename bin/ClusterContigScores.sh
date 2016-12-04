#! /bin/bash
# ClusterContigScores.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

export ProphageBlast=$1
export BlastxResults=$2
export PfamResults=$3
export ProphageOut=$4
export BlastxOut=$5
export PfamOut=$6
export PhageContigClusters=$7
export BacteriaContigClusters=$8
export OutputName=$9
export CrisprIn=${10}
export CrisprOut=${11}


###################
# Set Subroutines #
###################

AnnotateCollapseClusters () {
	export FileToAnnotate=$1
	export OutputAnnotate=$2

	# Replace occurences
	awk -F "\t" 'FNR==NR { a[$1] = $2; next } {print a[$1]"\t"a[$2]"\t"$3}' \
		./data/${OutputName}/ContClust.tsv \
		${FileToAnnotate} \
		| tail -n +2 \
		| sed 's/\t\t/\tNoClusterPhage\t/' \
		| sed 's/^\t/NoClusterBacteria\t/' \
		> ./data/${OutputName}/tmpAnnotations.tsv

	Rscript ./bin/CollapseGeneScores.R \
		-i ./data/${OutputName}/tmpAnnotations.tsv \
		-o ${OutputAnnotate}

	# # Remove the tmp file
	# rm ./data/${OutputName}/tmpAnnotations.tsv
}

export -f AnnotateCollapseClusters

################
# Run Analysis #
################
# Make output directory
mkdir ./data/${OutputName}

# Format the contig clustering table
sed 's/,/\t/' $7 | sed 's/\t/\tPhage_/' > ./data/${OutputName}/PhageContClust.tsv
sed 's/,/\t/' $8 | sed 's/\t/\tBacteria_/' > ./data/${OutputName}/BacteriaContClust.tsv

# Merge the list
cat \
	./data/${OutputName}/PhageContClust.tsv \
	./data/${OutputName}/BacteriaContClust.tsv \
	> ./data/${OutputName}/ContClust.tsv

# Run the subroutines
AnnotateCollapseClusters \
	${CrisprIn} \
	${CrisprOut}

AnnotateCollapseClusters \
	${ProphageBlast} \
	${ProphageOut}

AnnotateCollapseClusters \
	${BlastxResults} \
	${BlastxOut}

AnnotateCollapseClusters \
	${PfamResults} \
	${PfamOut}
