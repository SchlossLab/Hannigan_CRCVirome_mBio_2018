#! /bin/bash
# ClusterContigAbund.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

export ContigAbundance=$1
export ContigClusters=$2
export TransformedAbundanceFile=$3
export OutputName="tmpContigAbund"


###################
# Set Subroutines #
###################

AnnotateCollapseClusters () {
	FileToAnnotate=$1
	OutputAnnotate=$2

	# Replace occurences
	awk -F "\t" 'FNR==NR { a[$1] = $2; next } {print "Cluster_"a[$1]"\t"$3"\t"$2}' \
		./data/${OutputName}/ContClust.tsv \
		${FileToAnnotate} \
		> ./data/${OutputName}/tmpAnnotations.tsv

	Rscript ./bin/CollapseGeneScores.R \
		-i ./data/${OutputName}/tmpAnnotations.tsv \
		-o ./data/${OutputName}/tmpAnnotations2.tsv

	grep -v 'NA' ./data/${OutputName}/tmpAnnotations2.tsv > ${OutputAnnotate}

	# Remove the tmp file
	rm ./data/${OutputName}/tmpAnnotations.tsv
	rm ./data/${OutputName}/tmpAnnotations2.tsv
}

export -f AnnotateCollapseClusters

################
# Run Analysis #
################
# Make output directory
mkdir ./data/${OutputName}

# Format the contig clustering table to tab delimited
sed 's/,/\t/' ${ContigClusters} > ./data/${OutputName}/ContigClust.tsv

# Run the subroutines
# I know I know I should loop this
AnnotateCollapseClusters \
	${ContigAbundance} \
	${TransformedAbundanceFile}
