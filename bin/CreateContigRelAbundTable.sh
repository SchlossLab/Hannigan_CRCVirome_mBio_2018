#! /bin/bash
# CreateContigRelAbundTable.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#######################
# Set the Environment #
#######################
export BowtieReference=$1
export FastaSequences=$2
export Output='data/tmpbowtie'

mkdir -p ./${Output}
mkdir -p ./${Output}/bowtieReference


###################
# Set Subroutines #
###################
GetHits () {
	# 1 = Input Orfs
	# 2 = Bowtie Reference

	bowtie2 \
		-x ${2} \
		-q ${1} \
		-S ${1}-bowtie.sam \
		-p 8 \
		-L 25 \
		-N 1

	# Quantify alignment hits
	perl \
		./bin/calculate_abundance_from_sam.pl \
			${1}-bowtie.sam \
			${1}-bowtie.tsv
}

BowtieRun () {
	sampleid=$(echo ${1} | sed 's/_2.fastq//')
	GetHits \
		${1} \
		${2}

	# Remove the header
	sed -e "1d" ${1}-bowtie.tsv > ${1}-noheader

	awk -v name=${sampleid} '{ print $0"\t"name }' ${1}-noheader \
	| grep -v '\*' > ${1}-noheader-forcat
	# rm ${1}/${1}-noheader
}

# Export the subroutines
export -f GetHits
export -f BowtieRun

#############################
# Contig Relative Abundance #
#############################

echo Running bowtie...

BowtieRun ${FastaSequences} ${BowtieReference}

