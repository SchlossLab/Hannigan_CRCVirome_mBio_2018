#! /bin/bash
# CreateContigRelAbundTable.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#######################
# Set the Environment #
#######################
export ContigsFile=$1
export FastaSequences=$2
export MasterOutput=$3
export Output='data/tmpbowtie'

mkdir -p ./${Output}

###################
# Set Subroutines #
###################
GetHits () {
	# 1 = Input Orfs
	# 2 = Bowtie Reference

	mkdir -p ./${Output}/bowtieReference

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
		${FastaSequences}/${1} \
		./${Output}/bowtieReference/bowtieReference

	# Remove the header
	sed -e "1d" ${FastaSequences}/${1}-bowtie.tsv > ${FastaSequences}/${1}-noheader

	awk -v name=${sampleid} '{ print $0"\t"name }' ${FastaSequences}/${1}-noheader \
	| grep -v '\*' > ${FastaSequences}/${1}-noheader-forcat
	# rm ${FastaSequences}/${1}-noheader
}

# Export the subroutines
export -f GetHits
export -f BowtieRun

#############################
# Contig Relative Abundance #
#############################

echo Getting contig relative abundance table...

# Clear the file to prepare for appending to new file below
rm -f ${MasterOutput}

# Build bowtie reference
bowtie2-build \
	-q ${ContigsFile} \
	./${Output}/bowtieReference/bowtieReference

ls ${FastaSequences}/*_R2.fastq | sed "s/.*\///g" | xargs -I {} --max-procs=32 bash -c 'BowtieRun "$@"' _ {}

echo Catting files...

cat ${FastaSequences}/*-noheader-forcat > ${MasterOutput}
