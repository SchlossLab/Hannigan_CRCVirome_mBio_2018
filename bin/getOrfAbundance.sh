#! /bin/bash
# getOrfAbundance.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#######################
# Set the Environment #
#######################
export OrfFile=$1
export FastaSequences=$2
export diamondpath=/mnt/EXT/Schloss-data/bin/
export tmpfile=./data/tmporfabund

export WorkingDirectory=$(pwd)

###################
# Set Subroutines #
###################
GetOrfHits () {
	sampleid=$(echo ${1} | sed 's/_2.fastq//')
	echo Samplie ID is ${sampleid}

	# Use blast to get hits of ORFs to Uniprot genes
	echo Running Phage ORFs...
	${diamondpath}diamond blastx \
		-q ${FastaSequences}/${1} \
		-d ${tmpfile}/diamonddatabase \
		-a ${tmpfile}/${sampleid}-output.daa \
		-t ./ \
		--max-target-seqs 1 \
		--evalue 1e-15 \
		--id 0.90

	${diamondpath}diamond view \
		-a ${tmpfile}/${sampleid}-output.daa \
		-o ${FastaSequences}/${1}.diamondout
}

export -f GetOrfHits

################
# Run Analysis #
################

# Create diamond database
echo Creating Database...
${diamondpath}diamond makedb \
	--in "${OrfFile}" \
	-d ${tmpfile}/diamonddatabase

ls ${FastaSequences}/*_R2.fastq | sed "s/.*\///g" | xargs -I {} --max-procs=32 bash -c 'GetOrfHits "$@"' _ {}


