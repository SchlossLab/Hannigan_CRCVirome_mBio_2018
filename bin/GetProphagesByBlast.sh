#! /bin/bash
# GetProphagesByBlast.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Pennsylvania
# NOTE: One way to predict whether a phage is associated with a
# bacterium is to determine whether the phage integrates into that
# bacterium. The simplest way to do this is simply using blast to
# determine whether the phage or it's genes are found within the
# bacterial host.

#######################
# Set the Environment #
#######################
export WorkingDirectory=${4}
export Output='tmp'
export BlastPath=${5}

export PhageGenomes=${1}
export BacteriaGenomes=${2}
export OutputFile=${3}

# Make the output directory and move to the working directory
echo Creating output directory...
cd "${WorkingDirectory}" || exit
mkdir ./${Output}

BlastPhageAgainstBacteria () {
	# 1 = Phage Genomes
	# 2 = Bacterial Genomes

	echo Making blast database...
	${BlastPath}makeblastdb \
		-dbtype nucl \
		-in "${2}" \
		-out ./${Output}/BacteraGenomeReference

	echo Running blastn...
	${BlastPath}blastn \
    	-query "${1}" \
    	-out ./${Output}/PhageToBacteria.blastn \
    	-db ./${Output}/BacteraGenomeReference \
    	-evalue 1e10 \
    	-num_threads 8 \
    	-outfmt 6

    echo Formatting blast output...
    # Get the Spacer ID, Phage ID, and BitScore
	cut -f 1,2,12 ./${Output}/PhageToBacteria.blastn \
		| sed 's/_\d\+\t/\t/' \
		> "${3}"
}

export -f BlastPhageAgainstBacteria

BlastPhageAgainstBacteria \
	"${PhageGenomes}" \
	"${BacteriaGenomes}" \
	"${OutputFile}"

# Remove the tmp output file
rm -r ./${Output}
