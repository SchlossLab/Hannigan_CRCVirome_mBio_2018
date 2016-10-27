#! /bin/bash
# ContigAssembly.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export FileR1=${1}
export FileR2=${2}
export OutputDirectory=${3}

export megahitvar=/home/ghannig/bin/megahit/megahit

###################
# Set Subroutines #
###################

PairedAssembleContigs () {
	echo Output is "${3}"
	echo First file is "${1}"
	echo Second file is "${2}"
	python ${megahitvar} \
		-1 "${1}" \
		-2 "${2}" \
		--min-contig-len 1000 \
		--k-min 21 \
		--k-max 101\
		--k-step 20 \
		-t 4 \
		-o "${3}"
}

export -f PairedAssembleContigs

################
# Run Analysis #
################

echo PROGRESS: Getting only sequence pairs for contigs
python ./bin/get_trimmed_pairs.py \
	-f ${FileR1} \
	-s ${FileR2} \
	-o ${FileR1}.paired \
	-t ${FileR2}.paired

echo PROGRESS: Performing contig assembly
PairedAssembleContigs \
	${FileR1}.paired \
	${FileR2}.paired \
	${OutputDirectory}

echo PROGRESS: Cleaning up intermediate files
rm -f ${FileR1}.paired
rm -f ${FileR2}.paired
