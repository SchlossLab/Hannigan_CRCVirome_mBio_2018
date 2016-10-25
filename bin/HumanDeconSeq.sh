#! /bin/bash
# QualityProcess.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export Inputfile=$1
export Outputfilename=$2

# Dependencies
export DeconsSeq=/mnt/EXT/Schloss-data/bin/deconseq-standalone-0.4.3/deconseq.pl
export fastx=/home/ghannig/bin/fastq_quality_trimmer
export CutAdapt=/mnt/EXT/Schloss-data/bin/cutadapt-1.9.1/bin/cutadapt

###################
# Set Subroutines #
###################
# Some of these look simple, but this is an easy way to ensure the parameters are standardized
# across multiple calls of all subroutines.

runDeconSeq () {
	# This is set for human decontamination
	echo Input is "${1}"
	echo Output is "${2}"
	echo Output filename is "${3}"
	echo Contaminated reads "${4}"

	perl ${DeconsSeq} -f "${1}" -dbs hsref -out_dir "${2}"
	cp "${2}"/*clean.fq "${3}"
	cp "${2}"/*cont.fq "${4}"
	# rm -rf "${2}"
}

export -f runDeconSeq

############
# Run Data #
############

echo PROGRESS: Decontaminating human reads

runDeconSeq \
	${Inputfile} \
	${Inputfile}.output \
	${Outputfilename} \
	${Outputfilename}.cont
