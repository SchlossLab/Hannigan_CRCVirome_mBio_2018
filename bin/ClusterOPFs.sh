#! /bin/bash
# ClusterOPFs.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#######################
# Set the Environment #
#######################
export ContigsFile=${1}
export Output=${2}

export WorkingDirectory=$(pwd)

mkdir ./data/tmp-opfs

###################
# Set Subroutines #
###################

PredictOrfs () {
	# 1 = Contig Fasta File for Prodigal
	# 2 = Output File Name

	bash ./bin/ProdigalWrapperLargeFiles.sh \
		"${1}" \
		./data/tmp-opfs/tmp-genes.fa

    # Remove the block formatting
	perl \
	./bin/remove_block_fasta_format.pl \
		./data/tmp-opfs/tmp-genes.fa \
		"${2}"

	sed -i 's/\*//g' "${2}"
	sed -i 's/\/n//g' "${2}"
}

EstablishOpfs () {
	# 1 = Open Reading Frame fasta

	# Set MMseqs variables
	export MMDIR=/home/ghannig/bin/mmseqs2
	PATH=$MMDIR/bin:$PATH
	echo Path is "$PATH"

	cd ./data/tmp-opfs || exit

	# Create database
	mmseqs createdb ${1} DB

	mkdir ./tmp
    mmseqs clustlinear DB clu tmp

    # Convert to fasta
    mmseqs createseqfiledb DB clu clu_seq 
    mmseqs result2flat DB DB clu_seq clu_seq.fasta

    # Back out of the directory
    cd ../.. || exit
    cp ./data/tmp-opfs/clu_seq.fasta ${2}
}

# Export the subroutines
export -f PredictOrfs
export -f EstablishOpfs


################
# Predict ORFs #
################
echo Predicting ORFs...

# PredictOrfs \
# 	${ContigsFile} \
# 	./data/tmp-opfs/ContigOrfs.fa \
# 	|| exit

EstablishOpfs \
	ContigOrfs.fa \
	${Output}
