#! /bin/bash
# GetPairsByBlastx.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set the variables to be used in this script
export WorkingDirectory=${4}
export Output='tmp'

export SchlossBin=${5}

# Get the orfs that were already predicted in 'GerMicrobeOrfs.sh'
export PhageOrfs=${1}
export BacteriaOrfs=${2}
export OutputFile=${3}

# Make the output directory and move to the working directory
echo Creating output directory...
cd "${WorkingDirectory}" || exit
mkdir ./${Output}

echo Diamond is found in ${SchlossBin}...

GetPfamHits () {
	# 1 = Phage Orfs
	# 2 = Bacteria Orfs

	# Create diamond database
	echo Creating Bacteria Gene Database...
	${SchlossBin}diamond makedb \
		--in "${2}" \
		-d ./${Output}/DiamondReference

	# Use blast to get hits of ORFs to Uniprot genes
	${SchlossBin}diamond blastp \
		-q "${1}" \
		-d ./${Output}/DiamondReference \
		-a ./${Output}/Blastx.daa \
		-t ./

	${SchlossBin}diamond view \
		-a ./${Output}/Blastx.daa \
		-o "${OutputFile}"
}

export -f GetPfamHits

GetPfamHits \
	"${PhageOrfs}" \
	"${BacteriaOrfs}"
