#! /bin/bash
# PfamDomainInteractPrediction.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set the variables to be used in this script
export WorkingDirectory=${4}
export Output='tmp'

export PfamLocation=${6}
export PfamDatabase=${PfamLocation}Pfam-A-diamond
export PfamConversion=${PfamLocation}PfamAccToPF.tsv
export InteractionReference=${PfamLocation}PfamAccInteractions.tsv

export SchlossBin=${5}

# Get the orfs that were already predicted in 'GerMicrobeOrfs.sh'
export PhageOrfs=${1}
export BacteriaOrfs=${2}
export OutputFile=${3}

# Make the output directory and move to the working directory
echo Creating output directory...
cd "${WorkingDirectory}" || exit
mkdir ./${Output}

GetPfamHits () {
	# 1 = Database
	# 2 = Phage Orfs
	# 3 = Bacteria Orfs

	# Use blast to get hits of ORFs to Uniprot genes
	${SchlossBin}diamond blastp \
		-q "${2}" \
		-d "${1}" \
		-a ./${Output}/Phage.daa \
		-t ./

	${SchlossBin}diamond blastp \
		-q "${3}" \
		-d "${1}" \
		-a ./${Output}/Bacteria.daa \
		-t ./

	${SchlossBin}diamond view \
		-a ./${Output}/Phage.daa \
		-o ./${Output}/PhageBlast.txt

	${SchlossBin}diamond view \
		-a ./${Output}/Bacteria.daa \
		-o ./${Output}/BacteriaBlast.txt
}

OrfInteractionPairs () {
	# 1 = Phage Blast Results
	# 2 = Bacterial Blast Results
	# 3 = Interaction Reference
	# 4 = Acc to Pfam Conversion Table

	# Reverse the interaction reference for awk
	awk \
		'{ print $2"\t"$1 }' \
		"${3}" \
		> "${3}".inverse

	cat \
		"${3}" \
		"${3}".inverse \
		> ./${Output}/TotalInteractionRef.tsv

	# Get only the ORF IDs and corresponding interactions
	# Column 1 is the ORF ID, two is Uniprot ID
	cut -f 12,1,2 "${1}" | sed 's/\S\+|\(\S\+\)|\S\+$/\1/' | sed 's/\/.*\t/\t/' > ./${Output}/PhageBlastIdReference.tsv
	cut -f 12,1,2 "${2}" | sed 's/\S\+|\(\S\+\)|\S\+$/\1/' | sed 's/\/.*\t/\t/' > ./${Output}/BacteriaBlastIdReference.tsv

	# Convert the acc numbers to pfam IDs
	awk \
		'NR == FNR { a[$1] = $2; next } { print $1"\t"a[$2]"\t"$3 }' \
		"${4}" \
		./${Output}/PhageBlastIdReference.tsv \
	> ./${Output}/PhageBlastIdReferencePfams.tsv

	awk \
		'NR == FNR { a[$1] = $2; next } { print $1"\t"a[$2]"\t"$3 }' \
		"${4}" \
		./${Output}/BacteriaBlastIdReference.tsv \
	> ./${Output}/BacteriaBlastIdReferencePfams.tsv
	
	# Convert bacterial file to reference
	awk \
		'NR == FNR {a[$2] = $1"\t"$3; next} $1 in a { print $1"\t"$2"\t"a[$1] }' \
		./${Output}/PhageBlastIdReferencePfams.tsv \
		./${Output}/TotalInteractionRef.tsv \
		> ./${Output}/tmpMerge.tsv

	awk \
		'NR == FNR {a[$2] = $1"\t"$3; next} $2 in a { print $1"\t"$2"\t"$3"\t"$4"\t"a[$2] }' \
		./${Output}/BacteriaBlastIdReferencePfams.tsv \
		./${Output}/tmpMerge.tsv \
		| cut -f 3,4,5,6 \
		> "${OutputFile}"

	# This output can be used for input into perl script for adding
	# to the graph database.
}

export -f GetPfamHits
export -f OrfInteractionPairs

GetPfamHits \
	${PfamDatabase} \
	"${PhageOrfs}" \
	"${BacteriaOrfs}"

OrfInteractionPairs \
	./${Output}/PhageBlast.txt \
	./${Output}/BacteriaBlast.txt \
	${InteractionReference} \
	${PfamConversion}
