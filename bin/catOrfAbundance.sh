#! /bin/bash
# catOrfAbundance.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#######################
# Set the Environment #
#######################
export FastaSequences=$1
export MasterOutput=$2
export Output='data/tmporfabund'

mkdir -p ./${Output}

#################
# Cat Abundance #
#################
rm -f ${MasterOutput}

for file in $2/*.diamondout; do
	sampleID=$(echo $file | sed 's/.*\///g' | sed 's/_R2\.fastq\.diamondout//')
	echo Sample ID is ${sampleID}
	#Save a file with abundance \t ORF ID
	cut -f 2 ${file} \
		| sort \
		| uniq -c \
		| sed 's/^ *//' \
		| sed 's/ /\t/' \
		| sed "s/$/\t$sampleID/" \
		>> ${MasterOutput}
done
