#! /bin/bash
# getPairedCount.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#######################
# Set the Environment #
#######################
export SequencingFileDirectory=${1}
export Output=${2}
export tmpdir=./data/tmp-counts

mkdir -p ${tmpdir}

wc -l ${SequencingFileDirectory}/*_R1.fastq > ${tmpdir}/rawcounts.tsv
sed 's/^ *//g' ${tmpdir}/rawcounts.tsv \
	| grep -v total \
	| sed 's/ \+/\t/g' \
	| sed 's/\t.*\//\t/g' \
	| sed 's/_R1.fastq//' \
	| awk '{ print $1/2"\t"$2 }' \
	> ${Output}

# rm -rf ${tmpdir}
