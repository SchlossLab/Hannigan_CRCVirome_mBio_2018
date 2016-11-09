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
sed -i 's/^ *//g' ${tmpdir}/rawcounts.tsv
sed -i 's/ */\t/g' ${tmpdir}/rawcounts.tsv
sed -i 's/\t.*\//\t/g' ${tmpdir}/rawcounts.tsv
sed -i 's/_R1.fastq//' ${tmpdir}/rawcounts.tsv
awk '{ print $1/2"\t"$2 }' > ${Output}
# rm -rf ${tmpdir}
