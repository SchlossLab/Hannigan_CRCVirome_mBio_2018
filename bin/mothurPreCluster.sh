#! /bin/bash
# mothurPreCluster.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export WORKDIR=$1

export MothurPath=/mnt/EXT/Schloss-data/bin/mothur

################
# Run Analysis #
################

echo PROGRESS: Formatting files.

cp ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta ${WORKDIR}/unmatched.fasta
cp ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table ${WORKDIR}/unmatched.count_table
cp ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy ${WORKDIR}/unmatched.taxonomy

# Move Error analysis to error directory
mkdir -p ${WORKDIR}/error_analysis/

cp ${WORKDIR}/stability.*.error.* ${WORKDIR}/error_analysis/
cp ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table ${WORKDIR}/error_analysis/
cp ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta ${WORKDIR}/error_analysis/

echo PROGRESS: Splitting file for clustering.

# Split the files for clustering

$MothurPath "#cluster.split(fasta=${WORKDIR}/unmatched.fasta, count=${WORKDIR}/unmatched.count_table, taxonomy=${WORKDIR}/unmatched.taxonomy, splitmethod=classify, taxlevel=5, cutoff=0.1, cluster=F, processors=4)"
