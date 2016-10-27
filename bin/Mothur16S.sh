#! /bin/bash
# Mothur16S.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export WORKDIR=$1
export SampleDir=$2

export MothurPath=/mnt/EXT/Schloss-data/bin/mothur

################
# Run Analysis #
################

echo PROGRESS: Copying samples to new directory.
for file in "${SampleDir}"/*; do
	filename=$(echo ${file} | sed 's/.*\///g' | sed 's/\-/_/g')
	cp ${file} ${WORKDIR}/${filename}
done

echo PROGRESS: Running Mothur.

# Section from Mothur MiSeq SOP
# Adapted from Baxter, Sze, Schloss
$MothurPath "#make.file(inputdir=${WORKDIR});
	make.contigs(file=current, processors=8);
	summary.seqs(fasta=current);
	screen.seqs(fasta=current, group=current, summary=current, maxambig=0, maxlength=275);
	summary.seqs(fasta=current);
	unique.seqs(fasta=current);
	summary.seqs(fasta=current, name=current);
	count.seqs(name=current, group=current);
	summary.seqs(fasta=current, count=current);
	align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=data/references/silva.v4.align);
	summary.seqs(fasta=current, count=current);
	screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8);
	summary.seqs(fasta=current,count=current);
	filter.seqs(fasta=data/mothur16S/stability.trim.contigs.good.unique.good.align, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.uchime(fasta=current, count=current, dereplicate=t);
	remove.seqs(fasta=current, accnos=current);
	summary.seqs(fasta=current,count=current);
	classify.seqs(fasta=current, count=current, reference=data/references/trainset14_032015.pds.fasta, taxonomy=data/references/trainset14_032015.pds.tax, cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	get.groups(fasta=current, count=current, groups=mock1-mock2);
	seq.error(fasta=current, count=current, reference=$REF/HMP_MOCK.fasta, aligned=F)"

echo PROGRESS: Formatting files.

mv ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta ${WORKDIR}/unmatched.fasta
mv ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table ${WORKDIR}/unmatched.count_table
mv ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy ${WORKDIR}/unmatched.taxonomy

# Move Error analysis to error directory

mv ${WORKDIR}/stability.*.error.* ${WORKDIR}/error_analysis/
mv ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table ${WORKDIR}/error_analysis/
mv ${WORKDIR}/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta ${WORKDIR}/error_analysis/

echo PROGRESS: Splitting file for clustering.

# Split the files for clustering

$mothurRv "#cluster.split(fasta=${WORKDIR}/unmatched.fasta, count=${WORKDIR}/unmatched.count_table, taxonomy=${WORKDIR}/unmatched.taxonomy, splitmethod=classify, taxlevel=5, cutoff=0.1, cluster=F, processors=4)"

