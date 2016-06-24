#! /bin/bash
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# IMPORTANT: This is quick and dirty. Proper way will be
# the ORF voting method.

#PBS -N AnnotateContigs
#PBS -q first
#PBS -l nodes=1:ppn=1,mem=40gb
#PBS -l walltime=500:00:00
#PBS -j oe
#PBS -V
#PBS -A schloss_lab

# Use this shell script to subsample and prepare C diff
# infection dataset for analysis on PMACS.

# Set the variables to be used in this script
export WorkingDirectory=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data
export Output='ContigAnnotations'

export FileSource=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/ContigsAndAlignments/NexteraXT002Contigs.fa
export ReferenceFile=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ConjunctisViribus/data/phageSVAnospace.fa

export GitBin=/mnt/EXT/Schloss-data/ghannig/OpenMetagenomeToolkit/pakbin
export SeqtkPath=/home/ghannig/bin/seqtk/seqtk
export LocalBin=/home/ghannig/bin/

###################
# Set Subroutines #
###################

AnnotateBlast () {
	# 1 = Reference
	# 2 = Query
	# 3 = SampleID

	makeblastdb \
		-dbtype nucl \
		-in ${1} \
		-out ./${Output}/PhageReference
    
    blastn \
    	-query ${2} \
    	-out ./${Output}/${3}-blastn.tsv \
    	-db ./${Output}/PhageReference \
    	-outfmt 6 \
    	-evalue 1e-5 \
    	-max_target_seqs 1
}

export -f AnnotateBlast

################
# Run Analysis #
################

AnnotateBlast \
	${ReferenceFile} \
	${FileSource} \
	"NexteraXT002"
