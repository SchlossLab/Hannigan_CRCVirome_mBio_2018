#! /bin/bash
# QualityProcess.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#PBS -N QualityProcess
#PBS -q first
#PBS -l nodes=1:ppn=1,mem=40gb
#PBS -l walltime=500:00:00
#PBS -j oe
#PBS -V
#PBS -A schloss_lab

############################
# Load in Required Modules #
############################
# Setup the R environment
module load R/3.2.3

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export WorkingDirectory=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data
export Output='QualitySeqs'
export MappingFile=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/NexteraXT002Map.tsv

# Dependencies
export LocalBin=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/bin
export DeconsSeq=/mnt/EXT/Schloss-data/bin/deconseq-standalone-0.4.3/deconseq.pl
export fastx=/home/ghannig/bin/fastq_quality_trimmer
export CutAdapt=/mnt/EXT/Schloss-data/bin/cutadapt-1.9.1/bin/cutadapt

export RawSequenceDir=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/raw/NexteraXT002

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory} || exit
mkdir ./${Output}

###################
# Set Subroutines #
###################
# Some of these lo
