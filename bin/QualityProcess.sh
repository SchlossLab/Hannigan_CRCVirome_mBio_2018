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

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export WorkingDirectory=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data
export Output='QualitySeqs'

export MappingFile=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/NexteraXT001Map.tsv

export GitBin=/home/ghannig/git/SchlossLab/bin/
export SeqtkPath=/home/ghannig/bin/seqtk/seqtk
export CutAdapt=/mnt/EXT/Schloss-data/bin/cutadapt-1.9.1/bin/cutadapt
export DeconsSeq=/mnt/EXT/Schloss-data/bin/deconseq-standalone-0.4.3/deconseq.pl
export fastx=/home/ghannig/bin/fastq_quality_trimmer

export RawSequenceDir=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/rawFastq/NexteraXT001

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

###################
# Set Subroutines #
###################
# Some of these look simple, but this is an easy way to ensure the parameters are standardized
# across multiple calls of all subroutines.

runCutadaptWithMap () {
	echo Input fastq = ${1} #1 = full path to fastq file to be trimmed
	echo Mapping file = ${2} #2 = full path to mapping file mapping file
	echo Output file = ${3} #3 = full path for output directory
	export SAMPLEID=$(echo ${1} | sed 's/_L001.*//g')
	export THREEPRIME=$(awk --assign sampleid=$SAMPLEID '$2 == sampleid { print $22 }' ${2})
	export FIVEPRIME=$(awk --assign sampleid=$SAMPLEID '$2 == sampleid { print $19 }' ${2})
	python2.7 ${CutAdapt} --error-rate=0.1 --overlap=10 -a $THREEPRIME -a $FIVEPRIME ${1} > ${3}
}

runFastx () {
	${fastx} -t 33 -Q 33 -l 150 -i ${1} -o ${2}
}

runDeconSeq () {
	# This is set for mouse decontamination
	perl ${DeconsSeq} -f ${1} -dbs mmref -out_dir ${2}
}

# Get them subroutines
export -f runCutadaptWithMap
export -f runFastx
export -f runDeconSeq

############
# Run Data #
############
for name in $(awk '{ print $2 }' ${MappingFile}); do
	# Because we are dealing with both directions
	for primer in R1 R2; do
		echo Processing sample ${name} and primer ${primer}...
		mkdir ./${Output}/CutAdapt
		runCutadaptWithMap \
			${RawSequenceDir}/${name}*${primer}*.fastq \
			${MappingFile} \
			./${Output}/CutAdapt/${name}_${primer}.fastq

		mkdir ./${Output}/FastxTrim
		runFastx \
			./${Output}/CutAdapt/${name}_${primer}.fastq \
			./${Output}/FastxTrim/${name}_${primer}.fastq

		mkdir ./${Output}/DeconSeq
		runDeconSeq \
			./${Output}/FastxTrim/${name}_${primer}.fastq \
			./${Output}/DeconSeq/${name}_${primer}.fastq
	done
done
