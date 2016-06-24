#! /bin/bash
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#PBS -N ContigAssemblyAndAlignment
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
export Output='ContigsAndAlignments'

export FileSource=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/QualitySeqs/DeconSeq
export MappingFile=/mnt/EXT/Schloss-data/ghannig/Hannigan-2016-ColonCancerVirome/data/NexteraXT002Map.tsv

export GitBin=/mnt/EXT/Schloss-data/ghannig/OpenMetagenomeToolkit/pakbin
export SeqtkPath=/home/ghannig/bin/seqtk/seqtk
export LocalBin=/home/ghannig/bin/

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

###################################
# Get Forward Reads and Subsample #
###################################

Subsampler () {
	# Convert the fastq to fasta
	${LocalBin}idba-1.1.1/bin/fq2fa ${1} ./tmp.fa

	# Subsample the sequences
	${SeqtkPath} sample ./tmp.fa 60000 > ${2}

	# Remove samples that do no meet subsampling criteria
	export WordCount=$(wc -l ${2} | sed 's/ .*//')
	echo Count is ${WordCount}
	if [ ${WordCount} -ne 120000 ]
		then
		echo Remove sample ${2} because too few sequences...
		rm ${2}
	fi
	# Remove tmp file
	rm ./tmp.fa
}

AssembleContigs () {
	# Assemble the contigs
	${LocalBin}idba-1.1.1/bin/idba_ud \
		--pre_correction \
		-l ${1} \
		-o ./${Output}/TotalContigs

	# Pull out the contig fasta files
	mv ./${Output}/TotalContigs/contig.fa ./${Output}/${2}.fa

	# Cut down the names of the contigs
	sed -i 's/ length.*//' ./${Output}/${2}.fa

	# Calculate contig stats
	perl ${GitBin}/CalculateContigStats.pl \
		./${Output}/${2}.fa \
		./${Output}/${2}-ContigStats.tsv

	sed -i 's/>//g' ./${Output}/${2}-ContigStats.tsv

	# Make master contig ID list
	sed -n 1~2p ./${Output}/${2}.fa \
	| sed s'/>//g' \
	| sed '1 s/^/Contig_ID\n/' \
	> ./${Output}/${2}-MasterList.tsv
}

BowtieAlignment () {
	mkdir ./${Output}/bowtieReference

	bowtie2-build \
		-f ${1} \
		./${Output}/bowtieReference/bowtieReference

	bowtie2 \
		-x ./${Output}/bowtieReference/bowtieReference \
		-f ${2} \
		-S ./${Output}/tmp-bowtie.sam \
		-p 32 \
		-L 25 \
		-N 1

	# Quantify alignment hits
	perl ${GitBin}/calculate_abundance_from_sam.pl ./${Output}/tmp-bowtie.sam ${3}

	# Remove the intermediate files
	rm -r ./${Output}/bowtieReference/
	rm ./${Output}/tmp-bowtie.sam
}

CalculateRelativeAbundance () {
	# 1 = Treatment name
	# 2 = Contig stats table
	# 3 = Relative abundance tsv from sam

	# Important to get length sum information
	SUM=$(awk '{ SUM += $2 } END { print SUM }' ${2})
	echo Sum of contigs is ${SUM}...

	# Calculate rpkm
	awk \
		--assign=sum=$SUM \
		'FNR==NR { a[$1]=$2; next }	$1 in a	{ print $1"\t"$2"\t"a[$1]"\t"$2*1000000000/(a[$1]*sum) }' \
		${2} \
		${3} \
		> ./${Output}/${1}-rpkm.tsv

	cut -f 1,4 ./${Output}/${1}-rpkm.tsv \
	| sed "1 s/^/Contig_ID\t${1}\n/" \
	> ./${Output}/${1}-rpkm-parsed.tsv

	awk \
		'FNR==NR {a[$1]=$2;next}\
		{ print $1"\t"a[$1] }' \
		./${Output}/${1}-rpkm-parsed.tsv \
		./${Output}/NexteraXT002Contigs-MasterList.tsv \
	| sed '/[0-9]\t[0-9]/!s/$/0/' \
	| sed '1 s/\t0//' \
	| sed '1 s/0$//' \
	> ./${Output}/${1}-AbundanceOnMaster.tsv

	# Get only the abundance values for merging in the next step
	cut -f 2 \
		./${Output}/${1}-AbundanceOnMaster.tsv \
		> ./${Output}/${1}-AbundanceOnMasterForMerge.tsv

	# After this runs, I need to paste the files to complete the table
}

export -f Subsampler
export -f AssembleContigs
export -f BowtieAlignment
export -f CalculateRelativeAbundance

###########
Run Data #
###########
for name in $(awk '{ print $2 }' ${MappingFile}); do
	# Because we are dealing with both directions
	for primer in R1 R2; do
		echo Sumsampling sequences from ${name}_${primer}...
		mkdir ./${Output}/SubsampledFasta
		Subsampler \
			${FileSource}/${name}_${primer}_clean.fastq \
			./${Output}/SubsampledFasta/${name}_${primer}.fa
	done
done

echo Assembling contigs...
cat ./${Output}/SubsampledFasta/* > ./${Output}/TotalSeqs.fa
AssembleContigs \
	./${Output}/TotalSeqs.fa \
	"NexteraXT002Contigs" \

mkdir ./${Output}/BowtieOutput

for name in $(awk '{ print $2 }' ${MappingFile}); do
	# Because we are dealing with both directions
	for primer in R2; do
		echo Aligning reads from ${name}_${primer}...
		BowtieAlignment \
			./${Output}/NexteraXT002Contigs.fa \
			./${Output}/SubsampledFasta/${name}_${primer}.fa \
			./${Output}/BowtieOutput/${name}_${primer}.tsv

		CalculateRelativeAbundance \
			${name}_${primer} \
			./${Output}/NexteraXT002Contigs-ContigStats.tsv \
			./${Output}/BowtieOutput/${name}_${primer}.tsv
	done
done


paste \
	./${Output}/NexteraXT002Contigs-MasterList.tsv \
	./${Output}/*-AbundanceOnMasterForMerge.tsv \
	> ./${Output}/ContigRelativeAbundanceTable.tsv
