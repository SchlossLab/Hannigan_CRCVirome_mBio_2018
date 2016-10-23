# Makefile
# Hannigan-2016-ColonCancerVirome
# Geoffrey Hannigan
# Pat Schloss Lab

#
# Set General Variables 
#

SAMPLELIST := $(shell awk '{ print $$3 }' ./data/PublishedDatasets/metadatatable.tsv | sort | uniq)
DATENAME := $(shell date | sed 's/ /_/g' | sed 's/\:/\./g')

###################
# Quality Control #
###################

${SAMPLELIST}: ./data/QC/%_R1.fastq ./data/QC/%_R2.fastq: \
			./data/raw/NexteraXT003/%_R1.fastq.gz \
			./data/raw/NexteraXT003/%_R2.fastq.gz \
			./data/metadata/NexteraXT003Map.tsv \
			./bin/QualityProcess.sh
	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir ./data/QC
	bash ./bin/QualityProcess.sh \
		$<1 \
		./data/metadata/NexteraXT003Map.tsv \
		$@1
	bash ./bin/QualityProcess.sh \
		$<2 \
		./data/metadata/NexteraXT003Map.tsv \
		$@2
