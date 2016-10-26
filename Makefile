# Makefile
# Hannigan-2016-ColonCancerVirome
# Geoffrey Hannigan
# Pat Schloss Lab

###################
# Quality Control #
###################

SAMPLELIST_R1 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/_R1.fastq/' \
	| sed 's/^/data\/QC\//')

SAMPLELIST_R2 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/_R2.fastq/' \
	| sed 's/^/data\/QC\//')

DATENAME := $(shell date | sed 's/ /_/g' | sed 's/\:/\./g')

runqc: $(SAMPLELIST_R1) $(SAMPLELIST_R2)

$(SAMPLELIST_R1): data/QC/%_R1.fastq: data/raw/NexteraXT003/%_R1.fastq
	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/QC
	bash ./bin/QualityProcess.sh \
		$< \
		data/metadata/NexteraXT003Map.tsv \
		$@

$(SAMPLELIST_R2): data/QC/%_R2.fastq: data/raw/NexteraXT003/%_R2.fastq
	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/QC
	bash ./bin/QualityProcess.sh \
		$< \
		data/metadata/NexteraXT003Map.tsv \
		$@

#########################
# Human Decontamination #
#########################
DECON_R1 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/_R1.fastq/' \
	| sed 's/^/data\/HumanDecon\//')

DECON_R2 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/_R2.fastq/' \
	| sed 's/^/data\/HumanDecon\//')

humandeconseq: $(DECON_R1) $(DECON_R2)

$(DECON_R1): data/HumanDecon/%_R1.fastq: data/QC/%_R1.fastq
	echo $(shell date)  :  Performing HumanDecon and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/HumanDecon
	bash ./bin/HumanDeconSeq.sh \
		$< \
		$@

$(DECON_R2): data/HumanDecon/%_R2.fastq: data/QC/%_R2.fastq
	echo $(shell date)  :  Performing HumanDecon and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/HumanDecon
	bash ./bin/HumanDeconSeq.sh \
		$< \
		$@

###################
# Contig Assembly #
###################

CONTIGS_R1 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/.fastq/' \
	| sed 's/^/data\/contigs\//')

assemblecontigs: $(CONTIGS_R1)

$(CONTIGS_R1): data/contigs/%.fastq : data/HumanDecon/%_R1.fastq
	mkdir -p ./data/contigs
	bash ./bin/ContigAssembly.sh \
		$< \
		$(subst R1,R2,$<) \
		$@

#################
# Bacterial 16S #
#################
# Download the 16S reads from the SRA
mothurproc :
	mkdir -p ./data/mothur16S
	bash ./bin/Mothur16S.sh \
		./data/mothur16S \
		./data/raw/Zackular_16S

