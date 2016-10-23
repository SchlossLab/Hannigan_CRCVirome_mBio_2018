# Makefile
# Hannigan-2016-ColonCancerVirome
# Geoffrey Hannigan
# Pat Schloss Lab

#
# Set General Variables 
#

SAMPLELIST := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/_R1.fastq/' \
	| sed 's/^/data\/QC\//')
DATENAME := $(shell date | sed 's/ /_/g' | sed 's/\:/\./g')

objects = foo.o bar.o

all: $(objects)

$(objects): %.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

print: $(SAMPLELIST)

# $(SAMPLELIST) : %.txt :
# 	touch $@

###################
# Quality Control #
###################

$(SAMPLELIST): data/QC/%_R1.fastq: data/raw/NexteraXT003/%_R1.fastq.gz
	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/QC
	bash ./bin/QualityProcess.sh \
		$< \
		data/metadata/NexteraXT003Map.tsv \
		$@
