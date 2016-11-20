# Makefile
# Hannigan-2016-ColonCancerVirome
# Geoffrey Hannigan
# Pat Schloss Lab

############################################# METADATA ############################################
metadatafiles = ./data/metadata/NexteraXT003Map.tsv ./data/metadata/NexteraXT004Map.tsv

# Make a master metadata file
./data/metadata/MasterMeta.tsv : $(metadatafiles)
	cat ./data/metadata/NexteraXT003Map.tsv ./data/metadata/NexteraXT004Map.tsv > ./data/metadata/MasterMeta.tsv

######################################### QUALITY CONTROL #########################################

###################
# Quality Control #
###################
setfile1: ./data/metadata/MasterMeta.tsv
	$(eval SAMPLELIST_R1 := $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R1.fastq/' \
		| sed 's/^/data\/QC\//'))
	echo 'variable1 = $(SAMPLELIST_R1)' > $@

setfile2: ./data/metadata/MasterMeta.tsv
	$(eval SAMPLELIST_R2 = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R2.fastq/' \
		| sed 's/^/data\/QC\//'))
	echo 'variable2 = $(SAMPLELIST_R2)' > $@

DATENAME := $(shell date | sed 's/ /_/g' | sed 's/\:/\./g')

runqc:

include setfile1
include setfile2

runqc: $(variable1) $(variable2)

$(variable): data/QC/%_R1.fastq:
	echo Target is $@

$(variable1): data/QC/%_R1.fastq: data/raw/hiseqcat/%_R1.fastq
	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/QC
	bash ./bin/QualityProcess.sh \
		$< \
		data/metadata/MasterMeta.tsv \
		$@

$(variable2): data/QC/%_R2.fastq: data/raw/hiseqcat/%_R2.fastq
	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/QC
	bash ./bin/QualityProcess.sh \
		$< \
		data/metadata/MasterMeta.tsv \
		$@

#########################
# Human Decontamination #
#########################
setfile3: ./data/metadata/MasterMeta.tsv
	$(eval DECON_R1 = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R1.fastq/' \
		| sed 's/^/data\/HumanDecon\//'))
	echo 'variable3 = $(DECON_R1)' > $@

setfile4: ./data/metadata/MasterMeta.tsv
	$(eval DECON_R2 = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R2.fastq/' \
		| sed 's/^/data\/HumanDecon\//'))
	echo 'variable4 = $(DECON_R2)' > $@

humandeconseq:

include setfile3
include setfile4

humandeconseq: $(variable3) $(variable4)

$(variable3): data/HumanDecon/%_R1.fastq: data/QC/%_R1.fastq
	echo $(shell date)  :  Performing HumanDecon and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/HumanDecon
	bash ./bin/HumanDeconSeq.sh \
		$< \
		$@

$(variable4): data/HumanDecon/%_R2.fastq: data/QC/%_R2.fastq
	echo $(shell date)  :  Performing HumanDecon and contig alignment on samples $@ >> ${DATENAME}.makelog
	mkdir -p ./data/HumanDecon
	bash ./bin/HumanDeconSeq.sh \
		$< \
		$@

###################
# Contig Assembly #
###################
setfile5: ./data/metadata/MasterMeta.tsv
	$(CONTIGS_R1 = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/.fastq/' \
		| sed 's/^/data\/contigs\//'))
	echo 'variable5 = $(CONTIGS_R1)' > $@

setfile6: ./data/metadata/MasterMeta.tsv
	$(MOVE_CONTIGS = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/.fastq/' \
		| sed 's/^/data\/contigfastq\//'))
	echo 'variable6 = $(MOVE_CONTIGS)' > $@

assemblecontigs:

include setfile5

assemblecontigs: $(variable5)

movecontigs:

include setfile6

movecontigs: $(variable6)

$(variable5): data/contigs/%.fastq : data/HumanDecon/%_R1.fastq
	mkdir -p ./data/contigs
	bash ./bin/ContigAssembly.sh \
		$< \
		$(subst R1,R2,$<) \
		$@

$(variable6): data/contigfastq/%.fastq :
	mkdir -p data/contigfastq
	cp $(subst contigfastq,contigs,$@)/final.contigs.fa $@

############################################# CONTIGS #############################################
# Get a master file for bacterial and viral contigs
VIRUS_CONTIGS := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/.fastq/' \
	| sed 's/^/data\/makerecorddump\//')

BACTERIA_CONTIGS := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT004Map.tsv \
	| sort \
	| uniq \
	| sed 's/$$/.fastq/' \
	| sed 's/^/data\/makerecorddump\//')

contigpairs: $(VIRUS_CONTIGS) $(BACTERIA_CONTIGS)

$(VIRUS_CONTIGS): data/makerecorddump/%.fastq : data/contigfastq/%.fastq
	mkdir -p data/makerecorddump
	touch $@
	cat $< >> ./data/totalcontigsvirus.fa
	sed -i 's/ /_/g' ./data/totalcontigsvirus.fa

$(BACTERIA_CONTIGS): data/makerecorddump/%.fastq : data/contigfastq/%.fastq
	mkdir -p data/makerecorddump
	touch $@
	cat $< >> ./data/totalcontigsbacteria.fa
	sed -i 's/ /_/g' ./data/totalcontigsbacteria.fa

# ###############################
# # Contig Abundance Per Sample #
# ###############################
# ./data/ContigRelAbundForGraph.tsv : \
# 			./data/totalcontigs.fa
# 	bash ./bin/CreateContigRelAbundTable.sh \
# 		./data/totalcontigs.fa \
# 		./data/HumanDecon \
# 		./data/ContigRelAbundForGraph.tsv

# #####################
# # Contig Clustering #
# #####################
# ./data/ContigAbundForConcoct.tsv : ./data/ContigRelAbundForGraph.tsv
# 	Rscript ./bin/ReshapeAlignedAbundance.R \
# 		-i ./data/ContigRelAbundForGraph.tsv \
# 		-o ./data/ContigAbundForConcoct.tsv \
# 		-p 0.5

# ./data/ContigClusters/clustering_gt1000.csv : \
# 			./data/totalcontigs.fa \
# 			./data/ContigAbundForConcoct.tsv
# 	mkdir ./data/ContigClustersPhage
# 	concoct \
# 		--coverage_file ./data/ContigAbundForConcoct.tsv \
# 		--composition_file ./data/totalcontigs.fa \
# 		--clusters 500 \
# 		--kmer_length 4 \
# 		--length_threshold 1000 \
# 		--read_length 150 \
# 		--basename ./data/ContigClusters/ \
# 		--no_total_coverage \
# 		--iterations 50

# ############################
# # Cluster Contig Abundance #
# ############################
# ./data/ClusteredContigAbund.tsv : \
# 			./data/ContigRelAbundForGraph.tsv \
# 			./data/ContigClusters/clustering_gt1000.csv
# 	bash ./bin/ClusterContigAbund.sh \
# 		./data/ContigRelAbundForGraph.tsv \
# 		./data/ContigClusters/clustering_gt1000.csv \
# 		./data/ClusteredContigAbund.tsv

# #################
# # Identify OPFs #
# #################
# ./data/totalopfs.fa : ./data/totalcontigs.fa
# 	bash ./bin/ClusterOPFs.sh \
# 		$< \
# 		$@

# ############################
# # ORF Abundance Per Sample #
# ############################
# orfalign:
# 	bash ./bin/getOrfAbundance.sh \
# 		./data/tmp-opfs/ContigOrfsNoSpec.fa \
# 		./data/HumanDecon

# ./data/orfabund.tsv:
# 	bash ./bin/catOrfAbundance.sh \
# 		./data/HumanDecon \
# 		./data/orfabund.tsv

# #######################
# # OPF Abundance Table #
# #######################
# ./data/ClusteredOpfAbund.tsv : \
# 			./data/totalopfs.fa.tsv \
# 			./data/orfabund.tsv
# 	$(shell awk '{ print $$2"\t"$$1 }' ./data/totalopfs.fa.tsv > ./data/OpfClusters.tsv)
# 	$(shell awk '{ print $$2"\t"$$1"\t"$$3 }' ./data/orfabund.tsv > ./data/orfabundOrdered.tsv)
# 	bash ./bin/ClusterContigAbund.sh \
# 		./data/orfabundOrdered.tsv \
# 		./data/OpfClusters.tsv \
# 		./data/ClusteredOpfAbund.tsv

# #######################
# # OGU Alpha Diversity #
# #######################
# ./figures/diversity-alpha-ogu.pdf : \
# 			./data/ClusteredContigAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-alpha.R \
# 		--input ./data/ClusteredContigAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--out ./figures/diversity-alpha-ogu.pdf

# #######################
# # OPF Alpha Diversity #
# #######################
# ./figures/diversity-alpha-opf.pdf : \
# 			./data/ClusteredOpfAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-alpha.R \
# 		--input ./data/ClusteredOpfAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--out ./figures/diversity-alpha-opf.pdf

# #######################
# # OGU Beta Diversity #
# #######################
# ./figures/diversity-beta-ogu.pdf : \
# 			./data/ClusteredContigAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredContigAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--out ./figures/diversity-beta-ogu.pdf

# ./figures/diversity-betajaccard-ogu.pdf : \
# 			./data/ClusteredContigAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredContigAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--divmetric jaccard \
# 		--out ./figures/diversity-betajaccard-ogu.pdf

# #######################
# # OPF Beta Diversity #
# #######################
# ./figures/diversity-beta-opf.pdf : \
# 			./data/ClusteredOpfAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredOpfAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 50000 \
# 		--out ./figures/diversity-beta-opf.pdf

# ./figures/diversity-betajaccard-opf.pdf : \
# 			./data/ClusteredOpfAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredOpfAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 50000 \
# 		--divmetric jaccard \
# 		--out ./figures/diversity-betajaccard-opf.pdf

# ####################
# # Sequencing Depth #
# ####################
# ./data/ProjectSeqDepth.tsv :
# 	bash ./bin/getPairedCount.sh \
# 		./data/HumanDecon \
# 		./data/ProjectSeqDepth.tsv

# ./figures/qualitycontrol.pdf :
# 	Rscript ./bin/QualityControlStats.R \
# 		--input ./data/ProjectSeqDepth.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--out ./figures/qualitycontrol.pdf \
# 		--sdepth 100000

# ######################################## WHOME METAGENOME #########################################
# ###################
# # Quality Control #
# ###################

# SAMPLELIST_R1 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
# 	| sort \
# 	| uniq \
# 	| sed 's/$$/_R1.fastq/' \
# 	| sed 's/^/data\/QC\//')

# SAMPLELIST_R2 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
# 	| sort \
# 	| uniq \
# 	| sed 's/$$/_R2.fastq/' \
# 	| sed 's/^/data\/QC\//')

# DATENAME := $(shell date | sed 's/ /_/g' | sed 's/\:/\./g')

# runqc: $(SAMPLELIST_R1) $(SAMPLELIST_R2)

# $(SAMPLELIST_R1): data/QC/%_R1.fastq: data/raw/NexteraXT003/%_R1.fastq
# 	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
# 	mkdir -p ./data/QC
# 	bash ./bin/QualityProcess.sh \
# 		$< \
# 		data/metadata/NexteraXT003Map.tsv \
# 		$@

# $(SAMPLELIST_R2): data/QC/%_R2.fastq: data/raw/NexteraXT003/%_R2.fastq
# 	echo $(shell date)  :  Performing QC and contig alignment on samples $@ >> ${DATENAME}.makelog
# 	mkdir -p ./data/QC
# 	bash ./bin/QualityProcess.sh \
# 		$< \
# 		data/metadata/NexteraXT003Map.tsv \
# 		$@

# #########################
# # Human Decontamination #
# #########################
# DECON_R1 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
# 	| sort \
# 	| uniq \
# 	| sed 's/$$/_R1.fastq/' \
# 	| sed 's/^/data\/HumanDecon\//')

# DECON_R2 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
# 	| sort \
# 	| uniq \
# 	| sed 's/$$/_R2.fastq/' \
# 	| sed 's/^/data\/HumanDecon\//')

# humandeconseq: $(DECON_R1) $(DECON_R2)

# $(DECON_R1): data/HumanDecon/%_R1.fastq: data/QC/%_R1.fastq
# 	echo $(shell date)  :  Performing HumanDecon and contig alignment on samples $@ >> ${DATENAME}.makelog
# 	mkdir -p ./data/HumanDecon
# 	bash ./bin/HumanDeconSeq.sh \
# 		$< \
# 		$@

# $(DECON_R2): data/HumanDecon/%_R2.fastq: data/QC/%_R2.fastq
# 	echo $(shell date)  :  Performing HumanDecon and contig alignment on samples $@ >> ${DATENAME}.makelog
# 	mkdir -p ./data/HumanDecon
# 	bash ./bin/HumanDeconSeq.sh \
# 		$< \
# 		$@

# ###################
# # Contig Assembly #
# ###################

# CONTIGS_R1 := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
# 	| sort \
# 	| uniq \
# 	| sed 's/$$/.fastq/' \
# 	| sed 's/^/data\/contigs\//')

# MOVE_CONTIGS := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
# 	| sort \
# 	| uniq \
# 	| sed 's/$$/.fastq/' \
# 	| sed 's/^/data\/contigfastq\//')

# assemblecontigs: $(CONTIGS_R1)

# movecontigs: $(MOVE_CONTIGS)

# $(CONTIGS_R1): data/contigs/%.fastq : data/HumanDecon/%_R1.fastq
# 	mkdir -p ./data/contigs
# 	bash ./bin/ContigAssembly.sh \
# 		$< \
# 		$(subst R1,R2,$<) \
# 		$@

# $(MOVE_CONTIGS): data/contigfastq/%.fastq :
# 	mkdir -p data/contigfastq
# 	cp $(subst contigfastq,contigs,$@)/final.contigs.fa $@

# # Merge the contigs into a master file
# ./data/totalcontigs.fa :
# 	cat ./data/contigfastq/* > $@
# 	sed -i 's/ /_/g' $@

# ###############################
# # Contig Abundance Per Sample #
# ###############################
# ./data/ContigRelAbundForGraph.tsv : \
# 			./data/totalcontigs.fa
# 	bash ./bin/CreateContigRelAbundTable.sh \
# 		./data/totalcontigs.fa \
# 		./data/HumanDecon \
# 		./data/ContigRelAbundForGraph.tsv

# #####################
# # Contig Clustering #
# #####################
# ./data/ContigAbundForConcoct.tsv : ./data/ContigRelAbundForGraph.tsv
# 	Rscript ./bin/ReshapeAlignedAbundance.R \
# 		-i ./data/ContigRelAbundForGraph.tsv \
# 		-o ./data/ContigAbundForConcoct.tsv \
# 		-p 0.5

# ./data/ContigClusters/clustering_gt1000.csv : \
# 			./data/totalcontigs.fa \
# 			./data/ContigAbundForConcoct.tsv
# 	mkdir ./data/ContigClustersPhage
# 	concoct \
# 		--coverage_file ./data/ContigAbundForConcoct.tsv \
# 		--composition_file ./data/totalcontigs.fa \
# 		--clusters 500 \
# 		--kmer_length 4 \
# 		--length_threshold 1000 \
# 		--read_length 150 \
# 		--basename ./data/ContigClusters/ \
# 		--no_total_coverage \
# 		--iterations 50

# ############################
# # Cluster Contig Abundance #
# ############################
# ./data/ClusteredContigAbund.tsv : \
# 			./data/ContigRelAbundForGraph.tsv \
# 			./data/ContigClusters/clustering_gt1000.csv
# 	bash ./bin/ClusterContigAbund.sh \
# 		./data/ContigRelAbundForGraph.tsv \
# 		./data/ContigClusters/clustering_gt1000.csv \
# 		./data/ClusteredContigAbund.tsv

# #################
# # Identify OPFs #
# #################
# ./data/totalopfs.fa : ./data/totalcontigs.fa
# 	bash ./bin/ClusterOPFs.sh \
# 		$< \
# 		$@

# ############################
# # ORF Abundance Per Sample #
# ############################
# orfalign:
# 	bash ./bin/getOrfAbundance.sh \
# 		./data/tmp-opfs/ContigOrfsNoSpec.fa \
# 		./data/HumanDecon

# ./data/orfabund.tsv:
# 	bash ./bin/catOrfAbundance.sh \
# 		./data/HumanDecon \
# 		./data/orfabund.tsv

# #######################
# # OPF Abundance Table #
# #######################
# ./data/ClusteredOpfAbund.tsv : \
# 			./data/totalopfs.fa.tsv \
# 			./data/orfabund.tsv
# 	$(shell awk '{ print $$2"\t"$$1 }' ./data/totalopfs.fa.tsv > ./data/OpfClusters.tsv)
# 	$(shell awk '{ print $$2"\t"$$1"\t"$$3 }' ./data/orfabund.tsv > ./data/orfabundOrdered.tsv)
# 	bash ./bin/ClusterContigAbund.sh \
# 		./data/orfabundOrdered.tsv \
# 		./data/OpfClusters.tsv \
# 		./data/ClusteredOpfAbund.tsv

# #######################
# # OGU Alpha Diversity #
# #######################
# ./figures/diversity-alpha-ogu.pdf : \
# 			./data/ClusteredContigAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-alpha.R \
# 		--input ./data/ClusteredContigAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--out ./figures/diversity-alpha-ogu.pdf

# #######################
# # OPF Alpha Diversity #
# #######################
# ./figures/diversity-alpha-opf.pdf : \
# 			./data/ClusteredOpfAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-alpha.R \
# 		--input ./data/ClusteredOpfAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--out ./figures/diversity-alpha-opf.pdf

# #######################
# # OGU Beta Diversity #
# #######################
# ./figures/diversity-beta-ogu.pdf : \
# 			./data/ClusteredContigAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredContigAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--out ./figures/diversity-beta-ogu.pdf

# ./figures/diversity-betajaccard-ogu.pdf : \
# 			./data/ClusteredContigAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredContigAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 100000 \
# 		--divmetric jaccard \
# 		--out ./figures/diversity-betajaccard-ogu.pdf

# #######################
# # OPF Beta Diversity #
# #######################
# ./figures/diversity-beta-opf.pdf : \
# 			./data/ClusteredOpfAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredOpfAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 50000 \
# 		--out ./figures/diversity-beta-opf.pdf

# ./figures/diversity-betajaccard-opf.pdf : \
# 			./data/ClusteredOpfAbund.tsv \
# 			./data/metadata/NexteraXT003Map.tsv
# 	Rscript ./bin/diversity-beta.R \
# 		--input ./data/ClusteredOpfAbund.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--subsample 50000 \
# 		--divmetric jaccard \
# 		--out ./figures/diversity-betajaccard-opf.pdf

# ####################
# # Sequencing Depth #
# ####################
# ./data/ProjectSeqDepth.tsv :
# 	bash ./bin/getPairedCount.sh \
# 		./data/HumanDecon \
# 		./data/ProjectSeqDepth.tsv

# ./figures/qualitycontrol.pdf :
# 	Rscript ./bin/QualityControlStats.R \
# 		--input ./data/ProjectSeqDepth.tsv \
# 		--metadata ./data/metadata/NexteraXT003Map.tsv \
# 		--out ./figures/qualitycontrol.pdf \
# 		--sdepth 100000

# ############################################### 16S ###############################################

# #################
# # Bacterial 16S #
# #################
# # Download the 16S reads from the SRA
# mothurproc :
# 	mkdir -p ./data/mothur16S
# 	bash ./bin/Mothur16S.sh \
# 		./data/mothur16S \
# 		./data/raw/Zackular_16S

# precluster :
# 	bash ./bin/mothurPreCluster.sh \
# 		./data/mothur16S

# mothurcluster :
# 	bash ./bin/mothurCluster.sh \
# 		./data/mothur16S

# ###############################
# # 16S Classification Modeling #
# ###############################


