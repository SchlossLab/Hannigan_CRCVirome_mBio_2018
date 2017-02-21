# Makefile
# Hannigan-2016-ColonCancerVirome
# Geoffrey Hannigan
# Pat Schloss Lab

############################################# METADATA ############################################
metadatafiles = ./data/metadata/NexteraXT003Map.tsv ./data/metadata/NexteraXT004Map.tsv

# Make a master metadata file
./data/metadata/MasterMeta.tsv : $(metadatafiles)
	cat ./data/metadata/NexteraXT003Map.tsv ./data/metadata/NexteraXT004Map.tsv > ./data/metadata/MasterMeta.tsv

########################################## VIRUS DATABASE #########################################
./data/metadata/virus.txt:
	wget http://www.ebi.ac.uk/genomes/virus.txt -O $@

./data/metadata/phage.txt:
	wget http://www.ebi.ac.uk/genomes/phage.txt -O $@

# Merge the reference databases
./data/metadata/VirusPhageReference.tsv: ./data/metadata/virus.txt ./data/metadata/phage.txt
	cat ./data/metadata/virus.txt ./data/metadata/phage.txt > ./data/metadata/VirusPhageReference.tsv

# Download Viral Sequences
./data/metadata/VirusPhageReference.fa : ./data/metadata/VirusPhageReference.tsv
	bash ./bin/DownloadVirusesFromENA.sh \
		$< \
		$@

# Download Bacterial Sequences
./data/metadata/BacteriaReference.fa : ./data/metadata/bacteria.txt
	bash ./bin/DownloadVirusesFromENA.sh \
		$< \
		$@

 # Format reference fasta
./data/metadata/VirusPhageReferenceFormat.fa : ./data/metadata/VirusPhageReference.fa
	perl ./bin/remove_block_fasta_format.pl $< ./data/tmpreference.fa
	egrep -A 1 "complete genome" ./data/tmpreference.fa | egrep -v "\-\-" > $@
	sed -i 's/ /_/g' $@
	rm ./data/tmpreference.fa

# Note that the integrase gene had to be downloaded from the uniprot database manually.
# I used the same search query as my recent PeerJ manuscript.

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

####################
# Sequencing Depth #
####################
./data/ProjectSeqDepth.tsv :
	bash ./bin/getPairedCount.sh \
		./data/HumanDecon \
		./data/ProjectSeqDepth.tsv

./figures/qualitycontrol.pdf :
	Rscript ./bin/QualityControlStats.R \
		--input ./data/ProjectSeqDepth.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--out ./figures/qualitycontrol.pdf \
		--sdepth 100000

############################################# CONTIGS #############################################

###################
# Contig Assembly #
###################
setfile5: ./data/metadata/MasterMeta.tsv
	$(eval CONTIGS_R1 = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
		| sort \
		| uniq \
		| sed 's/$$/.fastq/' \
		| sed 's/^/data\/contigs\//'))
	echo 'variable5 = $(CONTIGS_R1)' > $@

setfile6: ./data/metadata/MasterMeta.tsv
	$(eval MOVE_CONTIGS = $(shell awk '{ print $$2 }' ./data/metadata/MasterMeta.tsv \
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

# Get a master file for bacterial and viral contigs
setfile7: ./data/metadata/MasterMeta.tsv
	$(eval VIRUS_CONTIGS := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
		| sort \
		| uniq \
		| sed 's/$$/.fastq/' \
		| sed 's/^/data\/makerecorddump\//'))
	echo 'variable7 = $(VIRUS_CONTIGS)' > $@

setfile8: ./data/metadata/MasterMeta.tsv
	$(eval BACTERIA_CONTIGS := $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT004Map.tsv \
		| sort \
		| uniq \
		| sed 's/$$/.fastq/' \
		| sed 's/^/data\/makerecorddump\//'))
	echo 'variable8 = $(BACTERIA_CONTIGS)' > $@

contigpairs:

include setfile7
include setfile8

cleancontigpairs:
	rm -f ./data/totalcontigsvirus.fa
	rm -f ./data/totalcontigsbacteria.fa
	rm -rf data/makerecorddump

contigpairs: $(variable7) $(variable8)

$(variable7): data/makerecorddump/%.fastq : data/contigfastq/%.fastq
	mkdir -p data/makerecorddump
	touch $@
	cat $< >> ./data/totalcontigsvirus.fa
	sed -i 's/ /_/g' ./data/totalcontigsvirus.fa

$(variable8): data/makerecorddump/%.fastq : data/contigfastq/%.fastq
	mkdir -p data/makerecorddump
	touch $@
	cat $< >> ./data/totalcontigsbacteria.fa
	sed -i 's/ /_/g' ./data/totalcontigsbacteria.fa

######################################### CONTIG ABUNDANCE ########################################

###############################
# Contig Abundance Per Sample #
###############################
# Viruses
setfile9: ./data/metadata/MasterMeta.tsv
	$(eval VIRUS_MOVE = $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R2.fastq/' \
		| sed 's/^/data\/virusseqsfastq\//'))
	echo 'variable9 = $(VIRUS_MOVE)' > $@

movevirusabund:

include setfile9

movevirusabund: $(variable9)

$(variable9): data/virusseqsfastq/%_R2.fastq : data/HumanDecon/%_R2.fastq
	mkdir -p data/virusseqsfastq
	cp $< $@

setfile9_1: ./data/metadata/MasterMeta.tsv
	$(eval BUILD_VIRUS_ABUND = $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT003Map.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R2.fastq-noheader-forcat/' \
		| sed 's/^/data\/virusseqsfastq\//'))
	echo 'variable9_1 = $(BUILD_VIRUS_ABUND)' > $@

virusabund:

include setfile9_1

virusabund: $(variable9_1)

./data/virusbowtieReference/bowtieReference.1.bt2 : ./data/totalcontigsvirus.fa
	mkdir -p ./data/virusbowtieReference
	bowtie2-build \
		-q ./data/totalcontigsvirus.fa \
		./data/virusbowtieReference/bowtieReference

$(variable9_1): data/virusseqsfastq/%_R2.fastq-noheader-forcat : data/virusseqsfastq/%_R2.fastq ./data/virusbowtieReference/bowtieReference.1.bt2
	qsub ./bin/CreateContigRelAbundTable.pbs -F './data/virusbowtieReference/bowtieReference $<'

# Bacteria
setfile10: ./data/metadata/MasterMeta.tsv
	$(eval BACTERIA_MOVE = $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT004Map.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R2.fastq/' \
		| sed 's/^/data\/bacteriaseqsfastq\//'))
	echo 'variable10 = $(BACTERIA_MOVE)' > $@

movebacteriaabund:

include setfile10

movebacteriaabund: $(variable10)

$(variable10): data/bacteriaseqsfastq/%_R2.fastq : data/HumanDecon/%_R2.fastq
	mkdir -p data/bacteriaseqsfastq
	cp $< $@

setfile10_1: ./data/metadata/MasterMeta.tsv
	$(eval BUILD_BACTERIA_ABUND = $(shell awk '{ print $$2 }' ./data/metadata/NexteraXT004Map.tsv \
		| sort \
		| uniq \
		| sed 's/$$/_R2.fastq-noheader-forcat/' \
		| sed 's/^/data\/bacteriaseqsfastq\//'))
	echo 'variable10_1 = $(BUILD_BACTERIA_ABUND)' > $@

bacteriaabund:

include setfile10_1

bacteriaabund: movebacteriaabund $(variable10_1)

./data/bacteriabowtieReference/bowtieReference.1.bt2 : ./data/totalcontigsbacteria.fa
	mkdir -p ./data/bacteriabowtieReference
	bowtie2-build \
		-q ./data/totalcontigsbacteria.fa \
		./data/bacteriabowtieReference/bowtieReference

$(variable10_1): data/bacteriaseqsfastq/%_R2.fastq-noheader-forcat : data/bacteriaseqsfastq/%_R2.fastq ./data/bacteriabowtieReference/bowtieReference.1.bt2
	qsub ./bin/CreateContigRelAbundTable.pbs -F './data/bacteriabowtieReference/bowtieReference $<'

# The results from each of these rule sets (bacteria and virus) need to be combined
./data/ContigRelAbundForGraphVirus.tsv : virusabund
	cat ./data/virusseqsfastq/*_R2.fastq-noheader-forcat > $@
	sed -i 's/_R2.fastq//g' $@
	sed -i 's/data\/virusseqsfastq\///g' $@

./data/ContigRelAbundForGraphBacteria.tsv : bacteriaabund
	cat ./data/bacteriaseqsfastq/*_R2.fastq-noheader-forcat > $@
	sed -i 's/_R2.fastq//g' $@
	sed -i 's/data\/bacteriaseqsfastq\///g' $@
	

######################################## CONTIG CLUSTERING ########################################

#####################
# Contig Clustering #
#####################
# Bacteria
# Lower sample number because of memory limitations
./data/ContigAbundForConcoctBacteria.tsv : ./data/ContigRelAbundForGraphBacteria.tsv
	Rscript ./bin/ReshapeAlignedAbundance.R \
		-i ./data/ContigRelAbundForGraphBacteria.tsv \
		-o ./data/ContigAbundForConcoctBacteria.tsv \
		-p 0.25

./data/ContigClustersBacteria/clustering_gt2000.csv : \
			./data/totalcontigsbacteria.fa \
			./data/ContigAbundForConcoctBacteria.tsv
	mkdir ./data/ContigClustersBacteria
	concoct \
		--coverage_file ./data/ContigAbundForConcoctBacteria.tsv \
		--composition_file ./data/totalcontigsbacteria.fa \
		--clusters 500 \
		--kmer_length 4 \
		--length_threshold 2000 \
		--read_length 150 \
		--basename ./data/ContigClustersBacteria/ \
		--no_total_coverage \
		--iterations 50

# Virus
./data/ContigAbundForConcoctVirus.tsv : ./data/ContigRelAbundForGraphVirus.tsv
	Rscript ./bin/ReshapeAlignedAbundance.R \
		-i ./data/ContigRelAbundForGraphVirus.tsv \
		-o ./data/ContigAbundForConcoctVirus.tsv \
		-p 0.5

./data/ContigClustersVirus/clustering_gt1000.csv : \
			./data/totalcontigsvirus.fa \
			./data/ContigAbundForConcoctVirus.tsv
	mkdir ./data/ContigClustersVirus
	concoct \
		--coverage_file ./data/ContigAbundForConcoctVirus.tsv \
		--composition_file ./data/totalcontigsvirus.fa \
		--clusters 500 \
		--kmer_length 4 \
		--length_threshold 1000 \
		--read_length 150 \
		--basename ./data/ContigClustersVirus/ \
		--no_total_coverage \
		--iterations 50

####################################### CLUSTERED ABUNDANCE #######################################

############################
# Cluster Contig Abundance #
############################
# Virus
# Get table of contig IDs and their lengths
./data/VirusContigLength.tsv : \
			./data/totalcontigsvirus.fa
	perl ./bin/ContigLengthTable.pl \
		-i ./data/totalcontigsvirus.fa \
		-o ./data/VirusContigLength.tsv

./data/CorrectedContigRelAbundForGraph.tsv : \
			./data/VirusContigLength.tsv
	perl ./bin/AbundLengthCorrection.pl \
		-i ./data/ContigRelAbundForGraphVirus.tsv \
		-l ./data/VirusContigLength.tsv \
		-o ./data/CorrectedContigRelAbundForGraphVirus.tsv \
		-f 1000

./data/VirusClusteredContigAbund.tsv : \
			./data/CorrectedContigRelAbundForGraphVirus.tsv
	bash ./bin/ClusterContigAbund.sh \
		./data/CorrectedContigRelAbundForGraphVirus.tsv \
		./data/ContigClustersVirus/clustering_gt1000.csv \
		./data/VirusClusteredContigAbund.tsv

# Bacteria
# Get table of contig IDs and their lengths
./data/BacteriaContigLength.tsv : \
			./data/totalcontigsbacteria.fa
	perl ./bin/ContigLengthTable.pl \
		-i ./data/totalcontigsbacteria.fa \
		-o ./data/BacteriaContigLength.tsv

./data/CorrectedContigRelAbundForGraphBacteria.tsv : \
			./data/BacteriaContigLength.tsv
	perl ./bin/AbundLengthCorrection.pl \
		-i ./data/ContigRelAbundForGraphBacteria.tsv \
		-l ./data/BacteriaContigLength.tsv \
		-o ./data/CorrectedContigRelAbundForGraphBacteria.tsv \
		-f 2000

./data/BacteriaClusteredContigAbund.tsv : \
			./data/CorrectedContigRelAbundForGraphBacteria.tsv
	bash ./bin/ClusterContigAbund.sh \
		./data/CorrectedContigRelAbundForGraphBacteria.tsv \
		./data/ContigClustersBacteria/clustering_gt2000.csv \
		./data/BacteriaClusteredContigAbund.tsv

################################### OPERATIONAL PROTEIN FAMILIES ##################################

#################
# Identify OPFs #
#################
idorfs: ./data/totalcontigorfsvirus.fa ./data/totalcontigorfsbacteria.fa

# Virus
./data/totalcontigorfsvirus.fa : ./data/totalcontigsvirus.fa
	bash ./bin/ClusterOPFs.sh \
		$< \
		$@

# Bacteria
./data/totalcontigorfsbacteria.fa : ./data/totalcontigsbacteria.fa
	bash ./bin/ClusterOPFs.sh \
		$< \
		$@

############################
# ORF Abundance Per Sample #
############################
orfalign:
	bash ./bin/getOrfAbundance.sh \
		./data/tmp-opfs/ContigOrfsNoSpec.fa \
		./data/HumanDecon

./data/orfabund.tsv:
	bash ./bin/catOrfAbundance.sh \
		./data/HumanDecon \
		./data/orfabund.tsv

#######################
# OPF Abundance Table #
#######################
./data/ClusteredOpfAbund.tsv : \
			./data/totalopfs.fa.tsv \
			./data/orfabund.tsv
	$(shell awk '{ print $$2"\t"$$1 }' ./data/totalopfs.fa.tsv > ./data/OpfClusters.tsv)
	$(shell awk '{ print $$2"\t"$$1"\t"$$3 }' ./data/orfabund.tsv > ./data/orfabundOrdered.tsv)
	bash ./bin/ClusterContigAbund.sh \
		./data/orfabundOrdered.tsv \
		./data/OpfClusters.tsv \
		./data/ClusteredOpfAbund.tsv

############################################ DIVERSITY ############################################

#######################
# OGU Alpha Diversity #
#######################
./figures/diversity-alpha-ogu.pdf : \
			./data/ClusteredContigAbund.tsv \
			./data/metadata/NexteraXT003Map.tsv
	Rscript ./bin/diversity-alpha.R \
		--input ./data/ClusteredContigAbund.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--subsample 100000 \
		--out ./figures/diversity-alpha-ogu.pdf

#######################
# OPF Alpha Diversity #
#######################
./figures/diversity-alpha-opf.pdf : \
			./data/ClusteredOpfAbund.tsv \
			./data/metadata/NexteraXT003Map.tsv
	Rscript ./bin/diversity-alpha.R \
		--input ./data/ClusteredOpfAbund.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--subsample 100000 \
		--out ./figures/diversity-alpha-opf.pdf

#######################
# OGU Beta Diversity #
#######################
./figures/diversity-beta-ogu.pdf \
./figures/diversity-beta-ogu-negative.pdf : \
			./data/ClusteredContigAbund.tsv \
			./data/metadata/NexteraXT003Map.tsv
	Rscript ./bin/diversity-beta.R \
		--input ./data/ClusteredContigAbund.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--subsample 100000 \
		--out ./figures/diversity-beta-ogu.pdf \
		--negout ./figures/diversity-beta-ogu-negative.pdf

./figures/diversity-betajaccard-ogu.pdf \
./figures/diversity-betajaccard-ogu-negative.pdf : \
			./data/ClusteredContigAbund.tsv \
			./data/metadata/NexteraXT003Map.tsv
	Rscript ./bin/diversity-beta.R \
		--input ./data/ClusteredContigAbund.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--subsample 100000 \
		--divmetric jaccard \
		--out ./figures/diversity-betajaccard-ogu.pdf \
		--negout ./figures/diversity-betajaccard-ogu-negative.pdf

#######################
# OPF Beta Diversity #
#######################
./figures/diversity-beta-opf.pdf : \
			./data/ClusteredOpfAbund.tsv \
			./data/metadata/NexteraXT003Map.tsv
	Rscript ./bin/diversity-beta.R \
		--input ./data/ClusteredOpfAbund.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--subsample 50000 \
		--out ./figures/diversity-beta-opf.pdf

./figures/diversity-betajaccard-opf.pdf : \
			./data/ClusteredOpfAbund.tsv \
			./data/metadata/NexteraXT003Map.tsv
	Rscript ./bin/diversity-beta.R \
		--input ./data/ClusteredOpfAbund.tsv \
		--metadata ./data/metadata/NexteraXT003Map.tsv \
		--subsample 50000 \
		--divmetric jaccard \
		--out ./figures/diversity-betajaccard-opf.pdf

############################################### 16S ###############################################

#################
# Bacterial 16S #
#################
# Download the 16S reads from the SRA
mothurproc :
	mkdir -p ./data/mothur16S
	bash ./bin/Mothur16S.sh \
		./data/mothur16S \
		./data/raw/Zackular_16S

precluster :
	bash ./bin/mothurPreCluster.sh \
		./data/mothur16S

mothurcluster :
	bash ./bin/mothurCluster.sh \
		./data/mothur16S

####################################### INTERACTION NETWORKS ######################################

##################################
# Score Contigs for Interactions #
##################################

VREF=./data/ViromeAgainstReferenceBacteria
# In this case the samples will get run against the bacteria reference genome set
ViromeRefRun : ${VREF}/BenchmarkCrisprsFormat.tsv \
	${VREF}/BenchmarkProphagesFormatFlip.tsv \
	${VREF}/MatchesByBlastxFormatOrder.tsv \
	${VREF}/PfamInteractionsFormatScoredFlip.tsv

${VREF}/BenchmarkCrisprsFormat.tsv \
${VREF}/BenchmarkProphagesFormatFlip.tsv \
${VREF}/MatchesByBlastxFormatOrder.tsv \
${VREF}/PfamInteractionsFormatScoredFlip.tsv : \
			./data/totalcontigsvirus.fa \
			./data/totalcontigsbacteria.fa \
			./bin/BenchmarkingModel.sh
	bash ./bin/BenchmarkingModel.sh \
		./data/totalcontigsvirus.fa \
		./data/totalcontigsbacteria.fa \
		${VREF}/BenchmarkCrisprsFormat.tsv \
		${VREF}/BenchmarkProphagesFormatFlip.tsv \
		${VREF}/MatchesByBlastxFormatOrder.tsv \
		${VREF}/PfamInteractionsFormatScoredFlip.tsv \
		"ViromeAgainstReferenceBacteria"

#########################
# Cluster Contig Scores #
#########################

# Annotate contig IDs with cluster IDs and further compress
clusterrun : ${VREF}/BenchmarkProphagesFormatFlipClustered.tsv \
	${VREF}/MatchesByBlastxFormatOrderClustered.tsv \
	${VREF}/PfamInteractionsFormatScoredFlipClustered.tsv

${VREF}/BenchmarkProphagesFormatFlipClustered.tsv \
${VREF}/MatchesByBlastxFormatOrderClustered.tsv \
${VREF}/PfamInteractionsFormatScoredFlipClustered.tsv :
	bash ./bin/ClusterContigScores.sh \
		${VREF}/BenchmarkProphagesFormatFlip.tsv \
		${VREF}/MatchesByBlastxFormatOrder.tsv \
		${VREF}/PfamInteractionsFormatScoredFlip.tsv \
		${VREF}/BenchmarkProphagesFormatFlipClustered.tsv \
		${VREF}/MatchesByBlastxFormatOrderClustered.tsv \
		${VREF}/PfamInteractionsFormatScoredFlipClustered.tsv \
		./data/ContigClustersVirus/clustering_gt1000.csv \
		./data/ContigClustersBacteria/clustering_gt2000.csv \
		"ViromeAgainstReferenceBacteria" \
		${VREF}/BenchmarkCrisprsFormat.tsv \
		${VREF}/BenchmarkCrisprsFormatClustered.tsv

# Make a graph database from the experimental information
expnetwork :
	# Note that this resets the graph database and erases
	# the validation information we previously added.
	rm -r ../../bin/neo4j-enterprise-2.3.0/data/graph.db/
	mkdir ../../bin/neo4j-enterprise-2.3.0/data/graph.db/
	bash ./bin/CreateProteinNetwork \
		./data/ValidationSet/Interactions.tsv \
		${VREF}/BenchmarkCrisprsFormatClustered.tsv \
		${VREF}/BenchmarkProphagesFormatFlipClustered.tsv \
		${VREF}/PfamInteractionsFormatScoredFlipClustered.tsv \
		${VREF}/MatchesByBlastxFormatOrderClustered.tsv \
		"FALSE"

# Predict interactions between nodes
./data/PredictedRelationshipTable.tsv :
	bash ./bin/RunPredictionsWithNeo4j.sh \
		./data/metadata/rfinteractionmodel.RData \
		./data/PredictedRelationshipTable.tsv

# Add relationships
finalrelationships \
./figures/BacteriaPhageNetworkDiagram.pdf \
./figures/BacteriaPhageNetworkDiagram.png \
./figures/PhageHostHist.pdf \
./figures/PhageHostHist.png \
./figures/BacteriaEdgeCount.pdf \
./figures/BacteriaEdgeCount.png : \
		./data/PredictedRelationshipTable.tsv \
		./bin/AddRelationshipsWrapper.sh
	echo $(shell date)  :  Adding relationships to network and plotting total graph >> ${DATENAME}.makelog
	bash ./bin/AddRelationshipsWrapper.sh \
		./data/PredictedRelationshipTable.tsv

################################## CONTIG CLUSTER IDENTIFICATION ##################################
# Get ID for longest contig in each cluster
./data/contigclustersidentity/longestcontigsvirus.tsv : ./data/VirusContigLength.tsv
	mkdir -p ./data/contigclustersidentity
	Rscript ./bin/GetLongestContig.R \
		--input ./data/VirusContigLength.tsv \
		--clusters ./data/ContigClustersVirus/clustering_gt1000.csv \
		--toplength 1 \
		--out $@

# Align the contig seqs to the virus reference database
./data/contigclustersidentity/VirusRepsetIds.tsv :
	bash ./bin/IdentifyContigs.sh \
		./data/totalcontigsvirus.fa \
		./data/metadata/VirusPhageReferenceFormat.fa \
		./data/contigclustersidentity/longestcontigsvirus.tsv \
		./data/contigclustersidentity/VirusRepsetIds.tsv \
		"/home/ghannig/bin/ncbi-blast-2.4.0+/bin/"

./rtables/idcount.tsv :
	bc <<< "scale=8; 100 * `wc -l < ./data/contigclustersidentity/clustax.tsv` / `wc -l < ./data/contigclustersidentity/longestcontigsvirus.tsv`" > ./rtables/idcount.tsv

./rtables/phagecount.tsv :
	bc <<< "scale=8; 100 * `egrep "Caudovirales|[Pp]hage" ./data/contigclustersidentity/clustax.tsv | wc -l` / `wc -l < ./data/contigclustersidentity/clustax.tsv`" > ./rtables/phagecount.tsv

################################# BACTERIA CLUSTER IDENTIFICATION #################################
# Get ID for longest contig in each cluster
./data/contigclustersidentity/longestcontigsbacteria.tsv : ./data/BacteriaContigLength.tsv
	mkdir -p ./data/contigclustersidentity
	Rscript ./bin/GetLongestContig.R \
		--input ./data/BacteriaContigLength.tsv \
		--clusters ./data/ContigClustersBacteria/clustering_gt2000.csv \
		--toplength 1 \
		--out $@

# Align the contig seqs to the virus reference database
./data/contigclustersidentity/BacteriaRepsetIds.tsv :
	bash ./bin/IdentifyContigsBac.sh \
		./data/totalcontigsbacteria.fa \
		./data/metadata/BacteriaReference.fa \
		./data/contigclustersidentity/longestcontigsbacteria.tsv \
		./data/contigclustersidentity/BacteriaRepsetIds.tsv \
		"/home/ghannig/bin/ncbi-blast-2.4.0+/bin/"

################################# PHAGE REPLICATION CYCLE #################################
./data/tmpidlytic/phagebacteriablastout.tsv :
	bash ./bin/IdentLyticPhages.sh \
		./data/totalcontigsvirus.fa \
		./data/contigclustersidentity/longestcontigsvirus.tsv \
		"/home/ghannig/bin/ncbi-blast-2.4.0+/bin/"







