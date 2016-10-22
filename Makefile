# Makefile
# Hannigan-2016-ColonCancerVirome
# Geoffrey Hannigan
# Pat Schloss Lab

#PBS -N ContigAssemblyAndAlignment
#PBS -q first
#PBS -l nodes=1:ppn=1,mem=40gb
#PBS -l walltime=500:00:00
#PBS -j oe
#PBS -V
#PBS -A schloss_lab

# Run quality control processes
# I know it's hardcoded now but the refactoring can wait for another time
./QualitySeqs/SequenceCounts/16sHits.tsv ./QualitySeqs/DeconSeq ./figures:./data/NexteraXT002Map.tsv ./data/raw/NexteraXT002
	bash ./bin/QualityProcess.sh

./data/ContigsAndAlignments/ContigRelativeAbundanceTable.tsv:./data/QualitySeqs/DeconSeq ./data/NexteraXT002Map.tsv
	bash ./bin/ContigAssemblyAndAlignment.sh
