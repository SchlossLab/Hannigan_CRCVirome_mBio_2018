#! /bin/bash
# ProdigalWrapperLargeFiles.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan
# WARNING: This cannot take block formatted fasta files!
# WARNING: This contains a tmp directory removal step!

##############################
# Default Values and CL Args #
##############################
# Set pilerCR path
export ProdigalPath=/mnt/EXT/Schloss-data/bin
# Maximum input file size in bytes
export MaxFileSize=25000000 #25 MB
export FastaInput=$1
export OutputName=$2
export SplitSize=20
export Remove=TRUE
# Determine file size of input file
export FileSize
FileSize=$(wc -c "${FastaInput}" | sed 's/ .*//')
# Specify working directory
export WorkDir
WorkDir=$(pwd)
echo "We are working in ${WorkDir}"

#############
# Call ORFs #
#############
# Make a tmp directory after cleaning any existing tmp directories
rm -r ./tmp
mkdir ./tmp
# Split files if needed
if [[ "${FileSize}" -gt "${MaxFileSize}" ]]; then
	echo "Input larger than ${MaxFileSize} B."
	# Split the file
	split \
		--suffix-length=7 \
		--lines=${SplitSize} \
		"${FastaInput}" \
		./tmp/tmpProdigal-
else
	echo "File is small so does not need split."
	# Copy file to tmp for ease
	cp "${FastaInput}" ./tmp/
fi

# Now run pilerCR on the files
echo "Running Prodigal."
ls ./tmp/* | xargs -I {} --max-procs=8 ${ProdigalPath}/prodigal -q -c -i {} -o {}.genes -a {}.out -d {}.nucl -p meta

# Collect the results together
cat ./tmp/*.out > ./"${OutputName}"
cat ./tmp/*.nucl > ./"${OutputName}".nucleotide

# Finally remove the tmp directories
if [[ "${Remove}" = "FALSE" ]]; then
	echo "Keeping tmp dir..."
else
	echo "Removing tmp dir..."
	rm -r ./tmp
fi
