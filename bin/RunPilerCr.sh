#! /bin/bash
# RunPilerCr.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

##############################
# Default Values and CL Args #
##############################
# Set pilerCR path
export PilerPath=$3
# Maximum input file size in bytes
export MaxFileSize=500000000 #500 MB
export FastaInput=$1
export OutputName=$2
export SplitSize=50
export Remove=TRUE
# Determine file size of input file
FileSize=$(wc -c "${FastaInput}" | sed 's/ .*//')
# Specify working directory
WorkDir=$(pwd)
echo "We are working in ${WorkDir}"

###################
# Extract CRISPRs #
###################
# Make a tmp123 directory
rm -r ./tmp123
mkdir ./tmp123
# Split files if needed
if [[ "${FileSize}" -gt "${MaxFileSize}" ]]; then
	echo "Input larger than ${MaxFileSize}."
	# Split the file
	split \
		--lines=${SplitSize} \
		"${FastaInput}" \
		./tmp123/tmp123Piler-
else
	echo "File is small so does not need split."
	# Copy file to tmp123 for ease
	cp "${FastaInput}" ./tmp123/
fi

# Now run pilerCR on the files
ls ./tmp123/* | xargs -I {} --max-procs=32 ${PilerPath}pilercr -quiet -in {} -out {}.out

# Collect the results together
cat ./tmp123/*.out > ./"$OutputName"

# Finally remove the tmp123 directories
if [[ "${Remove}" = "FALSE" ]]; then
	echo "Not removing tmp123 dir..."
else
	echo "Removing tmp123 dir..."
	rm -r ./tmp123
fi
