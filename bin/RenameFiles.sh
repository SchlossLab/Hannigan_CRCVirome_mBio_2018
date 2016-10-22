#! /bin/bash
# RenameFiles.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

export metadata=$1

for file in $(ls ./data/raw/NexteraXT003/*); do
	filename=$(echo ${file} | sed 's/.*\///g' | sed 's/_.*//')
	echo File is ${file}
	echo File name is ${filename}
	# Get the alternate name to use for replacing
	altname=$(awk -v name=${file} '$13 == name {print $2}' ${metadata})
	echo Alternate file name is ${altname}
	# Now actually replace the file name
	newname=$(echo ${file} | sed "s/$filename/$altname/" | sed 's/_[ATGC].*R/R/' | sed 's/_001//')
	echo New file is ${newname}
	# mv ${file} ${newname}
done
