#! /bin/bash
# RenameFiles.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

for file in $(ls ./data/raw/NexteraXT003/*); do
	filename=$(echo ${file} | sed 's/_*//')
	echo File is ${file}
	echo File name is ${filename}
	# Get the alternate name to use for replacing
	altname=$(awk -v name=${file} 'name = $13 {print $2)')
	echo Alternate file name is ${altname}
	# Now actually replace the file name
	newname=$(echo ${file} | sed "s/$filename/$altname/" | sed 's/_[ATGC].*R/R/' | sed 's/_001//')
	echo New file is ${newname}
	# mv ${file} ${newname}
done
