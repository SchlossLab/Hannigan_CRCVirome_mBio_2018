#! /bin/bash
# downloadTaxRef.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

mkdir -p ./tmp-database-download

wget "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
mv nucl_gb.accession2taxid.gz ./tmp-database-download/
gunzip ./tmp-database-download/nucl_gb.accession2taxid.gz


