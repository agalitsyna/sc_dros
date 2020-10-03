cd ../../data/GENOME/

# Obtain file with genome from UCSC Golden Path
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz
tar xvzf chromFa.tar.gz
cat chr2L.fa chr2R.fa chr3L.fa chr3R.fa chr4.fa chrX.fa > dm3.fa

# Build index for bwa
bwa index dm3.fa
