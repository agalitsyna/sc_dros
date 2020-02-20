#!/bin/bash

GENOME="dm3" 
GENOMEDIR="./data/GENOME"
INDIR="./data/FASTQ"

BAMDIR="./data/BAM"
PAIRDIR="./data/PAIR"
COOLDIR="./data/COOL"
STATSDIR="./data/STATS"

mkdir -p $BAMDIR
mkdir -p $PAIRDIR
mkdir -p $COOLDIR
mkdir -p $STATSDIR

# Listing all input files
FILELIST=(`ls $INDIR/*.fastq.gz`)
NFILES=$(( ${#FILELIST[@]}/2 ))

echo "Processing fastq.gz files, number of paired-end files: $NFILES"

for (( i=0; i<$NFILES; i++ ))
do
  
  j1=$(($i*2))
  
  dir="${FILELIST[$j1]}"
  pref=${dir##*/}
  pref=${pref%_*}

  file1="${pref}_R1.fastq.gz"
  file2="${pref}_R2.fastq.gz"
  
  echo $pref $file1 $file2

  # Mapping with bwa-mem
  bwa mem -t 15 -v 3 -SP ${GENOMEDIR}/${GENOME}.fa ${INDIR}/$file2 ${INDIR}/$file1 | samtools view -bS > ${BAMDIR}/${pref}.bam

  # Parsing with pairtools ORBITA (with rfrag annotation and ligation juction filtering)
  pairtools parse ${BAMDIR}/${pref}.bam -c $GENOMEDIR/${GENOME}.reduced.chrom.sizes --output ${PAIRDIR}/${pref}.pair -f $GENOMEDIR/rfrags_${GENOME}_DpnII.txt --walks-policy ligation_junctions --drop-sam --add-columns rfrag --output-stats ${STATSDIR}/${pref}.txt

  # Population Hi-C data
  if [[ $pref =~ 'Population' ]]
  then

    # Straightforward parsing of population data, analogous to distiller, see https://github.com/mirnylab/distiller-nf
    pairtools parse ${BAMDIR}/${pref}.bam -c $GENOMEDIR/${GENOME}.reduced.chrom.sizes --output ${PAIRDIR}/${pref}.pair.full -f $GENOMEDIR/rfrags_${GENOME}_DpnII.txt --walks-policy mask --drop-sam --output-stats ${STATSDIR}/${pref}.txt.full
    pairtools sort --nproc 15 -o ${PAIRDIR}/${pref}.pair.full.sorted ${PAIRDIR}/${pref}.pair.full

    ## Restriction fragments annotation for the calculation of statistics, for debug only
    #pairtools restrict ${PAIRDIR}/${pref}.pair.full --frags $GENOMEDIR/rfrags_${GENOME}_DpnII.txt --output ${PAIRDIR}/${pref}.pair.full.restricted
  fi

done

