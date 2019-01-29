#!/bin/bash
#PBS -l walltime=300:00:00,mem=1gb
#PBS -t 0-999
#PBS -d.

TAD_PATH="../DATA/TAD/"
IMG_PATH="../DATA/IMG/"
BED_PATH="../DATA/TAD_modularity/"

source ~/activate_env.sh

ARRAY=()
for pref in A8 A9 B15 B19 B26 B3 B31 sc19 sc29 A2 A3 A5 B16 B6 sc1 sc16 sc21 sc23 sc24 A6
do
  sem -j 20 python 05_find_optimal.py ${TAD_PATH}/${pref}.10.cool.subsample.${frac}.${num}.mod-12-4.segmentation_all.pickle ${pref} subsample_${frac}_${num}_minTAD4 ${IMG_PATH} ${BED_PATH}
done

sem --wait
