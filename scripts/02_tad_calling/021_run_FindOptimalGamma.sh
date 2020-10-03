#!/bin/bash

for pref in Cell1 Cell2 Cell3 Cell4 Cell5 Cell6 Cell7 Cell8 Cell9 Cell10 Cell11 Cell12 Cell13 Cell14 Cell15 Cell16 Cell17 Cell18 Cell19 Cell20
do
  TAD_PATH="../../data/TAD/"
  IMG_PATH="../../data/IMG/"
  BED_PATH="../../data/TAD_modularity/"
  python 021_FindOptimalGamma.py ${TAD_PATH}/${pref}.10.mod-3.segmentation_all.pickle ${pref} minTAD3_3 ${IMG_PATH} ${BED_PATH}
done
