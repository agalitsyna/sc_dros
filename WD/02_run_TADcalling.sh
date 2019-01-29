sem -j 20 "python 02_TADcalling.sh ../DATA/COOL/Dros.10.cool.full ../DATA/TAD/population.10.tmp_minsize4 modularity balanced 12 4"

# Single-cell Hi-C Drosophila
COOL_PATH="../DATA/COOL/"
TAD_PATH="../DATA/TAD/"
for pref in A9 B15 B19 B26 B3 B31 sc19 sc29 A2 A3 A5 A6 A8 B16 B6 sc1 sc16 sc21 sc23 sc24
do
  sem -j 35 "python 02_TADcalling.sh ${COOL_PATH}/${pref}${mode}.cool ${TAD_PATH}/${pref}${mode}.10.mod-12-4 modularity non-balanced 12 4"
done

# Models datasets
COOL_PATH="../DATA_MODEL/COOL/"
TAD_PATH="../DATA_MODEL/TAD/"
for mode in "" "_sh"
  do
    for pref in A9 B15 B19 B26 B3 B31 sc19 sc29 A2 A3 A5 A6 A8 B16 B6 sc1 sc16 sc21 sc23 sc24
    do
      sem -j 35 "python 02_TADcalling.sh ${COOL_PATH}/${pref}${mode}.cool ${TAD_PATH}/${pref}${mode}.10.mod-12-4 modularity non-balanced 12 4"
  done
done

# Shuffled datasets
COOL_PATH="../DATA_SHUF/COOL/"
TAD_PATH="../DATA_SHUF/TAD/"
for num in $(seq 0 2) 
  do
    for pref in A6 A8 A9 B15 B19 B26 B3 B31 sc19 sc29 A2 A3 A5 B16 B6 sc1 sc16 sc21 sc23 sc24 
    do
      file="${TAD_PATH}/${pref}.10.random.${num}.mod-12-4.segmentation_all.pickle"
      if [ -f "$file" ]
      then
        echo $file
      else
        echo $file not found
        sem -j 35 "python 02_TADcalling.sh ${COOL_PATH}/${pref}.10.cool.random.${num} ${TAD_PATH}/${pref}.10.random.${num}.mod-12-4 modularity non-balanced 12 4"
      fi
  done
done

#SUBSAMPLE
COOL_PATH="../DATA_SUBSAMPLE/COOL/"
TAD_PATH="../DATA_SUBSAMPLE/TAD/"
for num in $(seq 0 2) 
  do
    for pref in A8 A9 B15 B19 B26 B3 B31 sc19 sc29 A2 A3 A5 B16 B6 sc1 sc16 sc21 sc23 sc24  #A6
    do
      for frac in 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100
      do
        sem -j 15 "python 02_TADcalling.sh ${COOL_PATH}/${pref}.10.cool.subsample.${frac}.${num} ${TAD_PATH}/${pref}.10.cool.subsample.${frac}.${num}.mod-12-4 modularity non-balanced 12 4"
      done
  done
done

# FLYAMER
COOL_PATH="../DATA_FL/COOL/"
TAD_PATH="../DATA_FL/TAD/"
for pref in oS181 oN65 #oS31 oS174 oS175 oSH124 oSH33 oS167 oNH188 oN152 oS150 oS20 oSH9 oN25 oN30 oSH6 oS173 oS13 oS219 oSH201 oS3 oNH193 oSH34 oS144 oSH166 oN15 oSH16 oS194 oN155 oN159 oS35 oIH57 oSH42 oNH131 oS83 oNH186 oS130 oN5 oN47 oSH19 oS170 oSH11 oN218 oN153 oS177 oS8 oS51 oSH22 oS37 oS163 oN62 oN44 oNH176 oS29 oS191 oSH21 oSH162 oS99 oNH86 oI7 oSH192 oS14 oS165 oNH146 oI43 oNH138 oSH28 oS101 oN10 oIH50 oIH84 oSH61 oNH179 oSH164 oS161 oS122 oN55 oNH171 oS132 oSH143 oN182 oS68 oIH120 oSH197 oSH23 oN168 oN1 oN41 oS82 oNH140 oS172 oNH26 oSH52 oSH18 oS125 oN40 oNH145 oN189 oN59 oN48 oNH93 oSH56 oSH157 oN4 oS27 oN32 oS46 oN17 oS114 oS185 oS39 oSH135 oN134 oNH160 oN12 oN180 oSH49 oSH106 oSH2 
do
  sem -j 5 "python 02_TADcalling.sh ${COOL_PATH}/${pref}.10.cool ${TAD_PATH}/${pref}.10.mod-12-4 modularity non-balanced 12 4"
done

sem --wait
