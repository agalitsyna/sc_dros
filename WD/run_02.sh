for data_type in DATA DATA_SHUF DATA_SUBSAMPLE 
do
  mkdir -p ../${data_type}/PAIRIX
  mkdir -p ../${data_type}/COOL
  mkdir -p ../${data_type}/TXT/MTX/
  mkdir -p ../${data_type}/TXT/SPARSE/
done

exps="A6 B31 A8 A5 sc23 B26 B19 B16 sc16 B3 A2 A3 sc1 B6 A9 sc24 sc29 sc19 sc21 B15"
for pref in $exps
do 
  echo "$pref"
	sem -j 10 --id pairsam "python 02_pairsam2cooler.py ../DATA/PAIR/${pref}_\*.pairsam $pref"
  #python 02_pairsam2cooler_saving_marginals.py "\"../DATA/PAIR/${pref}_*.pairsam\"" "$pref"
done
sem --wait --id pairsam
