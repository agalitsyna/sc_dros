for data_type in data data_subsampled data_shuffled 
do
  mkdir -p ./${data_type}/PAIRIX
  mkdir -p ./${data_type}/COOL
  mkdir -p ./${data_type}/TXT/MTX/
  mkdir -p ./${data_type}/TXT/SPARSE/
done

exps="Cell1 Cell2 Cell3 Cell4 Cell5 Cell6 Cell7 Cell8 Cell9 Cell10 Cell11 Cell12 Cell13 Cell14 Cell15 Cell16 Cell17 Cell18 Cell19 Cell20"
for pref in $exps
do 
  echo "$pref"
  
    # Extracting JJ and PP pairs for QC
    for file in ../DATA/PAIR/${pref}_*.pairsam
    do
      grep -E 'JJ|#|PP' $file > ${file}.JJPP
      grep -E 'JJ|#' $file > ${file}.JJ
      grep -E 'PP|#' $file > ${file}.PP
    done
	
    # Running pairsam to cooler conversion
    python 012_pairsam2cooler.py ../DATA/PAIR/${pref}_\*.pairsam $pref
    python 012_pairsam2cooler_savemarg.py $pref 
done
