GENOME="dm3"
GENOMEDIR="../DATA/GENOME/"
PAIRDIR="../DATA/PAIR/"
COOLDIR="../DATA/COOL/"
STATSDIR="../DATA/STATS/"

pairsamtools merge --nproc 1 -o ${PAIRDIR}/Dros.pairsam.full.sorted ${PAIRDIR}/Dros_HiC_9_S14.pairsam.full.sorted ${PAIRDIR}/Dros_HiC_10_S15.pairsam.full.sorted

pairsamtools dedup \
    --max-mismatch 1 \
    --mark-dups \
    --output ${PAIRDIR}/Dros.pairsam.full.sorted.nodups.gz \
    --output-stats ${STATSDIR}/Dros.txt.full.dedup \
    ${PAIRDIR}/Dros.pairsam.full.sorted

cooler csort -c1 2 -c2 4 -p1 3 -p2 5 ${PAIRDIR}/Dros.pairsam.full.sorted.nodups.gz $GENOMEDIR/${GENOME}.reduced.chrom.sizes

for res in 1000 10000 20000 100000
do
  res_short=$(( res/1000 ))
  echo ${COOLDIR}/Dros.${res_short}.cool.full ${PAIRDIR}/Dros.pairsam.full.sorted.nodups.gz
  cooler cload pairix --nproc 16 --assembly ${GENOME} \
    $GENOMEDIR/${GENOME}.reduced.chrom.sizes:${res} ${PAIRDIR}/Dros.pairsam.full.sorted.nodups.blksrt.gz ${COOLDIR}/Dros.${res_short}.cool.full

  cooler balance --nproc 4 ${COOLDIR}/Dros.${res_short}.cool.full
done
