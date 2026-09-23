#!/usr/bin/env bash
# Escalado del tiempo secuencial con el tamaño de malla (apartado 2.3)
# Uso: bash scripts/bench_sizes.sh [pasos] [repeticiones]
set -e
STEPS=${1:-20}; REPS=${2:-3}
SIZES="32 48 64 80 96 112 128"
mkdir -p results
make -s clean; make -s CFLAGS="-O2"
OUT=results/sizes.csv
echo "N,pasos,rep,total_s,fuentes_s,vel_step_s,dens_step_s" > $OUT
for n in $SIZES; do
  for r in $(seq 1 $REPS); do
    line=$(./fluid $n $STEPS | grep '^CSV' | cut -d, -f2-)
    IFS=, read -r nn st tot src vel den <<< "$line"
    echo "$nn,$st,$r,$tot,$src,$vel,$den" >> $OUT
    echo "N=$n rep=$r total=${tot}s"
  done
done
echo "Resultados en $OUT"
