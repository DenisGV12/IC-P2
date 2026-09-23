#!/usr/bin/env bash
# Comparación de opciones de compilación (apartados 2.4 y 2.5)
# Uso: bash scripts/bench_flags.sh [N] [pasos] [repeticiones]
set -e
N=${1:-96}; STEPS=${2:-30}; REPS=${3:-3}
mkdir -p results
OUT=results/flags.csv
echo "config,cflags,rep,total_s,vel_step_s,dens_step_s,checksum_dens" > $OUT
CONFIGS=(
  "O0|-O0"
  "O1|-O1"
  "O2|-O2"
  "O2_novec|-O2 -fno-tree-vectorize"
  "O2_ftz|-O2 -mdaz-ftz"
  "O3|-O3"
  "O3_novec|-O3 -fno-tree-vectorize"
  "O3_native|-O3 -march=native"
  "O3_native_fastmath|-O3 -march=native -ffast-math"
  "Ofast_native_unroll|-Ofast -march=native -funroll-loops"
)
for c in "${CONFIGS[@]}"; do
  name=${c%%|*}; flags=${c#*|}
  gcc $flags -Wall -Wno-misleading-indentation main.c solver3d.c -o fluid_$name -lm
  echo "== $name: gcc $flags"
  for r in $(seq 1 $REPS); do
    res=$(./fluid_$name $N $STEPS)
    csv=$(echo "$res" | grep '^CSV' | cut -d, -f2-)
    IFS=, read -r nn st tot src vel den <<< "$csv"
    ck=$(echo "$res" | grep checksum | sed 's/.*dens=\([^ ]*\).*/\1/')
    echo "$name,\"$flags\",$r,$tot,$vel,$den,$ck" >> $OUT
    echo "   rep $r: ${tot}s  checksum=$ck"
  done
done
echo "Resultados en $OUT"
