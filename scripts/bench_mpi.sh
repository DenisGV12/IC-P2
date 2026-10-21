#!/usr/bin/env bash
# Práctica 3 — Tarea 3: mediciones del solver MPI frente al secuencial.
#
# Uso:   bash scripts/bench_mpi.sh "<tallas N>" "<nº procesos>" <pasos> [repeticiones]
# Ej.:   bash scripts/bench_mpi.sh "32 64 96 128" "1 2 3 4 6 8" 30 3
# Estrategias de halo (MODES, por defecto "0 1"): 0 = Sendrecv bloqueante, 1 = no bloqueante solapado.
#
# Genera (mediana de las repeticiones, criterio: tiempo total):
#   results/mpi.csv      N,pasos,P,modo,total,vel_step,dens_step,comm_max,comm_media,halo_max,gather_max,n_halo,n_gather,MB_halo,MB_gather
#   results/mpi_seq.csv  N,pasos,total,vel_step,dens_step   (secuencial ORIGINAL de la práctica 2, ./fluid)
#
# Varios ordenadores (sala de prácticas): mismo usuario, mismo directorio y ejecutable en todos los nodos.
#   HOSTFILE=hosts MPIFLAGS="-mca plm_rsh_no_tree_spawn 1 --oversubscribe" bash scripts/bench_mpi.sh ...
# (con HOSTFILE, P procesos se reparten según los slots del archivo; usad MPIFLAGS="--map-by node"
#  para forzar un proceso por máquina y medir el efecto de la red).
set -euo pipefail

SIZES=${1:-"32 64 96"}
PROCS=${2:-"1 2 3 4"}
STEPS=${3:-30}
REPS=${4:-3}
OUT=${OUT:-results}
MPIBIN=${MPIBIN:-./fluid_mpi}
SEQBIN=${SEQBIN:-./fluid}
MPIRUN=${MPIRUN:-mpirun}
MPIFLAGS=${MPIFLAGS:-}
MODES=${MODES:-"0 1"}
[ -n "${HOSTFILE:-}" ] && MPIFLAGS="-hostfile $HOSTFILE $MPIFLAGS"

mkdir -p "$OUT"
[ -x "$MPIBIN" ] || { echo "Falta $MPIBIN (make mpi)"; exit 1; }
[ -x "$SEQBIN" ] || { echo "Falta $SEQBIN (make)"; exit 1; }

median_row() {   # lee líneas CSV por stdin y devuelve la de mediana según la columna $1
	sort -t, -k"$1" -g | awk '{a[NR]=$0} END{print a[int((NR+1)/2)]}'
}

echo "N,pasos,P,modo,total,vel_step,dens_step,comm_max,comm_media,halo_max,gather_max,n_halo,n_gather,MB_halo,MB_gather" > "$OUT/mpi.csv"
echo "N,pasos,total,vel_step,dens_step" > "$OUT/mpi_seq.csv"

for n in $SIZES; do
	# secuencial original: salida "CSV,N,pasos,total,fuentes,vel,dens"
	row=$(for r in $(seq "$REPS"); do "$SEQBIN" "$n" "$STEPS" | grep '^CSV'; done | median_row 4)
	echo "$row" | awk -F, '{printf "%s,%s,%s,%s,%s\n",$2,$3,$4,$6,$7}' >> "$OUT/mpi_seq.csv"
	echo "[seq] N=$n  $(echo "$row" | awk -F, '{print $4" s"}')"

	for p in $PROCS; do
		if [ "$p" -gt "$n" ]; then echo "[mpi] N=$n P=$p omitido (P > N)"; continue; fi
		for m in $MODES; do
			[ "$p" -eq 1 ] && [ "$m" -ne 0 ] && continue     # con P=1 no hay comunicación: un solo modo
			row=$(for r in $(seq "$REPS"); do
				# shellcheck disable=SC2086
				$MPIRUN $MPIFLAGS -np "$p" "$MPIBIN" "$n" "$STEPS" 0.4 "$m" | grep '^CSV'
			done | median_row 6)
			echo "$row" | cut -d, -f2- >> "$OUT/mpi.csv"
			echo "[mpi] N=$n P=$p modo=$m  $(echo "$row" | awk -F, '{print $6" s  (comm media "$10" s)"}')"
		done
	done
done
echo "Listo: $OUT/mpi.csv y $OUT/mpi_seq.csv  (siguiente: make charts-mpi)"
