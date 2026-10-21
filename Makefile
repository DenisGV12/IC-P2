# Práctica 2 IC — Solver de fluidos 3D (Jos Stam), versión sin interfaz
CC      := gcc
CFLAGS  ?= -O2
WARN    := -Wall -Wextra -Wno-misleading-indentation
LDLIBS  := -lm
SRCS    := main.c solver3d.c
BIN     := fluid

# Caso de prueba de referencia: N=96, 30 pasos, -O2 (~17 s)
N       ?= 96
STEPS   ?= 30

# --- Práctica 3 (MPI) ---
MPICC     ?= mpicc
MPI_BIN   := fluid_mpi
MPI_SRCS  := main_mpi.c solver3d_mpi.c
NP        ?= 4
SIZES     ?= 32 64 96 128
PROCS     ?= 1 2 3 4
REPS      ?= 3

.PHONY: all run clean vec asm bench-sizes bench-flags charts mpi run-mpi bench-mpi charts-mpi gprof

all: $(BIN)

$(BIN): $(SRCS)
	$(CC) $(CFLAGS) $(WARN) $(SRCS) -o $@ $(LDLIBS)

run: $(BIN)
	./$(BIN) $(N) $(STEPS)

# Informe de autovectorización
vec:
	mkdir -p results
	$(CC) -O3 -march=native $(WARN) -fopt-info-vec-optimized=results/vec_ok.txt $(SRCS) -o $(BIN)_vec $(LDLIBS)
	$(CC) -O3 -march=native $(WARN) -fopt-info-vec-missed=results/vec_missed.txt $(SRCS) -o $(BIN)_vec $(LDLIBS)
	@echo "== Bucles vectorizados:"; sort -u results/vec_ok.txt | grep 'loop vectorized'

# Ensamblador del solver
asm:
	mkdir -p results
	$(CC) -O2 -S -fverbose-asm solver3d.c -o results/solver3d_O2.s
	$(CC) -O3 -march=native -S -fverbose-asm solver3d.c -o results/solver3d_O3_native.s
	@echo "Instrucciones vectoriales empaquetadas (ps) en O3+native:"
	@grep -cE '\bv?[a-z]+ps\b' results/solver3d_O3_native.s || true

# Mediciones
bench-sizes:
	bash scripts/bench_sizes.sh

bench-flags:
	bash scripts/bench_flags.sh $(N) $(STEPS)

# Gráficas desde results/*.csv (requiere python3-matplotlib)
charts:
	python3 scripts/plot_results.py results results/fig

# Versión MPI: mpirun -np $(NP) ./fluid_mpi N pasos
mpi: $(MPI_BIN)

$(MPI_BIN): $(MPI_SRCS) slab.h
	$(MPICC) $(CFLAGS) $(WARN) $(MPI_SRCS) -o $@ $(LDLIBS)

run-mpi: $(MPI_BIN)
	mpirun -np $(NP) ./$(MPI_BIN) $(N) $(STEPS)

# Mediciones MPI frente al secuencial: make bench-mpi SIZES="32 64 96" PROCS="1 2 4" STEPS=30
# (varios ordenadores: HOSTFILE=hosts MPIFLAGS="-mca plm_rsh_no_tree_spawn 1" make bench-mpi)
bench-mpi: $(BIN) $(MPI_BIN)
	bash scripts/bench_mpi.sh "$(SIZES)" "$(PROCS)" $(STEPS) $(REPS)

charts-mpi:
	python3 scripts/plot_mpi.py results results/fig

# Tarea 1: perfil gprof del solver secuencial, -O0 y -O2 (informes en results/gprof_*.txt)
gprof:
	mkdir -p results
	for o in O0 O2; do \
	  mkdir -p results/prof_$$o && cd results/prof_$$o && \
	  $(CC) -$$o -g -pg ../../main.c ../../solver3d.c -o fluid_prof -lm && \
	  ./fluid_prof $(N) $(STEPS) > /dev/null && \
	  gprof -b -p -q ./fluid_prof gmon.out > ../gprof_$$o.txt && cd ../..; \
	done
	@echo "Informes: results/gprof_O0.txt results/gprof_O2.txt"

clean:
	rm -f $(BIN) $(BIN)_vec $(BIN)_O* $(MPI_BIN) *.o
