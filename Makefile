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

.PHONY: all run clean vec asm bench-sizes bench-flags charts

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

clean:
	rm -f $(BIN) $(BIN)_vec $(BIN)_O* *.o
