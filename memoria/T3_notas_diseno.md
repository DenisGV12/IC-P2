# Práctica 3 · Tarea 3 — notas para la memoria (solver de fluidos 3D con MPI)

Borrador de los apartados «Identificación», «Implementación» y «Evaluación». Los huecos `[…]` se rellenan con los CSV de `scripts/bench_mpi.sh` (medidos en la sala de prácticas / vuestro equipo; no con estos borradores).

## 1. Identificación de la parte a paralelizar

- Por paso, el solver llama **6 veces a `lin_solve`** (3 `diffuse` de velocidad, 2 `project`, 1 `diffuse` de densidad) y cada llamada hace `LINEARSOLVERTIMES = 10` barridos de toda la malla → **60 barridos por paso**, más 4 `advect`, 2 `project` (bucles de divergencia/gradiente) y `add_source`. Todo es O(N³) por paso.
- gprof de la P2 (N=96, 30 pasos, `-O2`): `diffuse` ≈ 48 % del tiempo, camino crítico ≈ 53 %, speed-up máximo ≈ ×1,9. Esa cota venía de la **dependencia de Gauss-Seidel lexicográfico** (cada celda usa vecinos ya actualizados en el mismo barrido). La solución de la P2 era el esquema rojo-negro; aquí se aplica: elimina la dependencia dentro del barrido y la cota de ×1,9 deja de aplicar (a costa de cambiar el orden numérico, ver §4).
- Se paraleliza **todo el paso** (no solo `lin_solve`): si solo se repartiera el solver lineal, `advect`/`project` quedarían secuenciales y limitarían la ganancia (Amdahl).

## 2. Diseño

| Decisión | Elección | Motivo |
|---|---|---|
| Descomposición | Planos `k` contiguos (slabs), reparto equilibrado ±1 plano | `k` es la dimensión más lenta de `IX(i,j,k)`: cada plano es un bloque contiguo → los halos son un único mensaje sin empaquetar |
| Memoria | Cada proceso reserva la malla completa (mismos índices globales) y **solo calcula sus planos** | Mantiene el código idéntico al secuencial (`IX`, `set_bnd`, `advect`). Limitación: la memoria no escala con P (8 campos × (N+2)³ × 4 B por proceso) |
| `lin_solve` | Gauss-Seidel **rojo-negro**: medio barrido de un color, intercambio de halos, medio barrido del otro color, halos, `set_bnd` | Las celdas de un color solo dependen del otro color |
| `advect` | `MPI_Allgatherv` in-place del campo a advectar antes del bucle | Es semi-lagrangiano: lee `d0` en posiciones arbitrarias (hasta `dt·N·|u|` celdas lejos), no basta con un halo |
| `set_bnd` | Sin comunicación: cada proceso rellena el contorno de sus planos; rank 0 y último rank las caras z | Todos los valores de contorno salen de caras del mismo plano o de los planos k=1 / k=O |
| `project` | Local; solo necesita el halo de `w` (div) y de `p` (gradiente) | `u`,`v` solo se derivan en i,j (dentro del plano) |
| Escenario | Cada proceso pone a cero y rellena solo sus planos | Evita una parte secuencial replicada |
| Bucles | Orden k-j-i (i contiguo) en `lin_solve`, `advect`, `project` | Lo que se concluye de la Tarea 0 (recorrido de matrices): el original recorre con `k` en el bucle interno (paso de un plano entero) |

**Dos estrategias de comunicación de halos** (se elige con el 4.º argumento, `modo`):
- **A (modo 0)**: `MPI_Sendrecv` bloqueante tras cada medio barrido.
- **B (modo 1)**: se calculan primero los planos frontera, se lanzan `MPI_Isend/Irecv`, se calculan los planos interiores mientras viajan los mensajes y se espera con `MPI_Waitall`. Es correcto porque los interiores no leen los planos fantasma y los planos enviados ya no se modifican. Que el solape dé ganancia real depende de que la implementación MPI haga progreso asíncrono (con mensajes grandes de rendez-vous suele no ocurrir): es parte del análisis de la cuestión 4.

## 3. Funciones MPI usadas

| Función | Dónde | Papel |
|---|---|---|
| `MPI_Init/Finalize`, `MPI_Comm_rank/size` | main | Entorno, identificación (modelo SPMD) |
| `MPI_Sendrecv` | `halo_exchange` (A) | Punto a punto: intercambio de planos frontera con vecinos; `MPI_PROC_NULL` en los extremos |
| `MPI_Isend/Irecv/Waitall` | `halo_start/finish` (B) | Punto a punto no bloqueante con solape cálculo/comunicación |
| `MPI_Allgatherv` (`MPI_IN_PLACE`) | `allgather_field` | Colectiva: replica el campo para `advect` (trozos de tamaño ligeramente distinto) |
| `MPI_Barrier` | main | Sincroniza el inicio/fin de la región medida |
| `MPI_Reduce` (MAX, SUM) | main | Tiempos (el paralelo es el del proceso más lento), tiempos de comunicación, volúmenes |
| `MPI_Wtime` | main, halos | Medición |

Volumen de comunicación (por proceso interior y paso), con `plane = (N+2)²` floats:
- Halos: `lin_solve` hace 1 + 2·10 = 21 intercambios × 6 llamadas + 1 de `w` = **127 intercambios/paso**, cada uno envía 2 planos → `127 · 2 · (N+2)² · 4 B`. No depende de P (sí el nº de mensajes por proceso: 2 vecinos).
- Allgather: **4 por paso** (u0, v0, w0 y densidad). Cada proceso recibe ≈ `(P−1)/P · (N+2)³ · 4 B` por campo: **no decrece al aumentar P** → es el término que limita la escalabilidad.
- Contadores medidos: `n_halo`, `n_gather`, `MB_halo`, `MB_gather` en `results/mpi.csv`.

## 4. Verificación de corrección (reproducible)

El resultado MPI es **idéntico para cualquier P** y para ambos modos (cada celda solo depende de sus vecinos del otro color). Referencia: el código original con *solo* `lin_solve` cambiado a rojo-negro (secuencial). Prueba hecha con N=32, 30 pasos:

| Versión | checksum dens | checksum vel |
|---|---|---|
| MPI P=1,2,3,4,5 (ambos modos) | 6,234267e+02 | 5,336001e+01 |
| Secuencial rojo-negro | 6,234267e+02 | 5,336001e+01 |
| Secuencial original (lexicográfico) | 7,446310e+02 | 7,513547e+01 |

También coincide con reparto desigual (N=33, P=4 y 7) y con 1 plano por proceso (N=8, P=8). La diferencia con el original **es esperada**: con 10 iteraciones el solver de Poisson no converge y rojo-negro y lexicográfico dan iterados distintos (ambos válidos). Hay que decirlo en la memoria y aclarar que el checksum de la versión paralela se compara con la rojo-negro, no con la original.

## 5. Cómo medir (para que las gráficas sean defendibles)

1. `make && make mpi`, y copiar `fluid_mpi` (y `hosts`) al mismo directorio en todos los nodos, mismo usuario (Tarea 2).
2. Un solo nodo (multinúcleo): `make bench-mpi SIZES="32 64 96 128" PROCS="1 2 4 8" STEPS=30 REPS=3`.
3. Varios nodos: `HOSTFILE=hosts MPIFLAGS="-mca plm_rsh_no_tree_spawn 1" make bench-mpi …` y, para aislar el efecto de la red, repetir con `MPIFLAGS="--map-by node"` (un proceso por máquina).
4. `make charts-mpi` → `results/fig/p3_*.png`. Caracterizar la máquina con `lscpu` y `free -h` (nº de procesadores/núcleos, cachés, memoria) para el apartado «Sobre los resultados».
5. Mantened las máquinas sin otra carga durante la medición y no mezcléis resultados de equipos distintos en la misma curva (guardad un `results_<equipo>/` por cada uno y pasad la carpeta a `plot_mpi.py`).

**Dos referencias de speed-up, ambas en las gráficas:**
- `T_seq_original / T_P`: ganancia total respecto a la entrega de la P2. **Incluye** el efecto del cambio de orden de bucles (mejor acceso a memoria) y de rojo-negro, además del paralelismo.
- `T_{MPI,P=1} / T_P`: ganancia estrictamente paralela (mismo algoritmo). Es la que corresponde a «speed-up» y «eficiencia paralela» de las cuestiones 1, 2 y 6.

## 6. Guion para las cuestiones de «Evaluación»

1. **Speed-up** (`p3_speedup.png`, `p3_heatmap.png`): en función de P y N (N = parámetro de tamaño; pasos solo escala linealmente).
2. **Escalabilidad** (`p3_speedup.png`, `p3_eficiencia.png`): ¿mejora siempre? Esperable: no, a N pequeño la comunicación domina; límites: allgather constante en P, latencia de los 127 halos/paso, desequilibrio ±1 plano cuando P no divide a N, saturación del ancho de banda de memoria dentro de un nodo.
3. **Comunicación** (`p3_comunicacion.png`): `comm_media / total`. Ojo: incluye también las esperas por sincronización (un `Sendrecv` espera al vecino más lento), no solo la transferencia; se puede separar halo/allgather con `halo_max` y `gather_max`.
4. **Estrategias** (`p3_estrategias.png`): A (bloqueante) frente a B (solape), y, si queréis, 1 proceso por nodo frente a varios por nodo. Justificad con los tiempos de comunicación y con N.
5. **Ganancia máxima teórica**: parte secuencial del cálculo ≈ nula (todas las fases se reparten), por lo que el límite de Amdahl del cálculo es P; el límite real lo marcan la comunicación (allgather O(N³) por proceso independiente de P), la latencia y el desequilibrio. Comparad con la cota ×1,9 de la P2 y explicad por qué rojo-negro la rompe.
6. **Eficiencia secuencial frente a paralela** (`p3_eficiencia.png`): `S/P` para cada (N, P); indicad la mejor configuración y a partir de qué N compensa paralelizar.

## 7. Limitaciones a mencionar
- Memoria replicada (no distribuida del todo): válido para estas mallas (N=160 → ≈136 MB/proceso), pero no escalaría en memoria. Una versión totalmente distribuida guardaría solo slab + halos y reemplazaría el allgather por un intercambio acotado por el desplazamiento máximo del `advect`.
- Descomposición 1D: el volumen de halo no baja con P; 2D/3D lo reduciría pero complica `set_bnd` y los halos no contiguos.
- Con `diff = visc = 0` (valores del escenario) `diffuse` equivale a la identidad (a=0 → x = x0) pero sigue ejecutando los 10 barridos y sus halos: se mantiene así para que el código sea general y comparable con el secuencial.
