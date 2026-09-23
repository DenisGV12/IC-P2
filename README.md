# Práctica 2 IC — Solver de fluidos 3D (Jos Stam), versión sin interfaz

Código base: https://github.com/BlainMaguire/3dfluid (MIT, commit `1e58a59`), basado en
Jos Stam, "Real-Time Fluid Dynamics for Games", GDC 2003.

## Ficheros
- `main.c` — programa principal sin interfaz: parámetros, guion de fuentes/fuerzas, medición y checksums.
- `solver3d.c` — solver del original, sin cambios algorítmicos.
- `Makefile` — `make`, `make run`, `make clean` + objetivos de análisis.
- `scripts/` — mediciones (`bench_sizes.sh`, `bench_flags.sh`) y gráficas (`plot_results.py`).
- `memoria/` — memoria, figuras y datos medidos.

## Uso
```bash
make                 # compila con gcc -O2
make run             # caso de referencia: N=96, 30 pasos (~17 s)
./fluid N pasos [dt]
make vec             # informe de autovectorización -> results/
make asm             # ensamblador del solver -> results/
make bench-sizes     # escalado con N -> results/sizes.csv
make bench-flags     # opciones de compilación -> results/flags.csv
make charts          # gráficas -> results/fig/
make clean
```

## Changelog (modificaciones respecto al original)
- `main.c` reescrito sin GLUT/OpenGL: sin ventana, render, teclado ni ratón.
- Parámetros por línea de comandos: N (malla N³), pasos, dt.
- Teclas X/D/W/S sustituidas por un guion determinista (`scenario_events`): fuente en el centro
  cada 5 pasos y fuerza alternando ejes x/y/z cada 15 pasos.
- Medición con `clock_gettime(CLOCK_MONOTONIC)`, total y por fase; checksums para validar resultados.
- `solver3d.c`: finales de línea CRLF→LF y `LINEARSOLVERTIMES` configurable con `-D`.
- Makefile y scripts de medición nuevos.
