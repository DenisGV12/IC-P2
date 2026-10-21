#!/usr/bin/env python3
"""Práctica 3 — gráficas de la Tarea 3 a partir de results/mpi.csv y results/mpi_seq.csv.
Uso: python3 scripts/plot_mpi.py [carpeta_results] [carpeta_salida]
Genera: p3_speedup.png, p3_eficiencia.png, p3_comunicacion.png, p3_tiempos.png, p3_heatmap.png
(speed-up respecto al secuencial ORIGINAL de la P2, y respecto a la propia versión MPI con P=1).
Si hay mediciones con modo=1 (solape) también genera p3_estrategias.png."""
import csv, sys, os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = sys.argv[1] if len(sys.argv) > 1 else "results"
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(res, "fig")
os.makedirs(out, exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

mpi_all = list(csv.DictReader(open(os.path.join(res, "mpi.csv"))))
mpi = [r for r in mpi_all if int(r["modo"]) == 0 or int(r["P"]) == 1]   # estrategia A (bloqueante) para las figuras principales
seq = {int(r["N"]): float(r["total"]) for r in csv.DictReader(open(os.path.join(res, "mpi_seq.csv")))}
steps = int(mpi[0]["pasos"])

T, COMM, HALO, GATH = {}, {}, {}, {}
for r in mpi:
    n, p = int(r["N"]), int(r["P"])
    T[(n, p)] = float(r["total"]); COMM[(n, p)] = float(r["comm_media"])
    HALO[(n, p)] = float(r["halo_max"]); GATH[(n, p)] = float(r["gather_max"])
Ns = sorted({n for n, _ in T}); Ps = sorted({p for _, p in T})
have_p1 = all((n, 1) in T for n in Ns)

def sp_seq(n, p):  return seq[n] / T[(n, p)]
def sp_p1(n, p):   return T[(n, 1)] / T[(n, p)]

# 1) speed-up
fig, axs = plt.subplots(1, 2 if have_p1 else 1, figsize=(11 if have_p1 else 6, 4.2), squeeze=False)
ax = axs[0][0]
for n in Ns:
    ps = [p for p in Ps if (n, p) in T]
    ax.plot(ps, [sp_seq(n, p) for p in ps], "o-", label=f"N={n}")
ax.set_title(f"Speed-up frente al secuencial original ({steps} pasos)")
ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("T_seq / T_P"); ax.legend()
if have_p1:
    ax = axs[0][1]
    for n in Ns:
        ps = [p for p in Ps if (n, p) in T]
        ax.plot(ps, [sp_p1(n, p) for p in ps], "o-", label=f"N={n}")
    ax.plot(Ps, Ps, "--", color="gray", label="Ideal (S=P)")
    ax.set_title("Speed-up frente a la versión MPI con P=1")
    ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("T_1 / T_P"); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(out, "p3_speedup.png")); plt.close(fig)

# 2) eficiencia
fig, axs = plt.subplots(1, 2 if have_p1 else 1, figsize=(11 if have_p1 else 6, 4.2), squeeze=False)
ax = axs[0][0]
for n in Ns:
    ps = [p for p in Ps if (n, p) in T]
    ax.plot(ps, [sp_seq(n, p) / p for p in ps], "o-", label=f"N={n}")
ax.axhline(1, ls="--", color="gray"); ax.set_ylim(bottom=0)
ax.set_title("Eficiencia S_seq / P"); ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("Eficiencia"); ax.legend()
if have_p1:
    ax = axs[0][1]
    for n in Ns:
        ps = [p for p in Ps if (n, p) in T]
        ax.plot(ps, [sp_p1(n, p) / p for p in ps], "o-", label=f"N={n}")
    ax.axhline(1, ls="--", color="gray"); ax.set_ylim(bottom=0)
    ax.set_title("Eficiencia paralela S_1 / P"); ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("Eficiencia"); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(out, "p3_eficiencia.png")); plt.close(fig)

# 3) cálculo frente a comunicación/sincronización
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
nbig = Ns[-1]; ps = [p for p in Ps if (nbig, p) in T]
comp = [T[(nbig, p)] - COMM[(nbig, p)] for p in ps]; com = [COMM[(nbig, p)] for p in ps]
x = range(len(ps))
axs[0].bar(x, comp, label="Cálculo"); axs[0].bar(x, com, bottom=comp, label="Comunicación + sincronización")
axs[0].set_xticks(list(x)); axs[0].set_xticklabels(ps)
axs[0].set_title(f"Reparto del tiempo, N={nbig} (media de procesos)"); axs[0].set_xlabel("Procesos MPI (P)"); axs[0].set_ylabel("Tiempo (s)"); axs[0].legend()
for n in Ns:
    pp = [p for p in Ps if (n, p) in T and p > 1]
    axs[1].plot(pp, [100 * COMM[(n, p)] / T[(n, p)] for p in pp], "o-", label=f"N={n}")
axs[1].set_title("% del tiempo en comunicación + sincronización"); axs[1].set_xlabel("Procesos MPI (P)"); axs[1].set_ylabel("%"); axs[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(out, "p3_comunicacion.png")); plt.close(fig)

# 4) tiempos absolutos (escala log) con el secuencial como referencia
fig, ax = plt.subplots(figsize=(6, 4.2))
for n in Ns:
    ps = [p for p in Ps if (n, p) in T]
    l, = ax.plot(ps, [T[(n, p)] for p in ps], "o-", label=f"N={n}")
    ax.axhline(seq[n], ls=":", color=l.get_color())
ax.set_yscale("log"); ax.set_title("Tiempo total (línea de puntos: secuencial original)")
ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("Tiempo (s)"); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(out, "p3_tiempos.png")); plt.close(fig)

# 5) mapa de calor N x P
fig, ax = plt.subplots(figsize=(6, 4))
grid = [[sp_seq(n, p) if (n, p) in T else float("nan") for p in Ps] for n in Ns]
im = ax.imshow(grid, aspect="auto", origin="lower", cmap="viridis")
ax.set_xticks(range(len(Ps))); ax.set_xticklabels(Ps); ax.set_yticks(range(len(Ns))); ax.set_yticklabels(Ns)
for i, n in enumerate(Ns):
    for j, p in enumerate(Ps):
        if (n, p) in T: ax.text(j, i, f"{grid[i][j]:.1f}", ha="center", va="center", color="w", fontsize=8)
ax.grid(False); ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("N"); ax.set_title("Speed-up frente al secuencial original")
fig.colorbar(im); fig.tight_layout(); fig.savefig(os.path.join(out, "p3_heatmap.png")); plt.close(fig)
# 6) comparación de estrategias de halo: A (bloqueante) frente a B (solapado)
B = {(int(r["N"]), int(r["P"])): float(r["total"]) for r in mpi_all if int(r["modo"]) == 1}
if B:
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for n in Ns:
        ps = [p for p in Ps if (n, p) in B and (n, p) in T]
        l, = ax.plot(ps, [T[(n, p)] / B[(n, p)] for p in ps], "o-", label=f"N={n}")
    ax.axhline(1, ls="--", color="gray")
    ax.set_title("Ganancia de B (solape) sobre A (bloqueante)"); ax.set_xlabel("Procesos MPI (P)"); ax.set_ylabel("T_A / T_B"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "p3_estrategias.png")); plt.close(fig)
print("Gráficas en", out)
