#!/usr/bin/env python3
"""Genera las gráficas de la memoria a partir de results/sizes.csv y results/flags.csv.
Uso: python3 scripts/plot_results.py [carpeta_results] [carpeta_salida]
Requiere matplotlib (sudo apt install python3-matplotlib)."""
import csv, sys, os, statistics as st
from collections import defaultdict, OrderedDict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = sys.argv[1] if len(sys.argv) > 1 else "results"
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(res, "fig")
os.makedirs(out, exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

def load(path):
    if not os.path.exists(path): return []
    return list(csv.DictReader(open(path)))

# ---------- tamaños ----------
rows = load(os.path.join(res, "sizes.csv"))
if rows:
    g = defaultdict(list)
    for r in rows: g[int(r["N"])].append(r)
    Ns = sorted(g); steps = int(rows[0]["pasos"])
    med = lambda n, k: st.median(float(r[k]) for r in g[n])
    tot = [med(n, "total_s") for n in Ns]
    vel = [med(n, "vel_step_s") for n in Ns]
    den = [med(n, "dens_step_s") for n in Ns]
    lo  = [min(float(r["total_s"]) for r in g[n]) for n in Ns]
    hi  = [max(float(r["total_s"]) for r in g[n]) for n in Ns]

    # 1) tiempo vs N con referencia N^3
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(Ns, tot, yerr=[[t-l for t, l in zip(tot, lo)], [h-t for t, h in zip(tot, hi)]],
                marker="o", capsize=3, label="Medido (mediana, min-max)")
    ref = [tot[0] * (n / Ns[0])**3 for n in Ns]
    ax.plot(Ns, ref, "--", color="gray", label=f"Escalado ideal ∝ N³ (desde N={Ns[0]})")
    ax.set_xlabel("N (malla N×N×N)"); ax.set_ylabel("Tiempo (s)")
    ax.set_title(f"Tiempo secuencial frente a tamaño ({steps} pasos, -O2)")
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(out, "tiempo_vs_N.png")); plt.close(fig)

    # 2) coste por celda y paso
    fig, ax = plt.subplots(figsize=(6, 4))
    cpc = [t / (n**3 * steps) * 1e6 for t, n in zip(tot, Ns)]
    ax.plot(Ns, cpc, marker="o", color="tab:red")
    ax2 = ax.twiny(); ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(Ns); ax2.set_xticklabels([f"{(n+2)**3*32/2**20:.0f}" for n in Ns], fontsize=8)
    ax2.set_xlabel("Memoria de los 8 arrays (MB)")
    ax.set_xlabel("N"); ax.set_ylabel("µs por celda y paso")
    ax.set_title("Coste por celda: crece al dejar de caber en caché", pad=28)
    fig.tight_layout(); fig.savefig(os.path.join(out, "coste_por_celda.png")); plt.close(fig)

    # 3) desglose por fase
    fig, ax = plt.subplots(figsize=(6, 4))
    pv = [v / t * 100 for v, t in zip(vel, tot)]; pd = [d / t * 100 for d, t in zip(den, tot)]
    x = range(len(Ns))
    ax.bar(x, pv, label="vel_step"); ax.bar(x, pd, bottom=pv, label="dens_step")
    ax.set_xticks(list(x)); ax.set_xticklabels(Ns); ax.set_ylim(0, 105)
    ax.set_xlabel("N"); ax.set_ylabel("% del tiempo total"); ax.set_title("Reparto del tiempo por fase")
    ax.legend(loc="lower right"); fig.tight_layout(); fig.savefig(os.path.join(out, "fases.png")); plt.close(fig)

# ---------- flags ----------
rows = load(os.path.join(res, "flags.csv"))
if rows:
    g = OrderedDict()
    for r in rows: g.setdefault(r["config"], []).append(float(r["total_s"]))
    names = list(g); meds = [st.median(v) for v in g.values()]
    base = meds[names.index("O2")] if "O2" in names else meds[0]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.barh(names[::-1], meds[::-1], color=["tab:orange" if n == "O2" else "tab:blue" for n in names[::-1]])
    for b, m in zip(bars, meds[::-1]):
        ax.text(b.get_width(), b.get_y() + b.get_height()/2, f" {m:.1f} s (×{base/m:.2f})", va="center", fontsize=8)
    ax.set_xlabel("Tiempo (s), mediana"); ax.set_title("Opciones de compilación (ganancia respecto a -O2)")
    ax.set_xlim(0, max(meds) * 1.3); fig.tight_layout(); fig.savefig(os.path.join(out, "flags.png")); plt.close(fig)

print("Gráficas en", out)
