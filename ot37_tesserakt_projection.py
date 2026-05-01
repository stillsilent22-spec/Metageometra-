"""
OT-37: Tesseract Fraktal-Projektion auf S³ — Hausdorff-Dimensionstest
======================================================================
Testet: D_H(pi(T_4D)) = 0.77 = D_f,geo des HTM-Frameworks

Methode:
  1. IFS (Iterated Function System) auf dem 4D-Einheitstesserakt
     f_i(x) = r * R(Delta-theta) * x + t_i
     mit HTM-Torsionswinkel Delta-theta, Skalierung r, N_sub Kopien

  2. Projektion pi: R^4 -> S³:  pi(x) = x / ||x||

  3. Box-Counting auf S³ (via winkelbasierte Abdeckung in R^4)

  4. Parameter-Scan: r in [0.3..0.7], Delta-theta in [50..70] deg, N_sub in {4,8,16}

  5. Falsifikationstest: Liegt D_H ≈ 0.77 bei Delta-theta ≈ 58.7 deg?

Metageometra V19.0 / HTM Framework — Kevin Hannemann 2026
"""

import os, sys, math, time
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

# ── HTM-Konstanten ────────────────────────────────────────────────────────────
DF_GEO_TARGET = 0.77
THETA_HTM     = 58.7   # Grad

# ── IFS: Rotationsmatrix in der (x0,x1)-Ebene im R^4 ─────────────────────────
def rotation_4d(theta_deg):
    """4x4 Rotationsmatrix in der (x0,x1)-Ebene."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    R = np.eye(4)
    R[0, 0] =  c;  R[0, 1] = -s
    R[1, 0] =  s;  R[1, 1] =  c
    return R

# ── IFS-Translationsvektoren: Teilmenge der Tesserakt-Ecken ──────────────────
def get_translations(n_sub):
    """
    n_sub Translationsvektoren aus den 16 Ecken des Einheitstesserakts {±1}^4.
    Wählt gleichmäßig verteilte Ecken (4, 8 oder 16).
    """
    all_corners = np.array([[s0, s1, s2, s3]
                             for s0 in (-1, 1) for s1 in (-1, 1)
                             for s2 in (-1, 1) for s3 in (-1, 1)],
                            dtype=float)  # shape (16, 4)
    if n_sub >= 16:
        return all_corners
    # Gleichmäßige Unterauswahl via Linspace-Index
    idx = np.round(np.linspace(0, 15, n_sub)).astype(int)
    return all_corners[idx]

# ── IFS Iteration mit Chaos-Game (random IFS) ─────────────────────────────────
def ifs_fractal(r, theta_deg, n_sub, n_points=30_000, n_warmup=200, seed=42):
    """
    Erzeugt Punkte des IFS-Fraktals im R^4 via Chaos-Game.
    Jede Kontraktion: f_i(x) = r * R * x + (1-r) * t_i
    (Translation skaliert mit (1-r) damit Fixpunkte bei t_i liegen)
    """
    rng = np.random.default_rng(seed)
    R   = rotation_4d(theta_deg)
    T   = get_translations(n_sub)          # (n_sub, 4)

    x = np.zeros(4)
    points = np.empty((n_points, 4))

    for i in range(n_warmup + n_points):
        k = rng.integers(n_sub)
        x = r * (R @ x) + (1.0 - r) * T[k]
        if i >= n_warmup:
            points[i - n_warmup] = x

    return points  # (n_points, 4)

# ── Projektion auf S³ ─────────────────────────────────────────────────────────
def project_s3(points):
    """Normiert jeden Punkt auf die Einheits-S³."""
    norms = np.linalg.norm(points, axis=1, keepdims=True)
    # Punkte nahe 0 vermeiden (sollten nicht vorkommen beim IFS)
    mask = (norms[:, 0] > 1e-12)
    return points[mask] / norms[mask], mask

# ── Box-Counting auf S³ (via R^4-Winkelboxen) ────────────────────────────────
def box_count_s3(pts_s3, eps_list):
    """
    Box-Counting für Punkte auf S³ ⊂ R^4.
    Boxgröße eps in Einheiten im Bereich [-1,1] je Achse.
    Nutzt np.unique für Effizienz.
    """
    counts = []
    for eps in eps_list:
        bins = np.floor((pts_s3 + 1.0) / eps).astype(np.int32)
        # np.unique auf 2D-Array: jede Zeile = ein Bin-Tuple
        unique_rows = np.unique(bins, axis=0)
        counts.append(len(unique_rows))
    return np.array(counts, dtype=float)

# ── Hausdorff-Dimension via Log-Log-Fit ────────────────────────────────────────
def hausdorff_dim(eps_list, counts):
    """
    D_H = Steigung von log(N) vs log(1/eps).
    Verwendet nur Punkte mit counts > 0 und plausiblen Werten.
    """
    eps_arr = np.array(eps_list)
    cnt_arr = counts.copy()
    mask = cnt_arr > 1
    if mask.sum() < 3:
        return np.nan, np.nan, np.nan
    log_inv_eps = np.log(1.0 / eps_arr[mask])
    log_N       = np.log(cnt_arr[mask])
    slope, intercept, r, p, se = stats.linregress(log_inv_eps, log_N)
    return slope, r**2, se

# ── Parameter-Scan ────────────────────────────────────────────────────────────
def run_scan():
    r_values      = np.arange(0.30, 0.71, 0.10)                    # 5 Werte
    theta_values  = np.arange(50.0, 71.0, 5.0)                     # 5 Werte
    nsub_values   = [4, 8, 16]

    # Box-Größen: 2 Dekaden von eps=0.04 bis eps=0.5
    eps_list = list(np.logspace(np.log10(0.040), np.log10(0.50), 16))

    results = []   # (r, theta, n_sub, D_H, R2, se)

    total = len(r_values) * len(theta_values) * len(nsub_values)
    done  = 0

    print(f"\n  Scan: {len(r_values)} r × {len(theta_values)} theta × {len(nsub_values)} N_sub "
          f"= {total} Kombinationen")
    print(f"  {'r':>5}  {'theta':>6}  {'N_sub':>5}  {'D_H':>6}  {'R²':>6}  {'Delta':>7}")
    print("  " + "-" * 48)

    t0 = time.time()
    for r in r_values:
        for theta in theta_values:
            for n_sub in nsub_values:
                pts4d = ifs_fractal(r, theta, n_sub, n_points=30_000)
                pts_s3, _ = project_s3(pts4d)
                counts = box_count_s3(pts_s3, eps_list)
                D_H, R2, se = hausdorff_dim(eps_list, counts)
                delta = D_H - DF_GEO_TARGET if not math.isnan(D_H) else float('nan')
                results.append((r, theta, n_sub, D_H, R2, se))
                done += 1
                print(f"  {r:5.2f}  {theta:6.1f}  {n_sub:5d}  "
                      f"{D_H:6.3f}  {R2:6.4f}  {delta:+7.3f}")

    elapsed = time.time() - t0
    print(f"\n  Scan abgeschlossen in {elapsed:.1f}s")
    return results, eps_list

# ── Hauptanalyse ──────────────────────────────────────────────────────────────
print("=" * 70)
print("  OT-37: Tesseract-Projektion auf S³ — Hausdorff-Dimensionstest")
print("=" * 70)
print(f"  Zielwert: D_f,geo = {DF_GEO_TARGET}  (HTM V19.0)")
print(f"  HTM-Torsionswinkel: Delta-theta = {THETA_HTM} deg")

results, eps_list = run_scan()

# ── Auswertung: HTM-Zone ──────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  HTM-KONSISTENZ-ANALYSE")
print("=" * 70)

# Alle Treffer mit D_H in [0.72, 0.82]
hits = [(r, th, ns, D, R2, se) for (r, th, ns, D, R2, se) in results
        if not math.isnan(D) and abs(D - DF_GEO_TARGET) <= 0.05]

# Treffer bei theta nahe 58.7 ± 5 deg
htm_hits = [(r, th, ns, D, R2, se) for (r, th, ns, D, R2, se) in hits
            if abs(th - THETA_HTM) <= 5.0]

print(f"\n  Gesamt-Treffer |D_H - 0.77| <= 0.05: {len(hits)} / {len(results)}")
print(f"  davon bei Delta-theta = 58.7 ± 5 deg: {len(htm_hits)}")

if htm_hits:
    print("\n  HTM-PARAMETER-TREFFER (theta nahe 58.7 deg, D_H nahe 0.77):")
    print(f"  {'r':>5}  {'theta':>6}  {'N_sub':>5}  {'D_H':>6}  {'R²':>6}  {'Delta':>7}")
    print("  " + "-" * 48)
    for (r, th, ns, D, R2, se) in sorted(htm_hits, key=lambda x: abs(x[3]-DF_GEO_TARGET)):
        print(f"  {r:5.2f}  {th:6.1f}  {ns:5d}  {D:6.3f}  {R2:6.4f}  {D-DF_GEO_TARGET:+7.3f}")
    verdict = "BESTANDEN"
    verdict_detail = (f"D_H = {htm_hits[0][3]:.3f} bei r={htm_hits[0][0]:.2f}, "
                      f"theta={htm_hits[0][1]:.1f} deg, N_sub={htm_hits[0][2]}")
else:
    print("\n  KEIN Treffer bei HTM-Winkel 58.7 ± 5 deg mit D_H ≈ 0.77")
    # Beste verfügbare Kombination bei 58.7 ± 10 deg
    near_theta = [(r, th, ns, D, R2, se) for (r, th, ns, D, R2, se) in results
                  if not math.isnan(D) and abs(th - THETA_HTM) <= 10.0]
    if near_theta:
        best = min(near_theta, key=lambda x: abs(x[3] - DF_GEO_TARGET))
        print(f"  Nächster Wert bei theta ± 10 deg: D_H = {best[3]:.3f} "
              f"(r={best[0]:.2f}, theta={best[1]:.1f}, N_sub={best[2]})")
    verdict = "FALSIFIZIERT"
    verdict_detail = "Kein (r, Delta-theta, N_sub) liefert D_H ≈ 0.77 bei Delta-theta ≈ 58.7 deg"

print(f"\n  OT-37 VERDICT: {verdict}")
print(f"  {verdict_detail}")

# ── Beste Parameter für Visualisierung ────────────────────────────────────────
# Suche Parameter die D_H am nächsten an 0.77 bringen, bevorzugt bei theta=58.7
if htm_hits:
    best_r, best_theta, best_nsub = htm_hits[0][0], htm_hits[0][1], htm_hits[0][2]
else:
    # Fallback: global bester
    best = min(results, key=lambda x: abs(x[3] - DF_GEO_TARGET) if not math.isnan(x[3]) else 99)
    best_r, best_theta, best_nsub = best[0], best[1], best[2]

print(f"\n  Beste Parameter für Visualisierung: r={best_r}, theta={best_theta}, N_sub={best_nsub}")

# ── Visualisierung: Stereografische Projektion R^4 -> R^3 ─────────────────────
print("\n  Erstelle Projektion...")
pts4d_vis = ifs_fractal(best_r, best_theta, best_nsub, n_points=30_000)
pts_s3_vis, _ = project_s3(pts4d_vis)

# Stereografische Projektion vom Nordpol (0,0,0,1):
# (x0,x1,x2,x3) -> (x0/(1-x3), x1/(1-x3), x2/(1-x3))
mask_stereo = np.abs(pts_s3_vis[:, 3]) < 0.98   # Nordpol ausschließen
pts_vis = pts_s3_vis[mask_stereo]
denom = 1.0 - pts_vis[:, 3]
X = pts_vis[:, 0] / denom
Y = pts_vis[:, 1] / denom
Z = pts_vis[:, 2] / denom

fig = plt.figure(figsize=(10, 8))
ax  = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z, s=0.3, alpha=0.3, c=Z, cmap='plasma', linewidths=0)
ax.set_title(f'OT-37: Fraktal-Tesserakt → S³ (stereografisch)\n'
             f'r={best_r}, Δθ={best_theta}°, N_sub={best_nsub}  '
             f'D_H={min(results, key=lambda x: abs(x[0]-best_r)+abs(x[1]-best_theta)+abs(x[2]-best_nsub))[3]:.3f}',
             fontsize=11)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.tick_params(labelsize=8)

plt.tight_layout()
vis_path = os.path.join(RESULTS, "OT_37_projection.png")
plt.savefig(vis_path, dpi=120)
plt.close()
print(f"  Visualisierung gespeichert: {vis_path}")

# ── Log-Log-Plot für bestes Ergebnis ──────────────────────────────────────────
best_entry = min(results, key=lambda x: (abs(x[0]-best_r) + abs(x[1]-best_theta) + abs(x[2]-best_nsub)))
pts4d_ll = ifs_fractal(best_r, best_theta, best_nsub, n_points=30_000)
pts_s3_ll, _ = project_s3(pts4d_ll)
counts_ll = box_count_s3(pts_s3_ll, eps_list)

eps_arr = np.array(eps_list)
mask_ll = counts_ll > 1
slope_ll, R2_ll, se_ll = hausdorff_dim(eps_list, counts_ll)

fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.scatter(np.log10(1.0/eps_arr[mask_ll]), np.log10(counts_ll[mask_ll]),
            color='steelblue', s=20, label='Box-Counting')
x_fit = np.array([np.log10(1.0/eps_arr[mask_ll]).min(),
                  np.log10(1.0/eps_arr[mask_ll]).max()])
y_fit = slope_ll * x_fit + np.log10(counts_ll[mask_ll]).mean() - slope_ll * np.log10(1.0/eps_arr[mask_ll]).mean()
ax2.plot(x_fit, y_fit, 'r-', label=f'Fit: D_H = {slope_ll:.3f}  (R²={R2_ll:.3f})')
ax2.axhline(y=np.nan, color='gray', linestyle='--')
ax2.set_xlabel('log₁₀(1/ε)')
ax2.set_ylabel('log₁₀ N(ε)')
ax2.set_title(f'OT-37 Box-Counting: r={best_r}, Δθ={best_theta}°, N_sub={best_nsub}')
ax2.legend()
plt.tight_layout()
loglog_path = os.path.join(RESULTS, "OT_37_boxcount.png")
plt.savefig(loglog_path, dpi=120)
plt.close()
print(f"  Log-Log-Plot gespeichert: {loglog_path}")

# ── Ergebnistext speichern ────────────────────────────────────────────────────
lines = [
    "=" * 70,
    "OT-37: Tesseract-Projektion auf S³ — Hausdorff-Dimensionstest (V19.0)",
    "=" * 70,
    "",
    f"Zielwert:             D_f,geo = {DF_GEO_TARGET}  (HTM Framework)",
    f"HTM-Torsionswinkel:   Delta-theta = {THETA_HTM} deg",
    f"Methode:              IFS Chaos-Game + S³-Projektion + Box-Counting",
    "",
    "-" * 70,
    "PARAMETER-SCAN ERGEBNISSE",
    "-" * 70,
    f"  {'r':>5}  {'theta':>6}  {'N_sub':>5}  {'D_H':>6}  {'R²':>6}  {'Delta':>7}",
    "  " + "-" * 48,
]
for (r, th, ns, D, R2, se) in results:
    delta_str = f"{D-DF_GEO_TARGET:+7.3f}" if not math.isnan(D) else "   nan"
    D_str     = f"{D:6.3f}" if not math.isnan(D) else "   nan"
    R2_str    = f"{R2:6.4f}" if not math.isnan(R2) else "   nan"
    lines.append(f"  {r:5.2f}  {th:6.1f}  {ns:5d}  {D_str}  {R2_str}  {delta_str}")

lines += [
    "",
    "-" * 70,
    "HTM-KONSISTENZ-ANALYSE",
    "-" * 70,
    "",
    f"  Gesamt-Treffer |D_H - 0.77| <= 0.05: {len(hits)} / {len(results)}",
    f"  davon bei Delta-theta = 58.7 +/- 5 deg: {len(htm_hits)}",
    "",
]

if htm_hits:
    lines += [
        "  HTM-PARAMETER-TREFFER:",
        f"  {'r':>5}  {'theta':>6}  {'N_sub':>5}  {'D_H':>6}  {'R²':>6}  {'Delta':>7}",
        "  " + "-" * 48,
    ]
    for (r, th, ns, D, R2, se) in sorted(htm_hits, key=lambda x: abs(x[3]-DF_GEO_TARGET)):
        lines.append(f"  {r:5.2f}  {th:6.1f}  {ns:5d}  {D:6.3f}  {R2:6.4f}  {D-DF_GEO_TARGET:+7.3f}")

lines += [
    "",
    "-" * 70,
    "VERDICT",
    "-" * 70,
    "",
    f"  OT-37: {verdict}",
    f"  {verdict_detail}",
    "",
    "  Interpretation:",
]

if verdict == "BESTANDEN":
    lines += [
        f"  Die S³-Projektion des IFS-Fraktaltesserakts erzeugt D_H ≈ 0.77",
        f"  bei Delta-theta nahe dem HTM-Wert 58.7 deg.",
        "  Dies stuetzt die geometrische Interpretation von D_f,geo = 0.77",
        "  als Hausdorff-Dimension des projizierten 4D-Fraktals.",
        "  WICHTIG: IFS-Parameter (r, N_sub) sind nicht eindeutig bestimmt.",
        "  Ergebnis ist konsistent, aber nicht zwingend eindeutig.",
    ]
else:
    lines += [
        "  Bei Delta-theta = 58.7 deg wird D_H ≈ 0.77 nicht reproduziert.",
        "  Die Verbindung D_f,geo = 0.77 <-> Tesseract-Projektion ist",
        "  mit diesen IFS-Parametern nicht haltbar.",
        "  Pfad C (Planck-Dichte) und Pfad B (tau=10^-32 s) bleiben davon unberührt.",
    ]

lines += [
    "",
    f"  Beste Visualisierung: r={best_r}, theta={best_theta}, N_sub={best_nsub}",
    f"  -> results/OT_37_projection.png",
    f"  -> results/OT_37_boxcount.png",
    "",
    "=" * 70,
]

result_text = "\n".join(lines)
print("\n" + result_text)

out = os.path.join(RESULTS, "OT_37_result.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write(result_text)
print(f"\n  Ergebnis: {out}")
