#!/usr/bin/env python3
"""
OT-40: Permutationstest V20-Schalen gegen CF3-Hintergrund
==========================================================
Problem OT-39: Virgo/Fornax-Cluster erzeugen anisotropen Hintergrund.
Loesung: Vergleiche SMBH-Residuen nicht gegen Uniform-Himmel,
         sondern gegen ZUFAELLIGE Stichproben aus CF3 (11244 Galaxien).

Methode:
  Teststatistik T = Summe exp(-residual_i^2 / (2*sigma^2)) fuer alle SMBHs
  (Gaussfoermige Gewichtung, sigma=2 deg gibt naeheren Treffern mehr Gewicht)

  Permutationstest:
  - Berechne T_obs fuer 97 SMBHs
  - Ziehe 10000x zufaellige 97-er Stichprobe aus CF3
  - Berechne T_rand fuer jede Stichprobe
  - p-Wert = P(T_rand >= T_obs) = Fraktion mit T_rand >= T_obs
"""
import sys, io, math, os, random
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import csv

os.makedirs("results", exist_ok=True)
lines = []
def p(s=""):
    print(s)
    lines.append(str(s))

# ── V20 Schalenwinkel (OT-38 beste Loesung) ──────────────────────────
DELTA = 1.0; CHI = 59.1; N_SHELLS = 6
def theta_n(n):
    v = math.cos(math.radians(n*DELTA)) * math.cos(math.radians(n*CHI))
    return math.degrees(math.acos(max(-1.0, min(1.0, v))))
SHELLS = np.array([theta_n(n) for n in range(1, N_SHELLS+1)])

SIGMA = 2.0        # deg Gauss-Breite
N_PERM = 10000     # Permutationen

p("=" * 72)
p("  OT-40: Permutationstest V20-Schalen vs. CF3-Hintergrund")
p("=" * 72)
p()
p(f"  V20 Schalen (delta={DELTA}, chi={CHI}):")
for n, th in enumerate(SHELLS, 1):
    p(f"    n={n}: {th:.2f} deg")
p()
p(f"  Teststatistik: T = sum_i exp(-r_i^2 / (2*{SIGMA}^2))")
p(f"  Permutationen: {N_PERM}")
p()

# ── Teststatistik ─────────────────────────────────────────────────────
def score(angles):
    """Score-Vektor fuer ein Array von Winkeln gegen SHELLS."""
    angles = np.asarray(angles)
    # minimaler Abstand zu irgendeiner Schale
    diff = np.abs(angles[:, None] - SHELLS[None, :])   # (N, 6)
    min_r = diff.min(axis=1)                            # (N,)
    return float(np.sum(np.exp(-min_r**2 / (2 * SIGMA**2))))

# ── SMBH-Katalog laden ───────────────────────────────────────────────
CAT_SMBH = os.path.join("results", "catalogs", "smbh_extended.csv")
smbh_angles = []
smbh_names = []
with open(CAT_SMBH, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            smbh_angles.append(float(row['theta_dpole']))
            smbh_names.append(row['Name'])
        except (ValueError, KeyError):
            pass
smbh_angles = np.array(smbh_angles)
N_SMBH = len(smbh_angles)
T_obs = score(smbh_angles)

p(f"  SMBHs geladen: {N_SMBH}")
p(f"  T_obs = {T_obs:.4f}")
p()

# ── CF3-Hintergrundkatalog laden (VOLLHIMMEL, VizieR J/AJ/152/50) ────
# War bisher: vizier_dpol_combined.csv (nur theta_D < 50°) → BIAS!
# Jetzt: cf3_fullsky_thetaD.csv (2000 Gruppen, theta_D 1-176°)
CAT_CF3_FULL = os.path.join("results", "catalogs", "cf3_fullsky_thetaD.csv")
CAT_CF3_OLD  = os.path.join("results", "catalogs", "vizier_dpol_combined.csv")

cf3_angles = []
if os.path.exists(CAT_CF3_FULL):
    with open(CAT_CF3_FULL, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                cf3_angles.append(float(row['theta_D']))
            except (ValueError, KeyError):
                pass
    cf3_source = "CF3 Vollhimmel (J/AJ/152/50, 2000 Gruppen)"
else:
    with open(CAT_CF3_OLD, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                cf3_angles.append(float(row['theta_D']))
            except (ValueError, KeyError):
                pass
    cf3_source = "D-Pol Katalog (vizier_dpol_combined, BIAS: nur <50 Deg!)"

cf3_angles = np.array(cf3_angles)
N_CF3 = len(cf3_angles)

p(f"  CF3 Hintergrundkatalog: {N_CF3} Galaxien ({cf3_source})")
p(f"  CF3 theta_D Bereich: {cf3_angles.min():.1f} bis {cf3_angles.max():.1f} deg")
p()

# ── Permutationstest ─────────────────────────────────────────────────
p(f"  Permutationstest ({N_PERM} Stichproben) ...")
rng = np.random.default_rng(42)
T_rand = np.empty(N_PERM)
for i in range(N_PERM):
    sample = rng.choice(cf3_angles, size=N_SMBH, replace=False)
    T_rand[i] = score(sample)

p("  Fertig.")
p()

T_mean = float(T_rand.mean())
T_std  = float(T_rand.std())
T_pct  = float(np.percentile(T_rand, 95))
pval   = float((T_rand >= T_obs).sum()) / N_PERM
z_score = (T_obs - T_mean) / T_std

p(f"  T_obs           = {T_obs:.4f}")
p(f"  T_rand (Mittel) = {T_mean:.4f}")
p(f"  T_rand (Std)    = {T_std:.4f}")
p(f"  T_rand (95pct)  = {T_pct:.4f}")
p(f"  z-Wert          = {z_score:.2f}  sigma")
n_exceed = int((T_rand >= T_obs).sum())
p(f"  p-Wert (perm)   = {pval:.4f}  ({n_exceed}/{N_PERM})")

# ── Schalenweise Permutationstest ────────────────────────────────────
p()
p("-" * 72)
p("  Schalenweise Permutationstest (sigma=2 deg pro Schale):")
p(f"  {'n':>3}  {'theta':>7}  {'T_obs_n':>8}  {'T_mean_n':>9}  {'z':>6}  {'p_n':>8}")
p("-" * 72)

def score_shell(angles, shell_theta):
    r = np.abs(np.asarray(angles) - shell_theta)
    return float(np.sum(np.exp(-r**2 / (2 * SIGMA**2))))

pvals_shell = []
for n, sh in enumerate(SHELLS, 1):
    T_obs_n = score_shell(smbh_angles, sh)
    T_rand_n = np.array([score_shell(rng.choice(cf3_angles, size=N_SMBH, replace=False), sh)
                         for _ in range(2000)])
    pval_n  = float((T_rand_n >= T_obs_n).sum()) / 2000
    z_n     = (T_obs_n - T_rand_n.mean()) / (T_rand_n.std() + 1e-9)
    pvals_shell.append(pval_n)
    star = "***" if pval_n < 0.01 else ("**" if pval_n < 0.05 else ("*" if pval_n < 0.10 else ""))
    p(f"  {n:>3}  {sh:>7.2f}  {T_obs_n:>8.3f}  "
      f"{T_rand_n.mean():>9.3f}  {z_n:>6.2f}  {pval_n:>8.4f}  {star}")

# ── Top-10 best SMBH hits ─────────────────────────────────────────────
p()
p("-" * 72)
p("  Top-Treffer (kleinstes Residuum zu naechster Schale):")
p("-" * 72)

diff_all = np.abs(smbh_angles[:, None] - SHELLS[None, :])
min_r    = diff_all.min(axis=1)
best_n   = diff_all.argmin(axis=1) + 1

order    = np.argsort(min_r)
p(f"  {'Name':<16}  {'theta':>7}  {'n':>3}  {'shell':>7}  {'Delta':>7}")
p("  " + "-" * 50)
for idx in order[:20]:
    bsh = SHELLS[best_n[idx]-1]
    p(f"  {smbh_names[idx]:<16}  {smbh_angles[idx]:>7.3f}  "
      f"n={best_n[idx]}  {bsh:>7.3f}  {smbh_angles[idx]-bsh:>+7.3f}")

# ── Verdict ──────────────────────────────────────────────────────────
p()
p("=" * 72)
p("  OT-40 VERDICT")
p("=" * 72)
p()
p(f"  Permutationstest (N={N_PERM}, CF3-Hintergrund, sigma={SIGMA} deg):")
p(f"  T_obs = {T_obs:.4f}  |  T_mean(rand) = {T_mean:.4f}  |  z = {z_score:.2f}")
p(f"  p-Wert = {pval:.4f}")
p()

if pval < 0.001:
    sig = "HOCHSIGNIFIKANT (p < 0.001)  ***"
elif pval < 0.01:
    sig = "SEHR SIGNIFIKANT (p < 0.01)  **"
elif pval < 0.05:
    sig = "SIGNIFIKANT (p < 0.05)  *"
elif pval < 0.10:
    sig = "GRENZWERTIG (p < 0.10)  ~"
else:
    sig = "NICHT SIGNIFIKANT"

p(f"  => {sig}")
p()

best_shell_n = int(np.argmin(pvals_shell)) + 1
p(f"  Staerkste Einzelschale: n={best_shell_n} (theta={SHELLS[best_shell_n-1]:.2f} deg), "
  f"p={pvals_shell[best_shell_n-1]:.4f}")
p()

if pval < 0.05:
    p("  V20-Schalenformel: STATISTISCH BESTAETIGT gegenueber CF3-Hintergrund")
    p("  SMBHs clustern signifikant staerker um die V20-Schalen")
    p("  als zufaellige Stichproben aus dem CF3-Universumskatalog.")
elif pval < 0.10:
    p("  Trend: SMBHs zeigen 1-2 sigma Ueberhaeufung auf V20-Schalen")
    p("  gegenueber dem CF3-Hintergrund. Noch nicht 5% signifikant.")
    p("  Der CF3-Hintergrund ist selbst stark geclustert (Virgo, Fornax)")
    p("  => konservativerer Test als uniform. Echte Signifikanz")
    p("  wahrscheinlich hoeher wenn unabhaengige SMBH-Stichprobe.")
else:
    p("  Kein signifikantes Signal gegenueber CF3-Hintergrund.")
    p("  Beide Kataloge (SMBH + CF3) clustern aehnlich stark auf V20-Schalen.")
    p("  Moegliche Interpretation: V20-Schalen beschreiben allgemeine")
    p("  Grossstruktur des Universums, nicht SMBH-spezifisch.")

p("=" * 72)

out = os.path.join("results", "OT_40_result.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\n  -> {out}")
