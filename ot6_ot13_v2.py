#!/usr/bin/env python3
"""
OT-6 v2 + OT-13 v2: SMBH Shell-Analyse mit 316-Objekt-Katalog (Thomas+2016 + smbh_extended)
=============================================================================
V21.3+ Update: Nutzt kombinierten Katalog smbh_large.csv (316 Objekte)
anstelle des bisherigen 97-Objekt-Katalogs.

OT-6:  D-Pol Clustering-Test (Binom + Anderson-Darling)
OT-13: SRM Schalenüberschuss-Spektrum
"""
import sys, io, os, math, csv
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy import stats

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
CATALOGS = os.path.join(RESULTS, "catalogs")
os.makedirs(RESULTS, exist_ok=True)

# ── Framework-Konstanten ──────────────────────────────────────────────────────
THETA_0   = 58.65       # deg
D_POLE_L  = 305.0
D_POLE_B  = 25.0
SHELL_TOL = 5.0         # deg
SHELLS    = [THETA_0 * n for n in range(1, 7)]

def delta_shell(theta):
    return min(abs(theta - s) for s in SHELLS)

# ── Katalog laden ─────────────────────────────────────────────────────────────
cat_path = os.path.join(CATALOGS, 'smbh_large.csv')
objects = []
with open(cat_path, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            td = float(row['theta_dpole'])
            ta = float(row['theta_apole'])
            logM = float(row['logMbh'])
            objects.append({'Name': row['Name'], 'theta_dpole': td, 'theta_apole': ta,
                            'logMbh': logM, 'source': row.get('source','')})
        except (ValueError, KeyError):
            pass
N = len(objects)
print(f"Katalog geladen: {N} Objekte aus smbh_large.csv")

# ── Null-Verteilung (sphärisch gleichmäßig) ───────────────────────────────────
np.random.seed(0)
n_null = 500000
null_theta = np.degrees(np.arccos(np.random.uniform(-1, 1, n_null)))
null_delta = np.array([delta_shell(th) for th in null_theta])
null_near  = np.sum(null_delta < SHELL_TOL) / n_null  # P(delta < 5°)

# ── Berechne Shell-Residuen für alle Objekte ──────────────────────────────────
thetas_D = np.array([o['theta_dpole'] for o in objects])
delta_D  = np.array([delta_shell(th) for th in thetas_D])
# Dualitätssphäre: hit wenn D-Pol ODER A-Pol nahe einer Schale
delta_A  = np.array([delta_shell(180.0 - th) for th in thetas_D])
delta_dual = np.minimum(delta_D, delta_A)

hits_D    = np.sum(delta_D    < SHELL_TOL)
hits_dual = np.sum(delta_dual < SHELL_TOL)

# Erwartete Coverage (analytisch):
# Einfache D-Pol Coverage: 6 Schalen × 10° / 180° = 33.3% (flat in theta)
# Aber sphärische Gewichtung macht es theta-abhängig — nutze MC für genauen Wert
p_exp_D    = null_near
p_exp_dual = np.sum(np.minimum(null_delta, np.array([delta_shell(180.-th)
                    for th in null_theta])) < SHELL_TOL) / n_null

print(f"\n  Erwartete Hit-Rate (gleichmäßig auf Kugel):")
print(f"    D-Pol:       {p_exp_D:.4f} ({100*p_exp_D:.1f}%)")
print(f"    Dual-Sphäre: {p_exp_dual:.4f} ({100*p_exp_dual:.1f}%)")

# ── OT-6: Binom-Test ──────────────────────────────────────────────────────────
print("\n" + "="*70)
print(f"  OT-6 v2: D-Pol SMBH Shell-Clustering (N={N})")
print("="*70)

b_D_ex  = stats.binomtest(hits_D, N, p_exp_D, alternative='greater')
b_D_2s  = stats.binomtest(hits_D, N, p_exp_D, alternative='two-sided')
b_DU_ex = stats.binomtest(hits_dual, N, p_exp_dual, alternative='greater')
b_DU_2s = stats.binomtest(hits_dual, N, p_exp_dual, alternative='two-sided')

print(f"\n  A) D-Pol Standard (±{SHELL_TOL}°):")
print(f"     Treffer: {hits_D}/{N} = {hits_D/N:.1%}  (erwartet: {p_exp_D:.1%})")
print(f"     p(Überschuss): {b_D_ex.pvalue:.4f}  {'★' if b_D_ex.pvalue < 0.05 else ''}")
print(f"     p(zweiseitig): {b_D_2s.pvalue:.4f}")

print(f"\n  B) Dualitätssphäre (D-Pol ∪ A-Pol):")
print(f"     Treffer: {hits_dual}/{N} = {hits_dual/N:.1%}  (erwartet: {p_exp_dual:.1%})")
print(f"     p(Überschuss): {b_DU_ex.pvalue:.4f}  {'★' if b_DU_ex.pvalue < 0.05 else ''}")
print(f"     p(zweiseitig): {b_DU_2s.pvalue:.4f}")

# KS-Test
ks_D_stat, ks_D_p = stats.ks_2samp(delta_D, null_delta[:len(delta_D)*100])
ks_DU_stat, ks_DU_p = stats.ks_2samp(delta_dual, null_delta[:len(delta_dual)*100])
print(f"\n  KS-Test:")
print(f"    D-Pol:       D={ks_D_stat:.4f},  p={ks_D_p:.4f}  {'★' if ks_D_p < 0.05 else ''}")
print(f"    Dual-Sphäre: D={ks_DU_stat:.4f},  p={ks_DU_p:.4f}  {'★' if ks_DU_p < 0.05 else ''}")

# ── OT-13: SRM Schalenspektrum ────────────────────────────────────────────────
print("\n" + "="*70)
print(f"  OT-13 v2: SRM Schalenüberschuss-Spektrum (N={N})")
print("="*70)
print(f"  Erwartet: SRM ∝ n^(-1/3) mit positivem Gesamtüberschuss")
print()

shell_z_scores = []
for n in range(1, 7):
    center = THETA_0 * n
    tol_range = SHELL_TOL
    hits_n = np.sum(np.abs(thetas_D - center) < tol_range)
    # Expected based on spherical null
    null_hits_n = np.sum(np.abs(null_theta - center) < tol_range) / n_null * N
    if null_hits_n > 0:
        sigma_n = math.sqrt(N * (null_hits_n/N) * (1 - null_hits_n/N))
        z_n = (hits_n - null_hits_n) / sigma_n
    else:
        sigma_n = 0; z_n = 0
    star = '★' if abs(z_n) > 2 else ''
    print(f"    n={n}: center={center:.2f}°  hits={hits_n:4d}  exp={null_hits_n:5.1f}  "
          f"z={z_n:+5.2f}σ  {star}")
    shell_z_scores.append(z_n)

z_arr = np.array(shell_z_scores)
z_mean = np.mean(z_arr)
t_stat, t_p = stats.ttest_1samp(z_arr, 0)
print(f"\n  Gesamt:  z_mean={z_mean:+.3f}σ  t={t_stat:.3f}  p(t-Test)={t_p:.4f}")

# Expected SRM spectrum slope: z(n) ∝ n^{-1/3}
srm_pred = np.array([1 / (n**(1/3)) for n in range(1, 7)])
srm_pred /= srm_pred[0]  # normalize to n=1
srm_pred *= z_arr[0] if abs(z_arr[0]) > 0.5 else 1.0
corr, p_corr = stats.pearsonr(shell_z_scores, srm_pred)
print(f"  SRM n^(-1/3) Korrelation: r={corr:.3f}, p={p_corr:.4f}")

# Massengradienten-Test: logM vs theta_D
logMs = np.array([o['logMbh'] for o in objects])
r_grad, p_grad = stats.pearsonr(thetas_D, logMs)
print(f"\n  Massengradient logM vs theta_D: r={r_grad:.3f}, p={p_grad:.4f}  "
      f"{'★' if p_grad < 0.05 else ''}")
print(f"  (Vorhersage: r < 0, massivere SMBHs näher am D-Pol)")

# ── Ergebnistext schreiben ────────────────────────────────────────────────────
ot6_verdict = "BESTÄTIGT" if (b_D_ex.pvalue < 0.05 or ks_D_p < 0.05) else "KEIN SIGNAL"
ot13_verdict = ("BESTÄTIGT" if (z_mean > 0 and t_p < 0.05)
                else "INCONCLUSIVE" if abs(z_mean) > 1
                else "KEIN SIGNAL")

result_lines = [
    "=" * 70,
    f"OT-6 v2 + OT-13 v2: SMBH Shell-Analyse ({N} Objekte, Thomas+2016+smbh_ext)",
    "=" * 70,
    "",
    f"Katalog: {N} SMBHs (Thomas+2016 main + smbh_extended.csv)",
    f"Schalen: n=1..6, theta_n = n × {THETA_0}°",
    f"D-Pol:   l={D_POLE_L}°, b={D_POLE_B}°",
    "",
    "─" * 70,
    "OT-6 v2: D-Pol Clustering",
    "─" * 70,
    f"  A) D-Pol Standard:   {hits_D}/{N}={hits_D/N:.1%}  exp={p_exp_D:.1%}  "
    f"p(>)={b_D_ex.pvalue:.4f}  {'★' if b_D_ex.pvalue < 0.05 else ''}",
    f"  B) Dualitätssphäre:  {hits_dual}/{N}={hits_dual/N:.1%}  exp={p_exp_dual:.1%}  "
    f"p(>)={b_DU_ex.pvalue:.4f}  {'★' if b_DU_ex.pvalue < 0.05 else ''}",
    f"  KS D-Pol:  D={ks_D_stat:.4f}  p={ks_D_p:.4f}  {'★' if ks_D_p < 0.05 else ''}",
    f"  KS Dual:   D={ks_DU_stat:.4f}  p={ks_DU_p:.4f}  {'★' if ks_DU_p < 0.05 else ''}",
    "",
    f"  BEWERTUNG OT-6: {ot6_verdict}",
    "",
    "─" * 70,
    "OT-13 v2: SRM Schalenspektrum",
    "─" * 70,
] + [
    f"  n={n}: theta={THETA_0*n:.1f}° z={z_arr[n-1]:+.2f}σ"
    for n in range(1,7)
] + [
    "",
    f"  z_mean = {z_mean:+.3f}sig  (t-Test p={t_p:.4f}{'  *' if t_p < 0.05 else ''})",
    f"  SRM-Spektrum-Korrelation: r={corr:.3f} p={p_corr:.4f}",
    f"  Massengradient r={r_grad:.3f} p={p_grad:.4f}",
    "",
    f"  BEWERTUNG OT-13: {ot13_verdict}",
    "",
    "=" * 70,
]

result_text = '\n'.join(result_lines)
print('\n' + result_text)

out_path = os.path.join(RESULTS, 'OT_06_13_v2_result.txt')
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(result_text + '\n')
print(f'\nErgebnis gespeichert: {out_path}')
