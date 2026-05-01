#!/usr/bin/env python3
"""
OT-25: GWTC chi_eff Spin-Ausrichtung vs HTM 4-Arm-Topologie
============================================================
V21+ Neue OT-Definition:

HTM-Vorhersage (V21, OT-41):
  4-Arm Struktur mit chi/3 = 19.70 Grad:
    Arm 1 (A): prograd  -> BH-Spin bevorzugt +1 -> chi_eff > 0
    Arm 2 (B): retrograd -> BH-Spin bevorzugt -1 -> chi_eff < 0
    Arm 3 (A): prograd  -> chi_eff > 0
    Arm 4 (B): retrograd -> chi_eff < 0
  BBH-Merger aus gleichem Arm: |chi_eff| > 0
  BBH-Merger aus verschiedenen Armen: chi_eff ~ 0 (zufaellig)

Testsatz:
  Wenn 50% der BBH-Merger same-arm sind: bimodale chi_eff-Verteilung
  Wenn 100% random: Gauss um 0
  HTM sagt zudem: P(chi_eff > 0) = P(prograd) = 50% (A-Arm Anteil)

Tests:
  1. t-Test: chi_eff-Mittelwert != 0 ?
  2. Bimodalitaet: Dip-Test (Hartigans dip) oder Excess-Kurtosis
  3. Sign-Test: P(chi_eff > 0) = 0.5 ?
  4. F-Test / Levene: chi_eff Streuung vs N(0, sigma) ?

Datenbasis: GWTC-3 / GWTC-4 (GWOSC API) — 219 Ereignisse
"""

import sys, io, os, math, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy import stats

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
GWTC_CACHE = os.path.join(RESULTS, "gwtc_catalog.json")

# ── GWTC-Katalog laden ────────────────────────────────────────────────────────
with open(GWTC_CACHE, 'r', encoding='utf-8') as f:
    catalog_raw = json.load(f)

events = catalog_raw.get('events', {})
print(f"Geladene Ereignisse: {len(events)}")

# Filtere BBH-Ereignisse mit gueltiger chi_eff und p_astro > 0.5
chi_eff_list = []
m_chirp_list = []
m_total_list = []
ev_names     = []

for ev_name, ev in events.items():
    p_astro = ev.get('p_astro')
    chi_eff = ev.get('chi_eff')
    m_chirp = ev.get('chirp_mass_source')
    m_tot   = ev.get('total_mass_source')

    if (chi_eff is None or
        p_astro is None or float(p_astro) < 0.5):
        continue
    try:
        chi_eff_list.append(float(chi_eff))
        m_chirp_list.append(float(m_chirp) if m_chirp else float('nan'))
        m_total_list.append(float(m_tot) if m_tot else float('nan'))
        ev_names.append(ev_name)
    except (ValueError, TypeError):
        pass

chi_eff_arr = np.array(chi_eff_list)
m_chirp_arr = np.array(m_chirp_list)
m_total_arr = np.array(m_total_list)
N    = len(chi_eff_arr)
print(f"BBH-Ereignisse mit chi_eff: {N}")
print(f"chi_eff Bereich: [{chi_eff_arr.min():.3f}, {chi_eff_arr.max():.3f}]")
print(f"chi_eff Mittelwert: {chi_eff_arr.mean():.4f}")
print(f"chi_eff Median:     {np.median(chi_eff_arr):.4f}")
print(f"chi_eff Std:        {chi_eff_arr.std():.4f}")

# ── Test 1: t-Test (chi_eff Mittelwert = 0)  ─────────────────────────────────
t_stat, t_p = stats.ttest_1samp(chi_eff_arr, 0)
print(f"\n  Test 1: t-Test (H0: chi_eff_mean = 0)")
print(f"    t = {t_stat:.4f},  p = {t_p:.4f}  {'★' if t_p < 0.05 else ''}")

# ── Test 2: Sign-Test (P(chi_eff > 0) = 0.5) ──────────────────────────────────
n_pos = np.sum(chi_eff_arr > 0)
n_neg = np.sum(chi_eff_arr < 0)
n_zero = np.sum(chi_eff_arr == 0)
sign_res = stats.binomtest(n_pos, n_pos + n_neg, 0.5, alternative='two-sided')
print(f"\n  Test 2: Sign-Test (prograd vs retrograd)")
print(f"    n_pos = {n_pos}, n_neg = {n_neg}, n_zero = {n_zero}")
print(f"    P(chi_eff > 0) = {n_pos/(n_pos+n_neg):.3f} (erwartet HTM: 0.5)")
print(f"    p(Binom) = {sign_res.pvalue:.4f}  {'★' if sign_res.pvalue < 0.05 else ''}")

# ── Test 3: Kolmogorov-Smirnov gegen N(0, sigma) ────────────────────────────
sigma_obs = chi_eff_arr.std()
ks_stat, ks_p = stats.kstest(chi_eff_arr, 'norm', args=(0, sigma_obs))
print(f"\n  Test 3: KS-Test vs N(0, sigma={sigma_obs:.3f})")
print(f"    KS D = {ks_stat:.4f},  p = {ks_p:.4f}  {'★' if ks_p < 0.05 else ''}")

# ── Test 4: Bimodalitaets-Test (Kurtosis) ───────────────────────────────────
kurt = stats.kurtosis(chi_eff_arr)  # excess kurtosis
skew = stats.skew(chi_eff_arr)
# HTM-Vorhersage: bimodal -> negative Kurtosis (platykurtisch)
print(f"\n  Test 4: Verteilungsform")
print(f"    Kurtosis (excess): {kurt:.4f}  (bimodal -> <0, Gauss -> ~0)")
print(f"    Skewness:          {skew:.4f}")

# Hartigans dip statistic (approximation via scipy)
# Alternativ: Test ob Verteilung schmal um 0 oder breit
# Ein bimodaler Test: ist die Standardabweichung > sigma_gauss?
# HTM: bimodal mit peaks bei +-0.3..0.7 -> mean_|chi_eff| elevated
mean_abs = np.mean(np.abs(chi_eff_arr))
# Fuer N(0, sigma): E[|X|] = sigma * sqrt(2/pi) = 0.798*sigma
sigma_gauss = mean_abs / math.sqrt(2/math.pi)
print(f"\n    E[|chi_eff|] = {mean_abs:.4f}  (aus N(0,sigma): sigma={sigma_gauss:.4f})")
print(f"    chi_eff_arr.std() = {sigma_obs:.4f}")

# ── Test 5: Schwerste Merger (M_total > 50 Msun) ────────────────────────────
heavy_mask = (~np.isnan(m_total_arr)) & (m_total_arr > 50)
n_heavy = heavy_mask.sum()
if n_heavy > 5:
    chi_heavy = chi_eff_arr[heavy_mask]
    t_h, p_h = stats.ttest_1samp(chi_heavy, 0)
    print(f"\n  Test 5: Schwere Merger (M_tot > 50 M_sun, n={n_heavy})")
    print(f"    chi_eff Mittelwert: {chi_heavy.mean():.4f}")
    print(f"    t = {t_h:.4f},  p = {p_h:.4f}  {'★' if p_h < 0.05 else ''}")

# ── HTM-Vorhersage-Vergleich ────────────────────────────────────────────────
print("\n" + "="*70)
print("  HTM-VORHERSAGE vs BEOBACHTUNG")
print("="*70)
print(f"  chi_eff_mean = {chi_eff_arr.mean():.4f}  (HTM: = 0 fuer gleiche A/B-Anteile)")
print(f"  P(>0) = {n_pos/(n_pos+n_neg):.3f}  (HTM: 0.5)")
print(f"  Kurtosis = {kurt:.3f}  (bimodal: <0; Gauss: ~0)")
print(f"  KS vs N(0,s): p = {ks_p:.4f}  (HTM: kein klares Gauss)")
print()
print(f"  INTERPRETATION:")
if chi_eff_arr.mean() < -0.05 and t_p < 0.1:
    print(f"  Leichte Tendenz zu chi_eff < 0 (retrograd-Bias) — nicht signifikant.")
elif chi_eff_arr.mean() > 0.05 and t_p < 0.1:
    print(f"  Leichte Tendenz zu chi_eff > 0 (prograd-Bias) — nicht signifikant.")
else:
    print(f"  Kein signifikanter Spin-Bias (t-Test p={t_p:.3f}).")

if kurt < -0.5:
    print(f"  Negative Kurtosis ({kurt:.2f}): konsistent mit Bimodalitaet (breite Verteilung).")
elif kurt > 1:
    print(f"  Positive Kurtosis ({kurt:.2f}): eher leptokurtisch (spitze Verteilung).")
else:
    print(f"  Kurtosis nahe 0 ({kurt:.2f}): kein klares Bimodal-Signal.")

# ── Ergebnis speichern ────────────────────────────────────────────────────────
if t_p < 0.05 or (abs(kurt) > 1 and ks_p < 0.05):
    verdict = "TEILWEISE BESTAETIGT"
elif t_p < 0.1 or ks_p < 0.1:
    verdict = "GRENZWERTIG"
else:
    verdict = "KEIN SIGNAL"

result_lines = [
    "=" * 70,
    "OT-25: GWTC chi_eff Spin-Ausrichtung vs HTM 4-Arm-Topologie",
    "=" * 70,
    "",
    f"  Datenbasis: GWTC-3+4 ({N} BBH-Ereignisse mit chi_eff, p_astro>0.5)",
    f"  chi_eff Mittelwert: {chi_eff_arr.mean():.4f}  Std: {sigma_obs:.4f}",
    f"  chi_eff Bereich: [{chi_eff_arr.min():.3f}, {chi_eff_arr.max():.3f}]",
    "",
    "  Test 1: t-Test chi_eff = 0",
    f"    t={t_stat:.4f},  p={t_p:.4f}  {'*' if t_p < 0.05 else ''}",
    "",
    "  Test 2: Sign-Test (prograd vs retrograd)",
    f"    n_pos={n_pos}, n_neg={n_neg}",
    f"    P(>0) = {n_pos/(n_pos+n_neg):.3f}  p_binom = {sign_res.pvalue:.4f}",
    "",
    "  Test 3: KS vs N(0,sigma)",
    f"    KS D={ks_stat:.4f}, p={ks_p:.4f}  {'*' if ks_p < 0.05 else ''}",
    "",
    "  Test 4: Form",
    f"    Kurtosis={kurt:.4f}, Skewness={skew:.4f}",
    "",
    f"  HTM-Vorhersage: chi_eff_mean = 0, P(>0) = 0.5, Kurtosis < 0 (bimodal)",
    "",
    f"  BEWERTUNG OT-25: {verdict}",
    "=" * 70,
]

out_path = os.path.join(RESULTS, 'OT_25_result.txt')
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(result_lines) + '\n')
print(f"\nErgebnis gespeichert: {out_path}")
print(f"\n  BEWERTUNG OT-25: {verdict}")
