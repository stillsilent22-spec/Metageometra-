"""
OT-15: Skaliertes SRM Halo (δ_Σ aus Torsionskompression)
==========================================================
Verbesserung über OT-14: Galaxiespezifischer Skalenradius

Skalierungsrelation (Metageometra V19):
  r_s(Galaxie) = r_s_global * (M_stellar / M_MW)^(1/3)
  r_s_global   = 42.3 kpc   (SRM-Vorhersage, kein freier Parameter)
  M_MW_stellar ≈ 5.0e10 M_sun

Einziger Freiheitsgrad: rho_0 pro Galaxie (wie in OT-14)
Vergleich: OT-14 (r_s = konstant 42.3 kpc) vs OT-15 (r_s skaliert)

Physik: Im Torsionskompressionsmodell skaliert die Shell-Dicke
        δ_Σ ∝ M^(1/3) wegen der Kugeloberflächen-Geometrie.
"""

import os, sys, io, math, urllib.request
import numpy as np
from scipy import optimize, integrate, interpolate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

# ── Konstanten ────────────────────────────────────────────────────────────────
G_KPC      = 4.3009e-6    # kpc (km/s)^2 / M_sun
RS_GLOBAL  = 42.3         # kpc — globaler SRM Skalenradius
SLOPE      = 1.205        # SRM-Exponent
UPS_D      = 0.5
UPS_B      = 0.7
VERR_FLOOR = 2.0
MIN_PTS    = 3
M_MW_STAR  = 5.0e10       # M_sun  — Milchstrasse Sternmasse

print("=" * 70)
print("  OT-15: Skaliertes SRM (r_s ∝ M^1/3) vs SPARC")
print("=" * 70)

# ── SRM Integraltabelle (einmal, für normiertes Profil) ───────────────────────
print("  Vorberechne SRM Integraltabelle...")
_x_tab = np.linspace(0, 25.0, 2501)
_J_tab = np.zeros(len(_x_tab))
for _i in range(1, len(_x_tab)):
    _J_tab[_i], _ = integrate.quad(
        lambda u: u**2 * (1.0 + u)**(-SLOPE), 0, _x_tab[_i], limit=200)
_J_interp = interpolate.CubicSpline(_x_tab, _J_tab, extrapolate=True)

def v2_srm(r_kpc, rho0, r_s):
    """V^2 des SRM-Halos fuer beliebiges r_s."""
    x = np.asarray(r_kpc) / r_s
    J = np.maximum(0.0, _J_interp(x))
    prefac = 4.0 * math.pi * G_KPC * r_s**3
    return rho0 * prefac * J / np.asarray(r_kpc)

def v_model(v2bary, v2halo):
    return np.sqrt(np.maximum(0.0, v2bary + v2halo))

def chi2_val(vobs, vmod, verr):
    return float(np.sum(((vobs - vmod) / verr) ** 2))

# ── Daten laden ───────────────────────────────────────────────────────────────
BASE_CDS   = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/"
TABLE2_CK  = os.path.join(RESULTS, "sparc_table2.dat")
TABLE1_CK  = os.path.join(RESULTS, "sparc_table1.dat")

def fetch_or_cache(url, cache_path):
    if os.path.exists(cache_path):
        print(f"  Cache: {os.path.basename(cache_path)}")
        with open(cache_path, "r", encoding="utf-8", errors="replace") as f:
            return f.read()
    print(f"  Download: {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "Python/MetageometraOT15"})
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read().decode("utf-8", errors="replace")
    with open(cache_path, "w", encoding="utf-8") as f:
        f.write(raw)
    return raw

raw_t2 = fetch_or_cache(BASE_CDS + "table2.dat", TABLE2_CK)
raw_t1 = fetch_or_cache(BASE_CDS + "table1.dat", TABLE1_CK)

# ── L3.6 aus table1 parsen (bytes 41-47, F7.3, GLsun = 10^9 L_sun) ────────────
# Format: Name(1-11) Type(13-14) Dist(16-21) eDist(23-27) fDist(29) i(31-34)
#         ei(36-39) L3.6(41-47) eL3.6(49-55) ...
# Zeilen können umbrechen (Lrecl=130) — CDS fixed-width aber softline-wrapped
l36_map = {}  # name -> L3.6 in 10^9 L_sun
for line in raw_t1.splitlines():
    if len(line) < 47:
        continue
    try:
        name = line[0:11].strip()
        if not name:
            continue
        l36 = float(line[40:47])
        if name and l36 > 0:
            l36_map[name] = l36
    except (ValueError, IndexError):
        pass

print(f"  L3.6 Werte fuer {len(l36_map)} Galaxien aus table1.dat")

# ── table2 parsen (same as OT-14) ──────────────────────────────────────────────
galaxies = {}
for line in raw_t2.splitlines():
    if len(line) < 59:
        continue
    try:
        name   = line[0:11].strip()
        dist   = float(line[12:18])
        r_kpc  = float(line[19:25])
        vobs   = float(line[26:32])
        evobs  = max(float(line[33:38]), VERR_FLOOR)
        vgas   = float(line[39:45])
        vdisk  = float(line[46:52])
        vbulge = float(line[53:59])
    except (ValueError, IndexError):
        continue
    if name not in galaxies:
        galaxies[name] = {"dist": dist, "r": [], "vobs": [],
                          "evobs": [], "vgas": [], "vdisk": [], "vbulge": []}
    galaxies[name]["r"].append(r_kpc)
    galaxies[name]["vobs"].append(vobs)
    galaxies[name]["evobs"].append(evobs)
    galaxies[name]["vgas"].append(vgas)
    galaxies[name]["vdisk"].append(vdisk)
    galaxies[name]["vbulge"].append(vbulge)

for name in galaxies:
    for key in ["r", "vobs", "evobs", "vgas", "vdisk", "vbulge"]:
        galaxies[name][key] = np.array(galaxies[name][key])

print(f"  {len(galaxies)} Galaxien in table2.dat")

# ── r_s Skalierung berechnen ──────────────────────────────────────────────────
# M_stellar = L3.6 * 1e9 * Upsilon_disk [M_sun]
# r_s(gal) = RS_GLOBAL * (M_stellar / M_MW_STAR)^(1/3)
rs_stat = {"used": 0, "fallback": 0, "min": 999, "max": 0}
def get_rs_gal(name):
    if name in l36_map:
        M_star = l36_map[name] * 1e9 * UPS_D   # M_sun
        if M_star <= 0:
            rs_stat["fallback"] += 1
            return RS_GLOBAL
        rs_gal = RS_GLOBAL * (M_star / M_MW_STAR) ** (1.0/3.0)
        rs_gal = max(0.5, min(1000.0, rs_gal))  # Sanity clip
        rs_stat["used"] += 1
        rs_stat["min"] = min(rs_stat["min"], rs_gal)
        rs_stat["max"] = max(rs_stat["max"], rs_gal)
        return rs_gal
    else:
        rs_stat["fallback"] += 1
        return RS_GLOBAL

# ── Fit-Schleife ──────────────────────────────────────────────────────────────
print("\n  Fitte Galaxien mit skaliertem SRM...", flush=True)

results = {}
for n_done, (name, data) in enumerate(galaxies.items()):
    r    = data["r"]
    vobs = data["vobs"]
    verr = data["evobs"]
    if len(r) < MIN_PTS:
        continue

    def sv2(v): return np.sign(v) * v**2
    v2bary = sv2(data["vgas"]) + UPS_D * data["vdisk"]**2 + UPS_B * data["vbulge"]**2

    r_s_gal = get_rs_gal(name)

    # OT-15 scaled SRM (1 param: rho_0)
    x_arr = r / r_s_gal
    J_arr = np.maximum(0.0, _J_interp(x_arr))
    prefac = 4.0 * math.pi * G_KPC * r_s_gal**3
    kernel = prefac * J_arr / r

    v2_need = np.maximum(0.0, vobs**2 - v2bary)
    k_pos = kernel > 0
    rho0_init = float(np.median(v2_need[k_pos] / kernel[k_pos])) if k_pos.any() else 0.0
    rho0_init = max(0.0, rho0_init)

    chi2_zero = chi2_val(vobs, v_model(v2bary, np.zeros_like(r)), verr)

    def chi2_srm15(log_r):
        rho0 = math.exp(log_r) if log_r > -30 else 0.0
        return chi2_val(vobs, v_model(v2bary, rho0 * kernel), verr)

    try:
        if rho0_init > 0:
            lr0 = math.log(rho0_init)
            res = optimize.minimize_scalar(chi2_srm15,
                bracket=(lr0 - 5, lr0, lr0 + 5), method="brent",
                options={"xtol": 1e-8, "maxiter": 500})
            rho0_best = math.exp(res.x) if res.x > -30 else 0.0
            chi2_15 = float(res.fun)
        else:
            rho0_best, chi2_15 = 0.0, chi2_zero
    except Exception:
        rho0_best, chi2_15 = rho0_init, chi2_zero

    if chi2_15 > chi2_zero or not math.isfinite(chi2_15):
        rho0_best, chi2_15 = 0.0, chi2_zero

    # Also compute OT-14 SRM (fixed r_s = 42.3 kpc) for reference
    J14 = np.maximum(0.0, _J_interp(r / RS_GLOBAL))
    k14 = 4.0 * math.pi * G_KPC * RS_GLOBAL**3 * J14 / r
    v2_ok14 = np.maximum(0.0, vobs**2 - v2bary)
    rho0_14 = float(np.median(v2_ok14[k14 > 0] / k14[k14 > 0])) if (k14 > 0).any() else 0.0
    rho0_14 = max(0.0, rho0_14)
    chi2_14_zero = chi2_zero
    def chi2_srm14(lr):
        r0 = math.exp(lr) if lr > -30 else 0.0
        return chi2_val(vobs, v_model(v2bary, r0 * k14), verr)
    try:
        if rho0_14 > 0:
            lr14 = math.log(rho0_14)
            res14 = optimize.minimize_scalar(chi2_srm14,
                bracket=(lr14-5, lr14, lr14+5), method="brent",
                options={"xtol":1e-8,"maxiter":500})
            rho0_14f = math.exp(res14.x) if res14.x > -30 else 0.0
            chi2_14 = float(res14.fun)
        else:
            rho0_14f, chi2_14 = 0.0, chi2_14_zero
    except Exception:
        rho0_14f, chi2_14 = rho0_14, chi2_14_zero
    if chi2_14 > chi2_14_zero or not math.isfinite(chi2_14):
        rho0_14f, chi2_14 = 0.0, chi2_14_zero

    results[name] = {
        "n": len(r), "r_s_gal": r_s_gal,
        "chi2_15": chi2_15, "chi2_14": chi2_14,
    }
    if (n_done + 1) % 25 == 0:
        print(f"    {n_done+1}/{len(galaxies)} fertig...", flush=True)

print(f"  {len(results)} Galaxien gefittet.")

# ── Statistiken ───────────────────────────────────────────────────────────────
names  = sorted(results.keys())
c15    = np.array([results[n]["chi2_15"] for n in names])
c14    = np.array([results[n]["chi2_14"] for n in names])
n_pts  = np.array([results[n]["n"] for n in names])
rs_arr = np.array([results[n]["r_s_gal"] for n in names])

dchi2   = c15 - c14   # >0 means OT-14 better, <0 means OT-15 better
n15_win = int((dchi2 < 0).sum())
n14_win = int((dchi2 > 0).sum())

# AIC: DAIC = AIC(15) - AIC(14) = (2*1 + chi2_15) - (2*1 + chi2_14) = dchi2
# Since BOTH have 1 free parameter, AIC = chi2 exactly (same n params)
# OT-15 only differs in having a different (physics-predicted) r_s!

rchi2_15 = c15 / np.maximum(1, n_pts - 1)
rchi2_14 = c14 / np.maximum(1, n_pts - 1)

total_15 = float(c15.sum())
total_14 = float(c14.sum())

print(f"\n  Gesamt-chi2: OT-14 = {total_14:.1f},  OT-15 = {total_15:.1f}")
print(f"  OT-15 verbessert OT-14: {n15_win}/{len(names)} Galaxien ({100*n15_win/len(names):.1f}%)")
print(f"  OT-14 besser: {n14_win}/{len(names)} Galaxien ({100*n14_win/len(names):.1f}%)")
print(f"  r_s_gal: {rs_arr.min():.2f} – {rs_arr.max():.2f} kpc (Median {np.median(rs_arr):.2f} kpc)")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 5), constrained_layout=True)
fig.suptitle("OT-15: Skaliertes SRM vs OT-14 (SPARC, n=175)", fontsize=12)

ax1 = fig.add_subplot(1, 3, 1)
ax1.hist(rs_arr, bins=25, color="steelblue", edgecolor="white")
ax1.axvline(RS_GLOBAL, color="red", ls="--", lw=1.5, label=f"OT-14: {RS_GLOBAL} kpc")
ax1.set_xlabel("r_s(gal) [kpc]", fontsize=10)
ax1.set_ylabel("Galaxien", fontsize=10)
ax1.set_title("Skalierter r_s pro Galaxie", fontsize=10)
ax1.legend(fontsize=9)

ax2 = fig.add_subplot(1, 3, 2)
ax2.hist(dchi2, bins=30, color="coral", edgecolor="white")
ax2.axvline(0, color="k", lw=1.5)
ax2.axvline(float(np.median(dchi2)), color="red", ls="--",
            label=f"Median={np.median(dchi2):.1f}")
ax2.set_xlabel("chi2_15 - chi2_14", fontsize=10)
ax2.set_ylabel("Galaxien", fontsize=10)
ax2.set_title(f"Deltachi2 (negativ = OT-15 besser)\nOT-15 besser: {n15_win}/{len(names)}", fontsize=9)
ax2.legend(fontsize=9)

ax3 = fig.add_subplot(1, 3, 3)
ax3.scatter(rs_arr, dchi2, c="steelblue", s=10, alpha=0.6)
ax3.axhline(0, color="k", lw=1)
ax3.set_xlabel("r_s(gal) [kpc]", fontsize=10)
ax3.set_ylabel("chi2_15 - chi2_14", fontsize=10)
ax3.set_title("Verbesserung vs Skalenradius", fontsize=10)

plot_path = os.path.join(RESULTS, "OT_15_plot.png")
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
print(f"  Plot: {plot_path}")

# ── Textbericht ───────────────────────────────────────────────────────────────
lines = []
lines += [
    "=" * 70,
    "OT-15: Skaliertes SRM Halo (r_s ∝ M^(1/3)) vs SPARC",
    "=" * 70, "",
    "Daten:  SPARC Lelli+2016 (CDS J/AJ/152/157)",
    f"        {len(results)} Galaxien gefittet, {int(n_pts.sum())} Datenpunkte",
    "",
    "SRM-Skalierung (Torsionskompression, V19):",
    f"  r_s(gal) = {RS_GLOBAL} kpc × (M_stellar / M_MW)^(1/3)",
    f"  M_MW_stellar = {M_MW_STAR:.1e} M_sun",
    f"  M_stellar = L_3.6 × 10^9 × {UPS_D} M_sun/L_sun",
    "",
    f"  r_s Range: {rs_arr.min():.2f} – {rs_arr.max():.2f} kpc",
    f"  r_s Median: {np.median(rs_arr):.2f} kpc",
    f"  (OT-14 hatte fest r_s = {RS_GLOBAL} kpc fuer ALLE Galaxien)",
    "",
    "-" * 70, "ERGEBNISSE", "-" * 70, "",
    f"  Gesamt-chi2  OT-14 (r_s fest):  {total_14:.1f}",
    f"  Gesamt-chi2  OT-15 (r_s skal.): {total_15:.1f}",
    f"  Delta-Gesamt: {total_15 - total_14:+.1f} ({100*(total_15-total_14)/total_14:+.1f}%)",
    "",
    f"  BEIDE Modelle haben 1 freien Parameter (nur rho_0)!",
    f"  -> Direkter Vergleich OHNE Parameterstrafe",
    "",
    f"  OT-15 verbessert OT-14: {n15_win}/{len(names)} = {100*n15_win/len(names):.1f}% der Galaxien",
    f"  OT-14 bleibt besser:    {n14_win}/{len(names)} = {100*n14_win/len(names):.1f}% der Galaxien",
    f"  Median Delta-chi2 (OT15-OT14): {float(np.median(dchi2)):+.2f}",
    f"  Median red.chi2 OT-14: {float(np.median(rchi2_14)):.3f}",
    f"  Median red.chi2 OT-15: {float(np.median(rchi2_15)):.3f}",
    "",
    "-" * 70, "INTERPRETATION", "-" * 70, "",
]

if total_15 < total_14:
    improv = 100 * (total_14 - total_15) / total_14
    lines += [
        f"  OT-15 VERBESSERT OT-14 um {improv:.1f}% im Gesamt-chi2.",
        "  Die Massenskalierung r_s ∝ M^(1/3) erhoeht die Passguete.",
        "  Das deutet auf eine physikalisch korrekte Skaling des SRM-Halos hin.",
        "  OT-15 BESTEHT (bedingt) — Massenskalierung statistisch gestützt.",
    ]
elif abs(total_15 - total_14) < 0.05 * total_14:
    lines += [
        "  OT-15 und OT-14 geben praktisch identische Ergebnisse (<5% Unterschied).",
        "  Das bedeutet: die Skalierung r_s ∝ M^(1/3) hat KEINEN messbaren Effekt.",
        "  Moeglicher Grund: SPARC-Radien << r_s auch nach Skalierung.",
        "  -> Massenskalierung auf SPARC-Daten nicht diskriminierbar.",
    ]
else:
    lines += [
        f"  OT-15 VERSCHLECHTERT OT-14 um {100*(total_15-total_14)/total_14:.1f}%.",
        f"  Die Massenskalierung r_s ∝ M^(1/3) passt schlechter als festes r_s.",
        "  Moeglicherweise ist das SRM-Skalierungsgesetz inkorrekt,",
        "  oder SPARC-Galaxien spannen nicht den noetigen Massenbereich auf.",
        "  OT-15 NEGATIV — Skalierung nicht bestaetigt.",
    ]

lines += ["", f"Plot: {plot_path}", "=" * 70]
result_text = "\n".join(lines)
print("\n" + result_text)
out = os.path.join(RESULTS, "OT_15_result.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write(result_text)
print(f"\n  Ergebnis: {out}")
