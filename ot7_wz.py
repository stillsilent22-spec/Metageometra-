"""
OT-7: w(z) HTM-Vorhersage vs DESI DR2 BAO Daten
=================================================
HTM-Formel (Metageometra V19):
  w_HTM(z) = -1 + (1 - D_f_geo/2) * D_f_eff * (1+z)^(3*D_f_eff)
           = -1 + (1 - 1.8/2) * 0.1 * (1+z)^(3*0.1)
           = -1 + 0.01 * (1+z)^0.3

Wobei aus der E-Mail:
  (1 - 1.8/2) = 0.1    [geometrische Dämpfung]
  0.1          = D_f_eff [effektive fraktale Dimension]
  3*0.1 = 0.3  [skalierender Exponent]

DESI DR2 Daten: github.com/CobayaSampler/bao_data/desi_bao_dr2
  arXiv:2503.14738 (Phys. Rev. D 112, 083515, 2025)

CPL best-fit (DESI DR2 + CMB, Tabelle V):
  w₀ = -0.838 ± 0.055,  wₐ = -0.685 ± 0.28

Kosmologische Fiducialwerte (Planck 2018 / DESI DR2):
  H₀ = 67.36 km/s/Mpc,  Ω_m = 0.3153,  Ω_DE = 0.6847
  r_s = 147.09 Mpc  (Schallhorizont)
"""

import os, math, urllib.request
import numpy as np
from scipy import integrate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS, exist_ok=True)

# ── Kosmologische Konstanten ──────────────────────────────────────────────────
H0      = 67.36          # km/s/Mpc
C_KMS   = 2.998e5        # km/s
OM      = 0.3153         # Omega_matter
ODE     = 1.0 - OM       # Omega_dark-energy (flat universe)
RS      = 147.09         # Mpc — Schallhorizont (Planck 2018)

# ── DESI DR2 BAO Daten (direkt aus CobayaSampler/bao_data) ───────────────────
# Versuche Download, fallback auf embedded Werte (arXiv:2503.14738, Tabellen 1-2)
MEAN_URL = ("https://raw.githubusercontent.com/CobayaSampler/bao_data/"
            "master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt")
COV_URL  = ("https://raw.githubusercontent.com/CobayaSampler/bao_data/"
            "master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt")

def try_download(url):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Python/metageometra-ot7"})
        with urllib.request.urlopen(req, timeout=15) as r:
            return r.read().decode("utf-8")
    except Exception as e:
        print(f"  Download fehlgeschlagen ({e}), nutze eingebettete Werte.")
        return None

def parse_mean_file(text):
    """Parst CobayaSampler mean-Datei: gibt list[(z, val, qty)] zurück."""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                z   = float(parts[0])
                val = float(parts[1])
                qty = parts[2]
                rows.append((z, val, qty))
            except ValueError:
                pass
    return rows

# Eingebettete DESI DR2 Werte (exakt aus dem mean.txt, oben bestätigt)
DESI_MEAN_EMBEDDED = """# [z] [value at z] [quantity]
0.295 7.94167639 DV_over_rs
0.510 13.58758434 DM_over_rs
0.510 21.86294686 DH_over_rs
0.706 17.35069094 DM_over_rs
0.706 19.45534918 DH_over_rs
0.934 21.57563956 DM_over_rs
0.934 17.64149464 DH_over_rs
1.321 27.60085612 DM_over_rs
1.321 14.17602155 DH_over_rs
1.484 30.51190063 DM_over_rs
1.484 12.81699964 DH_over_rs
2.330 38.98897396 DM_over_rs
2.330  8.63154567 DH_over_rs
"""

# Sigma (sqrt. diagonal der Kovarianzmatrix)
# Direkt aus der cov.txt (Diagonale, 13 Elemente)
DESI_COV_DIAG = np.array([
    5.78998687e-03,   # DV BGS
    2.83473742e-02,   # DM LRG1
    1.83928040e-01,   # DH LRG1
    3.23752442e-02,   # DM LRG2
    1.11469198e-01,   # DH LRG2
    2.61732816e-02,   # DM LRG3
    4.04183878e-02,   # DH LRG3
    1.05336516e-01,   # DM ELG2
    5.04233092e-02,   # DH ELG2
    5.83020277e-01,   # DM QSO
    2.68336193e-01,   # DH QSO
    1.02136194e-02,   # DM Lya
    2.82685779e-01,   # DH Lya
])
DESI_SIGMA = np.sqrt(DESI_COV_DIAG)

print("\n" + "="*70)
print("  OT-7: w(z) HTM vs DESI DR2 vs ΛCDM vs CPL")
print("="*70)

# Lade oder nutze embedded Daten
raw = try_download(MEAN_URL)
if raw:
    measurements = parse_mean_file(raw)
    print(f"  Daten heruntergeladen: {len(measurements)} Datenpunkte")
else:
    measurements = parse_mean_file(DESI_MEAN_EMBEDDED)
    print(f"  Eingebettete Daten: {len(measurements)} Datenpunkte")

# ── w(z)-Modelle ──────────────────────────────────────────────────────────────
D_F_GEO = 0.77
D_F_EFF = 0.1       # aus der HTM-Formel des Users  (0.01 Gesamtkoeffizient)

def w_lcdm(z):
    return -1.0

def w_cpl(z, w0=-0.838, wa=-0.685):
    """CPL parametrization: w(z) = w0 + wa * z/(1+z)"""
    return w0 + wa * z / (1.0 + z)

def w_htm(z):
    """HTM-Vorhersage: w(z) = -1 + 0.01*(1+z)^0.3"""
    return -1.0 + (1.0 - 1.8/2.0) * 0.1 * (1.0 + z)**0.3

# ── H(z) für variables w(z) ──────────────────────────────────────────────────
def integral_w(z_max, w_func):
    """∫₀^z (1+w(z'))/(1+z') dz'"""
    if z_max <= 0:
        return 0.0
    val, _ = integrate.quad(lambda z: (1.0 + w_func(z)) / (1.0 + z), 0, z_max)
    return val

def H_over_H0(z, w_func):
    matter = OM * (1.0 + z)**3
    de     = ODE * math.exp(3.0 * integral_w(z, w_func))
    return math.sqrt(matter + de)

def DM_over_rs(z, w_func):
    """DM(z)/rs = (c/H₀)/rs * ∫₀ᶻ dz'/E(z')"""
    prefac = C_KMS / (H0 * RS)
    val, _ = integrate.quad(lambda zp: 1.0 / H_over_H0(zp, w_func), 0, z)
    return prefac * val

def DH_over_rs(z, w_func):
    """DH(z)/rs = c/(H(z)·rs)"""
    return C_KMS / (H0 * RS * H_over_H0(z, w_func))

def DV_over_rs(z, w_func):
    """DV/rs = [(z·DM²·DH)^(1/3)] / rs"""
    dm = DM_over_rs(z, w_func) * RS   # in Mpc
    dh = DH_over_rs(z, w_func) * RS
    dv = (z * dm**2 * dh)**(1.0/3.0)
    return dv / RS

# ── Modell-Vorhersagen berechnen ──────────────────────────────────────────────
print("\n  Berechne Modellvorhersagen...")

models = {
    "LCDM": w_lcdm,
    "CPL":  w_cpl,
    "HTM":  w_htm,
}

def predict(meas_list, w_func):
    """Berechne Modellvorhersagen für alle Messpunkte."""
    preds = []
    for z, obs, qty in meas_list:
        if qty == "DV_over_rs":
            preds.append(DV_over_rs(z, w_func))
        elif qty == "DM_over_rs":
            preds.append(DM_over_rs(z, w_func))
        elif qty == "DH_over_rs":
            preds.append(DH_over_rs(z, w_func))
    return np.array(preds)

obs_vals = np.array([v for _, v, _ in measurements])

chi2 = {}
for name, wf in models.items():
    pred = predict(measurements, wf)
    resid = (obs_vals - pred) / DESI_SIGMA
    chi2[name] = float(np.sum(resid**2))
    print(f"  χ²({name:5s}) = {chi2[name]:.2f}  (dof=13, red.χ²={chi2[name]/13:.2f})")

print()
best = min(chi2, key=chi2.get)
for name in ["LCDM", "CPL", "HTM"]:
    delta = chi2[name] - chi2[best]
    print(f"  Δχ²({name:5s} - {best:4s}) = {delta:.2f}")

# ── Chi²-Detailtabelle ────────────────────────────────────────────────────────
print("\n  Residuen-Tabelle:")
print(f"  {'z':>6}  {'Qty':<12}  {'Obs':>8}  {'σ':>6}  "
      f"{'ΛCDM':>8}  {'CPL':>8}  {'HTM':>8}")
print("  " + "-"*68)
for i, (z, obs, qty) in enumerate(measurements):
    p_l = predict([(z, obs, qty)], w_lcdm)[0]
    p_c = predict([(z, obs, qty)], w_cpl)[0]
    p_h = predict([(z, obs, qty)], w_htm)[0]
    sig = DESI_SIGMA[i]
    print(f"  {z:>6.3f}  {qty:<12}  {obs:>8.4f}  {sig:>6.4f}  "
          f"{p_l:>8.4f}  {p_c:>8.4f}  {p_h:>8.4f}")

# ── Plot ──────────────────────────────────────────────────────────────────────
print("\n  Erstelle Plot...")
z_arr = np.linspace(0, 2.5, 300)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle("OT-7: HTM w(z) vs DESI DR2 (arXiv:2503.14738)", fontsize=13)

# Linkes Panel: w(z) Verläufe
ax1.axhline(-1, color="gray", ls="--", lw=1.2, label="ΛCDM (w = -1)")
ax1.plot(z_arr, [w_htm(z) for z in z_arr], "r-", lw=2.5, label="HTM")
ax1.plot(z_arr, [w_cpl(z) for z in z_arr], "b-", lw=2.0, label=f"CPL (w₀={-0.838}, wₐ={-0.685})")

# CPL 1σ-Band
wa_lo, wa_hi = -0.685 - 0.28, -0.685 + 0.28
w0_lo, w0_hi = -0.838 - 0.055, -0.838 + 0.055
ax1.fill_between(z_arr,
    [w_cpl(z, w0_lo, wa_lo) for z in z_arr],
    [w_cpl(z, w0_hi, wa_hi) for z in z_arr],
    alpha=0.15, color="blue", label="CPL 1σ-Band")

ax1.set_xlabel("Rotverschiebung z", fontsize=11)
ax1.set_ylabel("w(z)", fontsize=11)
ax1.set_xlim(0, 2.5)
ax1.set_ylim(-1.6, -0.5)
ax1.legend(fontsize=10)
ax1.set_title("Equation-of-State w(z)")
ax1.grid(True, alpha=0.3)

# Rechtes Panel: Chi²-Balkendiagramm + Residuen
z_meas  = [z for z, _, qty in measurements if qty in ("DM_over_rs", "DH_over_rs")]
obs_dm  = [v for z, v, qty in measurements if qty in ("DM_over_rs", "DH_over_rs")]
pred_lcdm = predict([(z, v, q) for z, v, q in measurements
                     if q in ("DM_over_rs", "DH_over_rs")], w_lcdm)
pred_htm  = predict([(z, v, q) for z, v, q in measurements
                     if q in ("DM_over_rs", "DH_over_rs")], w_htm)
sig_sub   = [DESI_SIGMA[i] for i, (z, v, q) in enumerate(measurements)
             if q in ("DM_over_rs", "DH_over_rs")]

# baue Inset: χ² Balken
names  = ["ΛCDM", "CPL", "HTM"]
chi2_v = [chi2["LCDM"], chi2["CPL"], chi2["HTM"]]
colors = ["gray", "blue", "red"]
bars = ax2.bar(names, chi2_v, color=colors, alpha=0.7, edgecolor="black")
ax2.set_ylabel("χ² (13 Messpunkte)", fontsize=11)
ax2.set_title("Modellvergleich χ² (DESI DR2 BAO)")
for bar, val in zip(bars, chi2_v):
    ax2.text(bar.get_x() + bar.get_width()/2, val + 0.2,
             f"{val:.1f}", ha="center", va="bottom", fontsize=12, fontweight="bold")
ax2.axhline(13, color="k", ls=":", lw=1, label="χ²=dof (perfekt)")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3, axis="y")

plt.tight_layout()
plot_path = os.path.join(RESULTS, "OT_07_plot.png")
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
print(f"  Plot gespeichert: {plot_path}")

# ── Textausgabe ───────────────────────────────────────────────────────────────
lines = []
lines.append("=" * 70)
lines.append("OT-07: w(z) HTM-Vorhersage vs DESI DR2 BAO")
lines.append("=" * 70)
lines.append("")
lines.append("Daten: DESI DR2 (arXiv:2503.14738, Phys. Rev. D 112, 083515, 2025)")
lines.append("       CobayaSampler/bao_data/desi_bao_dr2")
lines.append(f"       13 Messpunkte: DV/rs + DM/rs + DH/rs bei 7 Rotverschiebungen")
lines.append("")
lines.append("HTM-Formel:")
lines.append("  w_HTM(z) = -1 + (1 - 1.8/2) * 0.1 * (1+z)^(3*0.1)")
lines.append("           = -1 + 0.01 * (1+z)^0.3")
lines.append("")
lines.append("  z=0:    w = -0.990")
lines.append("  z=0.5:  w = {:.4f}".format(w_htm(0.5)))
lines.append("  z=1.0:  w = {:.4f}".format(w_htm(1.0)))
lines.append("  z=2.3:  w = {:.4f}".format(w_htm(2.3)))
lines.append("")
lines.append("CPL best-fit (DESI DR2 + CMB, Tab. V):")
lines.append("  w₀ = -0.838 ± 0.055,  wₐ = -0.685 ± 0.28")
lines.append("  Präferenz gg. ΛCDM: 3.1σ")
lines.append("")
lines.append("Chi²-Ergebnisse (13 Freiheitsgrade):")
lines.append(f"  χ²(ΛCDM) = {chi2['LCDM']:.2f}   red.χ² = {chi2['LCDM']/13:.2f}")
lines.append(f"  χ²(CPL)  = {chi2['CPL']:.2f}   red.χ² = {chi2['CPL']/13:.2f}")
lines.append(f"  χ²(HTM)  = {chi2['HTM']:.2f}   red.χ² = {chi2['HTM']/13:.2f}")
lines.append("")
lines.append(f"  Bestes Modell: {best}")
lines.append(f"  Δχ²(ΛCDM - {best}) = {chi2['LCDM'] - chi2[best]:.2f}")
lines.append(f"  Δχ²(HTM  - {best}) = {chi2['HTM']  - chi2[best]:.2f}")
lines.append(f"  Δχ²(CPL  - {best}) = {chi2['CPL']  - chi2[best]:.2f}")
lines.append("")
lines.append("INTERPRETATION:")
lines.append("  HTM w(z): sehr kleine Abweichung von ΛCDM (max ~1% bei z=2.3)")
lines.append("  Der CPL-Fit passt besser als ΛCDM (erwartet, da 3.1σ Signal)")

if chi2["HTM"] < chi2["LCDM"]:
    dchi = chi2["LCDM"] - chi2["HTM"]
    lines.append(f"  HTM verbessert ΛCDM um Δχ² = {dchi:.2f} -> schwache HTM-Signatur")
else:
    dchi = chi2["HTM"] - chi2["LCDM"]
    lines.append(f"  HTM liegt Δχ² = {dchi:.2f} schlechter als ΛCDM -> kein zusätzlicher Gewinn")

lines.append("")
lines.append("  Das DESI 3.1σ-Signal geht in Richtung w₀ > -1, wₐ < 0")
lines.append("  HTM-Vorhersage: w₀_eff ≈ -0.990 -> innerhalb 2σ von ΛCDM")
lines.append("  -> HTM-Abweichung zu klein um durch DESI DR2 aktuell zu trennen")
lines.append("")
lines.append(f"Plot: {plot_path}")
lines.append("=" * 70)

result_text = "\n".join(lines)
print(result_text)

out_path = os.path.join(RESULTS, "OT_07_result.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(result_text)
print(f"\n  Ergebnis gespeichert: {out_path}")
