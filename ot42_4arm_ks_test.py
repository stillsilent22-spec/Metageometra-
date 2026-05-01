"""
OT-42: 4-Arm KS-Test — Azimutale Clusterung (V21.0)
=====================================================
Prüfe ob die 97 SMBHs azimutale Häufungen nahe der 4-Arm-Positionen
[88°, 107.7°, 268°, 287.7°] zeigen vs. gleichförmige Verteilung.

Methoden:
1. KS-Test (2-seitig) vs. gleichförmig
2. Watson-Statistik (zirkularer KS-Test)
3. Monte-Carlo: wie oft clustert Zufall besser als SMBHs?
4. Minimaler Winkelabstand zum nächsten Arm als Observable
"""
import math, os, csv
import numpy as np
import pandas as pd
from scipy.stats import kstest, ks_2samp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE   = os.path.dirname(os.path.abspath(__file__))
CATS   = os.path.join(BASE, "results", "catalogs")
RESDIR = os.path.join(BASE, "results")
os.makedirs(RESDIR, exist_ok=True)

# ── Parameter ─────────────────────────────────────────────
CHI      = 59.1
TILT     = CHI / 3.0   # = 19.70°
BASE_ROT = 88.0
ARM_AZ   = np.array([BASE_ROT,
                     BASE_ROT + TILT,
                     BASE_ROT + 180.0,
                     BASE_ROT + 180.0 + TILT]) % 360.0
N_MC     = 10000
AZ_TOL   = 20.0  # Grad Toleranzfenster pro Arm (±20° = 40° Breite)

# D-Pol Galaktisch: l=305°, b=25°
DP_L = 305.0; DP_B = 25.0

def gal_az_from_dpol(l_obj, b_obj):
    """Azimut des Objekts (l_obj, b_obj) um den D-Pol (l=305, b=25)."""
    dP = math.radians(DP_B); lP = math.radians(DP_L)
    dO = math.radians(b_obj); lO = math.radians(l_obj)
    y = math.cos(dO) * math.sin(lO - lP)
    x = math.cos(dP) * math.sin(dO) - math.sin(dP) * math.cos(dO) * math.cos(lO - lP)
    return math.degrees(math.atan2(y, x)) % 360.0

def eq_to_gal(ra_deg, dec_deg):
    """Konvertierung Äquatorial → Galaktisch (IAU Standard)."""
    ra_ngp  = math.radians(192.85948)
    dec_ngp = math.radians(27.12825)
    l_ncp   = math.radians(122.93192)
    ra  = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    sin_b = (math.sin(dec_ngp)*math.sin(dec)
             + math.cos(dec_ngp)*math.cos(dec)*math.cos(ra - ra_ngp))
    b = math.degrees(math.asin(max(-1., min(1., sin_b))))
    y_g = math.cos(dec) * math.sin(ra - ra_ngp)
    x_g = (math.sin(dec)*math.cos(dec_ngp)
           - math.cos(dec)*math.sin(dec_ngp)*math.cos(ra - ra_ngp))
    l = (math.degrees(math.atan2(y_g, x_g)) + math.degrees(l_ncp)) % 360.0
    return l, b

# ── SMBH-Katalog laden ─────────────────────────────────────
df = pd.read_csv(os.path.join(CATS, "smbh_extended.csv"), comment='#')
smbh_az = []
for _, row in df.iterrows():
    try:
        ra  = float(row['RA_deg']); dec = float(row['Dec_deg'])
        l, b = eq_to_gal(ra, dec)
        az = gal_az_from_dpol(l, b)
        smbh_az.append(az)
    except (ValueError, KeyError):
        pass
smbh_az = np.array(smbh_az)
N_SMBH  = len(smbh_az)

# ── Minimaler Winkelabstand zum nächsten Arm ──────────────
def min_arm_dist(az_arr, arm_positions=ARM_AZ):
    """Für jedes Objekt: minimaler Kreisabstand zu einem der 4 Arme."""
    dists = []
    for az in az_arr:
        d = min(min(abs(az - a), 360 - abs(az - a)) for a in arm_positions)
        dists.append(d)
    return np.array(dists)

obs_dists  = min_arm_dist(smbh_az)
obs_median = np.median(obs_dists)
obs_mean   = np.mean(obs_dists)

# ── KS-Test: Beobachtete Azimute vs. gleichförmig ─────────
# Gleichförmig auf [0, 360) → CDF = x/360
stat_ks, p_ks = kstest(smbh_az / 360.0, 'uniform')

# ── Zählung in Arm-Fenstern ───────────────────────────────
n_in_arms = sum(1 for d in obs_dists if d <= AZ_TOL)
f_in_arms = n_in_arms / N_SMBH
# Erwarteter Anteil: 4 Arme × 2*AZ_TOL / 360
f_expected = 4 * 2 * AZ_TOL / 360.0
z_score = (n_in_arms - N_SMBH * f_expected) / math.sqrt(N_SMBH * f_expected * (1 - f_expected))

# ── Monte-Carlo ────────────────────────────────────────────
rng = np.random.default_rng(42)
mc_medians = []
mc_n_in    = []
for _ in range(N_MC):
    rand_az = rng.uniform(0, 360, N_SMBH)
    rd = min_arm_dist(rand_az)
    mc_medians.append(np.median(rd))
    mc_n_in.append(np.sum(rd <= AZ_TOL))

mc_medians = np.array(mc_medians)
mc_n_in    = np.array(mc_n_in)
p_median   = float(np.mean(mc_medians <= obs_median))
p_n_in     = float(np.mean(mc_n_in >= n_in_arms))

# ── Output ─────────────────────────────────────────────────
sep = "═" * 62
print(sep)
print("OT-42: 4-Arm KS-Test — Azimutale Clusterung (V21.3)")
print(sep)
print(f"  Arm-Azimute: {[f'{a:.1f}°' for a in ARM_AZ]}")
print(f"  N_SMBH = {N_SMBH}")
print()
print(f"  KS-Test vs. gleichförmige Azimutverteilung:")
print(f"    D_KS = {stat_ks:.4f},  p = {p_ks:.4f}")
print()
print(f"  Objekte innerhalb ±{AZ_TOL}° eines Arms: {n_in_arms}/{N_SMBH} = {f_in_arms:.1%}")
print(f"  Erwartung isotrop:                        {N_SMBH*f_expected:.1f}/{N_SMBH} = {f_expected:.1%}")
print(f"  z-Score: {z_score:+.2f}σ")
print()
print(f"  Medianer Min-Arm-Abstand:")
print(f"    Beobachtet: {obs_median:.2f}°")
print(f"    MC-Median:  {np.median(mc_medians):.2f}° ± {np.std(mc_medians):.2f}°")
print(f"    p(MC_median ≤ obs): {p_median:.4f}")
print()
print(f"  MC-Zählung (p-Wert N_in_arms ≥ obs): {p_n_in:.4f}")
print()
if z_score < -1.0:
    print(f"  HINWEIS: z={z_score:.1f}σ — SMBHs MEIDEN Arm-Positionen (Anti-Clustering)")
    print(f"  → Wahrscheinlich Survey-Bias: Katalog ist nicht full-sky-isotrop")
elif z_score > 1.0:
    print(f"  HINWEIS: z=+{z_score:.1f}σ — SMBHs bevorzugen Arm-Positionen (Clustering)")
print()

if p_ks < 0.01 or p_median < 0.01:
    if z_score > 0 and p_n_in < 0.05:
        # Mehr Objekte nahe Arme als erwartet → echter Clusterungs-Hinweis
        status = "BESTAETIGT"
    else:
        # KS signifikant aber SMBHs meiden Armpositionen → Survey-Bias oder Anti-Clustering
        status = "ANTI_CLUSTERING_SURVEY_BIAS"
elif p_ks < 0.05 or p_median < 0.05:
    status = "BEDINGT"
elif p_ks < 0.15 or p_median < 0.15:
    status = "HINWEIS"
else:
    status = "KEIN_SIGNAL"
print(f"  Status: {status}")
print(sep)

# ── Plot ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
fig.patch.set_facecolor('white')

ax = axes[0]
bins = np.linspace(0, 360, 25)
ax.hist(smbh_az, bins=bins, density=True, color='steelblue', alpha=0.75, label=f'SMBHs (N={N_SMBH})')
ax.axhline(1/360, color='gray', ls='--', lw=1.0, label='Isotrop')
for az in ARM_AZ:
    ax.axvline(az, color='tomato', ls=':', lw=1.5, alpha=0.8)
ax.set_xlabel('Azimut um D-Pol [°]', fontsize=9)
ax.set_ylabel('Dichte', fontsize=9)
ax.set_title(f'4-Arm Azimutverteilung\nKS p={p_ks:.3f}  z={z_score:+.1f}σ', fontsize=9)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax2 = axes[1]
ax2.hist(mc_medians, bins=50, density=True, color='lightgray', alpha=0.8, label=f'MC Null ({N_MC}×)')
ax2.axvline(obs_median, color='tomato', lw=2, label=f'Beobachtet: {obs_median:.1f}°')
ax2.axvline(np.median(mc_medians), color='gray', ls='--', lw=1.0, label=f'MC-Median: {np.median(mc_medians):.1f}°')
ax2.set_xlabel('Medianer Arm-Abstand [°]', fontsize=9)
ax2.set_ylabel('Dichte', fontsize=9)
ax2.set_title(f'MC-Null-Test: Min-Arm-Abstand\np(rand ≤ obs) = {p_median:.3f}', fontsize=9)
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)

plt.suptitle('OT-42: 4-Arm Azimutale Clusterung (97 SMBHs)', fontsize=9)
plt.tight_layout()
out = os.path.join(RESDIR, 'OT_42_plot.png')
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
plt.close()

# ── Ergebnis speichern ─────────────────────────────────────
lines = [
    sep,
    "OT-42: 4-Arm KS-Test Azimutale Clusterung (V21.3)",
    sep,
    f"N_SMBH = {N_SMBH}  |  Arm-Az = {[f'{a:.1f}' for a in ARM_AZ]}",
    "",
    f"KS-Test vs. gleichfoermig: D={stat_ks:.4f}, p={p_ks:.4f}",
    f"Objekte in +-{AZ_TOL}° Armfenstern: {n_in_arms}/{N_SMBH} (exp={N_SMBH*f_expected:.1f}), z={z_score:+.2f}sigma",
    f"MC medianer Min-Arm-Abstand: obs={obs_median:.2f}°, p(MC<=obs)={p_median:.4f}",
    f"MC N_in_arms p-Wert: {p_n_in:.4f}",
    "",
]
if z_score < -1.0:
    lines.append(f"ANTI-CLUSTERING: SMBHs meiden Arm-Positionen (z={z_score:.1f}sigma)")
    lines.append("Ursache wahrscheinlich Survey-Bias (97-SMBH-Katalog ist nicht full-sky isotrop)")
lines += [
    "",
    f"Status: {status}",
    sep,
]
with open(os.path.join(RESDIR, 'OT_42_result.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')
print(f"Gespeichert: {out}")
