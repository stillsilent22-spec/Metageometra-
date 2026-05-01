"""
OT-13: SRM Eigenmode-Spektrum — Schalenüberschuss-Methode (V21.3 Rev.3)

Observable: δ_n = (N_smbh_obs - N_smbh_exp) / sqrt(N_smbh_exp)
  N_smbh_exp = N_SMBH × (N_mq_in_band / N_mq_total)
  = SMBHs im Band RELATIV zu Milliquas-Hintergrunddichte
  = bereinigt um nicht-isotrope Survey-Abdeckung
Erwartung: δ_n ∝ n^{-D_eff} = n^{-1/3} (SRM Spektrum)
"""
import math, os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE   = os.path.dirname(os.path.abspath(__file__))
RESDIR = os.path.join(BASE, "results")
CATDIR = os.path.join(RESDIR, "catalogs")
os.makedirs(RESDIR, exist_ok=True)

C     = 2.998e8; T0 = 4.352e17
A0    = C / (2*math.pi*T0)
R_S   = 42.3 * 3.0857e19   # 42.3 kpc in meters
RHO_0 = 0.406
D_GEO  = math.log(2)/math.log(1/RHO_0)
D_EFF  = D_GEO / (3*D_GEO)  # = 1/3

SHELL_BW = 2.0   # Grad Halbbreite pro Schale
NSHELLS  = 6

def shell_theta(n, chi=59.1, delta=1.0):
    d = math.radians(delta); c = math.radians(chi)
    return math.degrees(math.acos(max(-1., min(1., math.cos(n*d)*math.cos(n*c)))))

shells = np.array([shell_theta(n) for n in range(1, NSHELLS+1)])

# D-Pol Richtungsvektor: l=305°, b=25° (galaktisch)
# Konvertierung galaktisch → äquatorial mittels IAU-Standardrotation
# NGP: RA=192.85948°, Dec=27.12825°, l_NCP=122.93192°
def dpole_eq_unit():
    ra_ngp = math.radians(192.85948)
    dec_ngp = math.radians(27.12825)
    l_ncp  = math.radians(122.93192)
    l_d = math.radians(305.0); b_d = math.radians(25.0)
    # Galactic to ecliptic unit vector
    sin_b = math.sin(b_d); cos_b = math.cos(b_d)
    sin_l = math.sin(l_d - l_ncp); cos_l = math.cos(l_d - l_ncp)
    sin_dec_ngp = math.sin(dec_ngp); cos_dec_ngp = math.cos(dec_ngp)

    # Rotation: galactic (l,b) → equatorial (RA, Dec)
    # z-component in galactic frame:
    z_g  = sin_b
    y_g  = cos_b * sin_l   # sin(l - l_ncp)
    x_g  = cos_b * cos_l   # cos(l - l_ncp)

    # z_eq = sin(Dec)
    z_eq = sin_dec_ngp * z_g + cos_dec_ngp * y_g  # not quite — use matrix
    # Standard matrix form
    # Equation: eq = R @ gal
    # R[0] = (-sin(l_ncp)*sin(ra_ngp) - cos(l_ncp)*cos(ra_ngp)*sin(dec_ngp),
    #          cos(l_ncp)*sin(ra_ngp) - sin(l_ncp)*cos(ra_ngp)*sin(dec_ngp),
    #          cos(ra_ngp)*cos(dec_ngp))  ... not needed here

    # Use direct formula: angle between D-pole and each object
    # Convert D-pole to equatorial unit vector using formula:
    # sin(dec_dp) = sin(dec_ngp)*sin(b) + cos(dec_ngp)*cos(b)*sin(l - l_ncp)
    dec_dp = math.asin(sin_dec_ngp * sin_b + cos_dec_ngp * cos_b * sin_l)
    ra_dp  = math.atan2(
        cos_b * cos_l,
        sin_b * cos_dec_ngp - cos_b * sin_dec_ngp * sin_l
    ) + ra_ngp
    # unit vector
    xd = math.cos(dec_dp) * math.cos(ra_dp)
    yd = math.cos(dec_dp) * math.sin(ra_dp)
    zd = math.sin(dec_dp)
    return xd, yd, zd

DPOLE_UV = dpole_eq_unit()

def theta_from_radec(ra_deg, dec_deg):
    """Winkel in Grad zwischen (RA, Dec) und D-Pol."""
    ra = math.radians(ra_deg); dec = math.radians(dec_deg)
    xp = math.cos(dec)*math.cos(ra)
    yp = math.cos(dec)*math.sin(ra)
    zp = math.sin(dec)
    dot = xp*DPOLE_UV[0] + yp*DPOLE_UV[1] + zp*DPOLE_UV[2]
    return math.degrees(math.acos(max(-1., min(1., dot))))

# ── Lade Milliquas Quasare (primäres AGN-Sample) ──────────────────────────────
mq_theta = []
mq_path  = os.path.join(RESDIR, "milliquas_sample.csv")
if os.path.exists(mq_path):
    with open(mq_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                ra  = float(row['RAJ2000'])
                dec = float(row['DEJ2000'])
                mq_theta.append(theta_from_radec(ra, dec))
            except (KeyError, ValueError):
                pass
mq_theta = np.array(mq_theta)

# ── Lade direkte SMBH-Messungen ───────────────────────────────────────────────
smbh_theta = []
smbh_path = os.path.join(CATDIR, "smbh_extended.csv")
if not os.path.exists(smbh_path):
    smbh_path = os.path.join(CATDIR, "smbh_catalog_combined.csv")
if os.path.exists(smbh_path):
    with open(smbh_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                smbh_theta.append(float(row['theta_dpole']))
            except (KeyError, ValueError):
                pass
smbh_theta = np.array(smbh_theta)

# Kombiniertes Sample: Milliquas als Hintergrund, SMBH als Signal
# all_theta nur für Plot; primäres Test-Sample = smbh_theta vs mq als BG
N_MQ   = len(mq_theta)
N_DIR  = len(smbh_theta)
N_SMBH = N_DIR  # Testsample: nur direkte SMBH-Messungen
all_theta = np.concatenate([mq_theta, smbh_theta]) if len(mq_theta) > 0 else smbh_theta

# ── Schalenüberschuss δ_n berechnen ─────────────────────────────────────────
# Methode Rev.3: Milliquas als empirischer Survey-Hintergrund
# f_mq(n) = Anteil Milliquas in ±bw um θ_n  → bereinigt um nicht-isotrope Abdeckung
# N_smbh_exp(n) = N_SMBH × f_mq(n)

delta_n_arr = []
N_obs_arr  = []
N_exp_arr  = []
n_arr = np.arange(1, NSHELLS + 1)

for i, (n, theta_deg) in enumerate(zip(n_arr, shells)):
    bw = SHELL_BW

    # SMBH-Zählung im Band (nur direkte Messungen als Signal)
    N_obs = float(np.sum((smbh_theta >= theta_deg - bw) & (smbh_theta <= theta_deg + bw)))

    # Hintergrund-Erwartung: Milliquas-Dichte in gleichem Band
    if N_MQ > 0:
        N_mq_band = float(np.sum((mq_theta >= theta_deg - bw) & (mq_theta <= theta_deg + bw)))
        f_mq = N_mq_band / N_MQ
        N_exp = N_SMBH * f_mq
    else:
        # Fallback isotrope Geometrie
        lo = max(0., theta_deg - bw); hi = min(180., theta_deg + bw)
        f_geom = (math.cos(math.radians(lo)) - math.cos(math.radians(hi))) / 2.0
        N_exp = N_SMBH * f_geom

    sigma = max(0.5, math.sqrt(max(0.01, N_exp)))
    delta = (N_obs - N_exp) / sigma

    N_obs_arr.append(N_obs)
    N_exp_arr.append(N_exp)
    delta_n_arr.append(delta)

delta_n_arr = np.array(delta_n_arr)
N_obs_arr   = np.array(N_obs_arr)
N_exp_arr   = np.array(N_exp_arr)

# ── Log-Log Fit δ_n vs n ─────────────────────────────────────────────────────
# Nur Schalen mit δ_n > 0 für log-log Fit (positive Überschüsse)
mask_pos = delta_n_arr > 0
n_fit  = n_arr[mask_pos].astype(float)
dn_fit = delta_n_arr[mask_pos]

if len(n_fit) >= 2:
    log_n = np.log(n_fit)
    log_d = np.log(dn_fit)
    slope, intercept = np.polyfit(log_n, log_d, 1)
else:
    slope, intercept = float('nan'), float('nan')

expected_slope = -D_EFF  # -0.3333
deviation = abs(slope - expected_slope) if not math.isnan(slope) else float('nan')

# Ebenfalls alte sin²-Methode für Vergleich (Methode 1)
EXPONENT_OLD = -2.0 / (2.0 - D_EFF)
N_MODES = 1000
r_grid = np.linspace(1e16, R_S, N_MODES)
def rho_SRM(r):
    gamma = 1 - D_EFF/2
    return (1 + gamma*(r/R_S))**EXPONENT_OLD
rho_arr = np.array([rho_SRM(r) for r in r_grid])
def compute_eigenmode(n_mode):
    psi_sq = np.sin(n_mode * math.pi * r_grid / R_S)**2
    return np.trapz(rho_arr * psi_sq, r_grid)
N_MAX = 10
w_raw = np.array([compute_eigenmode(n) for n in range(1, N_MAX+1)])
w_norm = w_raw / w_raw.sum()
slope_old, _ = np.polyfit(np.log(np.arange(1, N_MAX+1)), np.log(w_norm), 1)

print("══════════════════════════════════════════════════════")
print("OT-13 Rev.3: SRM Eigenmode-Spektrum — Survey-bereinigt (V21.3)")
print("  Methode: SMBH-Überschuss relativ zu Milliquas-Hintergrund")
print("══════════════════════════════════════════════════════")
print(f"  r_s         = {R_S:.4e} m = {R_S/3.086e19:.1f} kpc")
print(f"  D_eff       = {D_EFF:.4f}")
print(f"  SMBHs (Signal) = {N_DIR}")
print(f"  Milliquas (BG) = {N_MQ}  (Survey-Hintergrund)")
print()
print(f"  n  | θ_n [°]  | N_obs | N_exp  | δ_n")
print(f"  {'─'*2}+{'─'*9}+{'─'*7}+{'─'*8}+{'─'*8}")
for i, (n, th, no, ne, dn) in enumerate(zip(n_arr, shells, N_obs_arr, N_exp_arr, delta_n_arr)):
    print(f"  {n:>2} | {th:>7.2f}  | {int(no):>5} | {ne:>6.2f} | {dn:>+7.2f}")

print()
print(f"  Log-Log Fit (δ_n > 0): Steigung = {slope:.4f}")
print(f"  Erwartung (n^-D_eff) : Steigung = {expected_slope:.4f}")
print(f"  Abweichung            : {deviation:.4f}")
print()
print(f"  [Methode 1 sin²-Moden: Steigung = {slope_old:.4f} — falsche Observable, Ref.]")
print()

if not math.isnan(deviation) and deviation < 0.15:
    status = "BESTAETIGT"
elif not math.isnan(deviation) and deviation < 0.35:
    status = "BEDINGT"
else:
    status = "INCONCLUSIVE"
print(f"  Status: {status}")
print("══════════════════════════════════════════════════════")

# ── Plot ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
fig.patch.set_facecolor('white')

ax = axes[0]
colors = ['steelblue' if d > 0 else 'tomato' for d in delta_n_arr]
ax.bar(n_arr, delta_n_arr, color=colors, alpha=0.85)
# Show expected n^{-1/3} line scaled to n=1
if delta_n_arr[0] > 0:
    n_line = np.linspace(1, NSHELLS, 100)
    expected_line = delta_n_arr[0] * n_line**expected_slope
    ax.plot(n_line, expected_line, 'g--', lw=1.5, label=f'Erw.: n^{expected_slope:.3f}')
if not math.isnan(slope):
    n_line = np.linspace(1, NSHELLS, 100)
    fit_line = np.exp(intercept) * n_line**slope if np.any(mask_pos) else expected_line
    ax.plot(n_line, fit_line, 'r:', lw=1.5, label=f'Fit: n^{slope:.3f}')
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('Schale n', fontsize=9); ax.set_ylabel('δ_n = (N_obs - N_exp)/σ', fontsize=9)
ax.set_title(f'Schalenüberschuss δ_n\nSlope={slope:.3f} vs Erw.={expected_slope:.3f}', fontsize=9)
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_xticks(n_arr)

ax2 = axes[1]
# Log-log only for positive δ_n
if np.any(mask_pos):
    ax2.loglog(n_fit, dn_fit, 'bo', ms=8, label='δ_n > 0')
    ax2.loglog(n_line, fit_line, 'r--', lw=1.5, label=f'Fit: n^{slope:.3f}')
    if delta_n_arr[0] > 0:
        ax2.loglog(n_line, expected_line, 'g:', lw=1.5, label=f'n^{expected_slope:.3f}')
ax2.set_xlabel('n (log)', fontsize=9); ax2.set_ylabel('δ_n (log)', fontsize=9)
ax2.set_title('Log-Log: Spektrale Skalierung', fontsize=9)
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)

plt.suptitle('OT-13 Rev.3: SRM Survey-bereinigt (SMBH vs Milliquas BG)', fontsize=9)
plt.tight_layout()
out = os.path.join(RESDIR, 'eigenmode_plot.png')
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
plt.close()

lines = [
    "════════════════════════════════════════════════════════════",
    "OT-13 Rev.3: SRM Eigenmode-Spektrum — Survey-bereinigt (V21.3)",
    "════════════════════════════════════════════════════════════",
    f"r_s = {R_S:.4e} m = {R_S/3.086e19:.1f} kpc",
    f"D_eff = {D_EFF:.4f}",
    f"SMBHs (Signal) = {N_DIR}  |  Milliquas (BG) = {N_MQ}",
    f"Methode: SMBH-Überschuss vs Survey-Hintergrund (nicht isotrop)",
    "",
    "Schalenüberschuss δ_n = (N_obs - N_exp) / sqrt(N_exp):",
]
for i, (n, th, no, ne, dn) in enumerate(zip(n_arr, shells, N_obs_arr, N_exp_arr, delta_n_arr)):
    lines.append(f"  n={n}  θ={th:.2f}°  N_obs={int(no)}  N_exp={ne:.2f}  δ={dn:+.2f}")
lines += [
    "",
    f"Log-Log Fit Steigung: {slope:.4f}",
    f"Erwartung (n^-D_eff): {expected_slope:.4f}",
    f"Abweichung:           {deviation:.4f}",
    "",
    f"[Methode 1 sin²-Moden Steigung: {slope_old:.4f} — falsche Observable, nur Ref.]",
    "",
    f"Status: {status}",
    "════════════════════════════════════════════════════════════",
]
with open(os.path.join(RESDIR,'OT_13_result.txt'),'w',encoding='utf-8') as f:
    f.write('\n'.join(lines)+'\n')
print(f"Gespeichert: {out}")
