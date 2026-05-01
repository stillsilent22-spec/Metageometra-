"""
Minimalist visualization: Tesseract D-Pol-Projektion + V20-SMBH-Treffer
- Kleiner Punkt: D-Pol (Tesserakt-Projektionsachse)
- Duenne Kreise: HTM-Schalen theta_n
- Kreuze schwarz: 28 SMBH-Treffer (±2 deg)
- Kreuze hell: restliche 69 SMBHs
Koordinaten: flache Polare (theta * cos/sin PA) zentriert auf D-Pol
"""
import math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

BASE   = os.path.dirname(os.path.abspath(__file__))
CATS   = os.path.join(BASE, "results", "catalogs")
RESDIR = os.path.join(BASE, "results")

# ── V20-Parameter ──────────────────────────────────────────
chi_deg   = 59.1;  delta_deg = 1.0
chi_rad   = math.radians(chi_deg)
delta_rad = math.radians(delta_deg)

def shell_angle(n):
    a = math.cos(n*delta_rad)*math.cos(n*chi_rad)
    return math.degrees(math.acos(max(-1., min(1., a))))

shells = np.array([shell_angle(n) for n in range(1, 7)])

# ── D-Pol: gal(l=305, b=+25) -> äquat ──────────────────────
_RA_NGP  = math.radians(192.859508)
_Dec_NGP = math.radians(27.128336)
_l_NCP   = math.radians(122.932)
_l_dp    = math.radians(305.); _b_dp = math.radians(25.)
_sin_dec = (math.sin(_b_dp)*math.sin(_Dec_NGP)
            + math.cos(_b_dp)*math.cos(_Dec_NGP)*math.cos(_l_NCP - _l_dp))
_DP_Dec  = math.asin(max(-1., min(1., _sin_dec)))
_cos_dec = math.cos(_DP_Dec)
_sin_dra = -math.cos(_b_dp)*math.sin(_l_dp - _l_NCP) / _cos_dec
_cos_dra = (math.sin(_b_dp) - math.sin(_Dec_NGP)*_sin_dec) / (math.cos(_Dec_NGP)*_cos_dec)
_DP_RA   = math.atan2(_sin_dra, _cos_dra) + _RA_NGP

def dist_pa(ra_r, dec_r):
    """Angular distance and position angle from D-pole"""
    cos_a = (math.sin(_DP_Dec)*math.sin(dec_r)
             + math.cos(_DP_Dec)*math.cos(dec_r)*math.cos(_DP_RA - ra_r))
    th = math.acos(max(-1., min(1., cos_a)))
    dra = ra_r - _DP_RA
    y   = math.sin(dra)*math.cos(dec_r)
    x   = math.cos(_DP_Dec)*math.sin(dec_r) - math.sin(_DP_Dec)*math.cos(dec_r)*math.cos(dra)
    pa  = math.atan2(y, x)
    return math.degrees(th), pa

# ── SMBH-Katalog laden ─────────────────────────────────────
df = pd.read_csv(os.path.join(CATS, "smbh_extended.csv"), comment='#')
tol = 5.0  # ±5 deg

thetas, pas, hits = [], [], []
for _, row in df.iterrows():
    th, pa = dist_pa(math.radians(float(row['RA_deg'])),
                     math.radians(float(row['Dec_deg'])))
    delta  = float(np.min(np.abs(th - shells)))
    thetas.append(th); pas.append(pa); hits.append(delta < tol)

thetas = np.array(thetas); pas = np.array(pas); hits = np.array(hits)

# Belt ±15° um horizontale Achse
belt_mask = np.abs(np.sin(pas)) < np.sin(np.radians(15))   # sin(15°) ≈ 0.259

# 2D polar: x = theta*cos(pa), y = theta*sin(pa)
xs = thetas * np.cos(pas)
ys = thetas * np.sin(pas)

# ── Plot ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 5))
ax.set_aspect('equal')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Äusserster Referenzkreis 180°
phi = np.linspace(0, 2*np.pi, 720)
ax.plot(180*np.cos(phi), 180*np.sin(phi), lw=0.3, color='#ececec', zorder=0)

# Schalenkreise
for sh in sorted(shells):
    ax.plot(sh*np.cos(phi), sh*np.sin(phi), lw=0.5, color='#c0c0c0', zorder=1, alpha=0.7)

# Alle 97 SMBHs — kompletter Umkreis (alle Azimuthe)
# außerhalb Belt: winzige Punkte
ax.scatter(xs[~belt_mask & ~hits], ys[~belt_mask & ~hits],
           marker='.', s=5, c='#e0e0e0', zorder=2)
ax.scatter(xs[~belt_mask & hits],  ys[~belt_mask & hits],
           marker='+', s=16, c='#909090', linewidths=0.6, zorder=3)
# innerhalb Belt ±15°: gleiche Darstellung wie zuvor
ax.scatter(xs[belt_mask & ~hits], ys[belt_mask & ~hits],
           marker='+', s=14, c='#c0c0c0', linewidths=0.5, zorder=4)
ax.scatter(xs[belt_mask & hits],  ys[belt_mask & hits],
           marker='+', s=32, c='black', linewidths=0.9, zorder=5)

# Tesserakt / D-Pol — winziger Punkt im Zentrum
ax.plot(0, 0, '.', color='black', markersize=2.5, zorder=5)

ax.set_xlim(-188, 188)
ax.set_ylim(-188, 188)
ax.axis('off')

plt.tight_layout(pad=0.2)
out = os.path.join(RESDIR, 'viz_tesseract_smbh.png')
plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
print(f"OK: {out}  ({np.sum(hits)} Treffer markiert)")
