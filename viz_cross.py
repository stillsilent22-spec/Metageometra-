"""2D-Querschnitt: Dualitätssphäre + SMBHs, log-radialer Maßstab"""
import math, csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE     = os.path.dirname(os.path.abspath(__file__))
SMBH_CSV = os.path.join(BASE, "results", "catalogs", "smbh_extended.csv")
OUT      = os.path.join(BASE, "results", "viz_crosssection.png")

R_S = 4228.3

def radec_to_lb(ra, dec):
    ra = math.radians(ra); dc = math.radians(dec)
    RN = math.radians(192.85948); DN = math.radians(27.12825)
    b = math.asin(math.sin(dc)*math.sin(DN) + math.cos(dc)*math.cos(DN)*math.cos(ra-RN))
    y = math.cos(dc)*math.sin(ra-RN)
    x = math.cos(dc)*math.sin(DN)*math.cos(ra-RN) - math.sin(dc)*math.cos(DN)
    return (math.degrees(math.atan2(y,x)) + 122.93192) % 360.0, math.degrees(b)

smbhs = []
with open(SMBH_CSV, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try:
            l, b = radec_to_lb(float(r['RA_deg']), float(r['Dec_deg']))
            dist = float(r['dist_mpc']) if r.get('dist_mpc','').strip() not in ('','nan') else None
            if dist and dist > 0:
                smbhs.append({'l': l, 'b': b, 'dist': dist})
        except Exception:
            pass

# nur galaktische Scheibe: |b| <= 25°
smbhs = [s for s in smbhs if abs(s['b']) <= 25.0]

fig, ax = plt.subplots(figsize=(9, 9), facecolor='#0a0a1a')
ax.set_facecolor('#0a0a1a')
ax.set_aspect('equal')
ax.axis('off')

t = np.linspace(0, 2*math.pi, 1000)
log_RS = math.log10(R_S)

# Referenzkreise (log-Skala)
for r_ref, lbl in [(1, '1 Mpc'), (10, '10 Mpc'), (100, '100 Mpc'), (1000, '1000 Mpc')]:
    lr = math.log10(r_ref)
    ax.plot(lr*np.cos(t), lr*np.sin(t), color='white', lw=0.4, alpha=0.18, zorder=1)
    ax.text(0, lr + 0.04, lbl, color='white', fontsize=7.5, alpha=0.35,
            ha='center', va='bottom')

# Dualitätssphäre
ax.plot(log_RS*np.cos(t), log_RS*np.sin(t), color='white', lw=1.8, alpha=0.9, zorder=3)
ax.text(0, log_RS + 0.04, f'Dualitätssphäre  {R_S:.0f} Mpc',
        color='white', fontsize=9, ha='center', va='bottom', alpha=0.85)

# SMBHs als Kreuze
for s in smbhs:
    lr = math.log10(s['dist'])
    x = lr * math.cos(math.radians(s['l']))
    y = lr * math.sin(math.radians(s['l']))
    ax.plot(x, y, '+', color='cyan', ms=11, mew=1.4, zorder=5)

# Beobachter
ax.plot(0, 0, 'o', color='white', ms=4, zorder=10)

# Achsenbeschriftung
ax.text(0, -log_RS * 1.12, f'Horizontalschnitt  |b| ≤ 25°  ·  {len(smbhs)} SMBHs  ·  log-Skala',
        color='white', fontsize=8, ha='center', alpha=0.45)

lim = log_RS * 1.18
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

plt.tight_layout(pad=0.3)
plt.savefig(OUT, dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
print(f"Gespeichert: {OUT}")

