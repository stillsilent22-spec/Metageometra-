"""Einfache 2D-Himmelskarte: HTM-Schalen + SMBHs"""
import math, csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SMBH_CSV = os.path.join(BASE, "results", "catalogs", "smbh_extended.csv")
OUT_PNG  = os.path.join(BASE, "results", "viz_smbh.png")

THETA0 = 58.65
POLES  = [
    (305.0,  25.0, 'D'),
    (125.0, -25.0, 'A'),
    ( 35.0,  25.0, 'S1'),
    (215.0, -25.0, 'S2'),
    (215.0,  25.0, 'S3'),
    ( 35.0, -25.0, 'S4'),
]
POLE_COLORS = {'D':'red','A':'tomato','S1':'dodgerblue','S2':'steelblue','S3':'darkorange','S4':'gold'}

# ── SMBHs laden ─────────────────────────────────────────────────────────────────
def radec_to_lb(ra, dec):
    ra = math.radians(ra); dc = math.radians(dec)
    RA_NGP = math.radians(192.85948); DEC_NGP = math.radians(27.12825)
    b = math.asin(math.sin(dc)*math.sin(DEC_NGP) +
                  math.cos(dc)*math.cos(DEC_NGP)*math.cos(ra - RA_NGP))
    y = math.cos(dc)*math.sin(ra - RA_NGP)
    x = math.cos(dc)*math.sin(DEC_NGP)*math.cos(ra - RA_NGP) - math.sin(dc)*math.cos(DEC_NGP)
    l = (math.degrees(math.atan2(y, x)) + 122.93192) % 360.0
    return l, math.degrees(b)

def lb_xyz(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

smbhs = []
with open(SMBH_CSV, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try:
            ra = float(r['RA_deg']); dec = float(r['Dec_deg'])
            lM = float(r['logMbh'])
            l, b = radec_to_lb(ra, dec)
            # crossing density
            nc = 0
            for pl, pb, _ in POLES:
                for n in [1, 2, 3]:
                    pv = lb_xyz(pl, pb); sv = lb_xyz(l, b)
                    dot = max(-1.0, min(1.0, float(np.dot(pv, sv))))
                    if abs(math.degrees(math.acos(dot)) - n*THETA0) < 8.0:
                        nc += 1
            dist = float(r['dist_mpc']) if r.get('dist_mpc', '').strip() not in ('', 'nan') else None
            smbhs.append({'l': l, 'b': b, 'logM': lM, 'nc': nc, 'name': r['Name'].strip(), 'dist': dist})
        except Exception:
            continue

# ── Schalenringe zeichnen ────────────────────────────────────────────────────────
def shell_ring(pole_l, pole_b, theta_deg, n_pts=500):
    """Punkte auf dem Ring in abstand theta_deg vom Pol"""
    pv = lb_xyz(pole_l, pole_b)
    # Zwei senkrechte Vektoren zum Pol konstruieren
    ref = np.array([0,0,1]) if abs(pv[2]) < 0.9 else np.array([1,0,0])
    e1 = np.cross(pv, ref); e1 /= np.linalg.norm(e1)
    e2 = np.cross(pv, e1); e2 /= np.linalg.norm(e2)
    ct = math.cos(math.radians(theta_deg))
    st = math.sin(math.radians(theta_deg))
    ls, bs = [], []
    for phi in np.linspace(0, 2*math.pi, n_pts):
        v = ct*pv + st*(math.cos(phi)*e1 + math.sin(phi)*e2)
        r = math.sqrt(float(v@v))
        b_ = math.degrees(math.asin(max(-1, min(1, v[2]/r))))
        l_ = math.degrees(math.atan2(v[1], v[0])) % 360.0
        ls.append(l_); bs.append(b_)
    # Sprünge unterdrücken (Linienbruch bei Wrap-around)
    segments_l, segments_b = [[]], [[]]
    for i in range(len(ls)):
        if i > 0 and abs(ls[i] - ls[i-1]) > 90:
            segments_l.append([]); segments_b.append([])
        segments_l[-1].append(ls[i]); segments_b[-1].append(bs[i])
    return segments_l, segments_b

# ── Plot ─────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 7), facecolor='#0a0a1a')
ax.set_facecolor('#0a0a1a')

# Schalen-Ringe für alle Pole, n=1,2,3
shell_alpha = [0.35, 0.22, 0.12]
shell_lw    = [0.7,  0.5,  0.35]
for pl, pb, pn in POLES:
    pc = POLE_COLORS[pn]
    for ni, n in enumerate([1, 2, 3]):
        segs_l, segs_b = shell_ring(pl, pb, n * THETA0)
        for sl, sb in zip(segs_l, segs_b):
            ax.plot(sl, sb, '-', color=pc, alpha=shell_alpha[ni], lw=shell_lw[ni])

# Pole als Kreuz
for pl, pb, pn in POLES:
    pc = POLE_COLORS[pn]
    ax.plot(pl, pb, 'x', color=pc, ms=7, mew=1.5, alpha=0.8)
    ax.text(pl+1.5, pb+1.5, pn, color=pc, fontsize=7, alpha=0.9)

# Galaktisches Zentrum markieren
ax.axvline(0,   color='white', lw=0.3, alpha=0.2, ls='--')
ax.axhline(0,   color='white', lw=0.3, alpha=0.2, ls='--')
ax.axvline(180, color='white', lw=0.3, alpha=0.15, ls=':')

# SMBHs
nc_max = max(s['nc'] for s in smbhs)
cmap = plt.cm.plasma
for s in smbhs:
    sz  = max(10, (s['logM'] - 6.5) * 18)
    col = cmap(s['nc'] / max(nc_max, 1))
    ax.scatter(s['l'], s['b'], s=sz, c=[col], zorder=5, edgecolors='none', alpha=0.85)

# Labels nur für nc >= 3
labeled = sorted([s for s in smbhs if s['nc'] >= 3], key=lambda x: -x['nc'])
for s in labeled[:20]:
    ax.text(s['l']+1.5, s['b']+1.5, s['name'], color='white', fontsize=6.5,
            alpha=0.9, va='bottom')

# Legende
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, nc_max))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.018, pad=0.01)
cbar.set_label('nc (Kreuzungsdichte)', color='white', fontsize=9)
cbar.ax.yaxis.set_tick_params(color='white'); plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

sizes_leg = [mpatches.Patch(color='gray', label='Größe ∝ log(M_BH)'),
             mpatches.Patch(color='gray', label='Farbe = nc')]
ax.legend(handles=sizes_leg, loc='lower left', fontsize=8, facecolor='#111', edgecolor='gray',
          labelcolor='white', framealpha=0.6)

ax.set_xlim(0, 360); ax.set_ylim(-90, 90)
ax.set_xlabel('Galaktische Länge l [°]', color='white', fontsize=10)
ax.set_ylabel('Galaktische Breite b [°]', color='white', fontsize=10)
ax.set_title('SMBHs & HTM-Schalen (galaktische Koordinaten)', color='white', fontsize=12, pad=10)
ax.tick_params(colors='white'); ax.spines[:].set_color('#555')
ax.set_xticks(range(0, 361, 30)); ax.set_yticks(range(-90, 91, 30))

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
print(f"Gespeichert: {OUT_PNG}")

# ═══════════════════════════════════════════════════════════════════════════════
# ZWEITES BILD: 2D-Querschnitt durch die Dualitätssphäre
# ═══════════════════════════════════════════════════════════════════════════════
R_S = 4228.3
OUT_CROSS = os.path.join(BASE, "results", "viz_crosssection.png")

fig2, ax2 = plt.subplots(figsize=(10, 10), facecolor='#0a0a1a')
ax2.set_facecolor('#0a0a1a')
ax2.set_aspect('equal')

# Dualitätssphäre
theta_arr = np.linspace(0, 2*math.pi, 800)
ax2.plot(R_S*np.cos(theta_arr), R_S*np.sin(theta_arr),
         color='white', lw=2.0, alpha=0.9, label=f'Dualitätssphäre R_S={R_S:.0f} Mpc', zorder=3)

# Thales-Schalen R_n = R_S * n * θ₀/360
thales_colors = ['#ff6b35', '#ffd700', '#00bfff', '#7fff00', '#ff69b4', '#da70d6']
for n in range(1, 7):
    Rn = R_S * n * THETA0 / 360.0
    col = thales_colors[n % len(thales_colors)]
    ax2.plot(Rn*np.cos(theta_arr), Rn*np.sin(theta_arr),
             color=col, lw=1.0, alpha=0.5, ls='--', label=f'Thales n={n}: {Rn:.0f} Mpc')

# Beobachter
ax2.plot(0, 0, 'o', color='white', ms=6, zorder=10)
ax2.text(80, 80, 'Beobachter\n(Milchstraße)', color='white', fontsize=8.5, va='bottom')

# D-Pol Richtung einzeichnen
dpol_l, dpol_b = 305.0, 25.0
dpol_xyz = lb_xyz(dpol_l, dpol_b)
arrow_len = R_S * 1.08
ax2.annotate('', xy=(arrow_len*dpol_xyz[0], arrow_len*dpol_xyz[1]),
             xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
ax2.text(arrow_len*dpol_xyz[0]*1.02, arrow_len*dpol_xyz[1]*1.02,
         'D-Pol\n(l=305°, b=25°)', color='red', fontsize=8, ha='center')

# SMBHs mit bekannter Distanz — projiziert auf galaktische l-Ebene (b ignoriert)
smbhs_3d = [s for s in smbhs if s.get('dist') is not None]
if smbhs_3d:
    xs = [s['dist'] * math.cos(math.radians(s['l'])) for s in smbhs_3d]
    ys = [s['dist'] * math.sin(math.radians(s['l'])) for s in smbhs_3d]
    szs = [max(20, (s['logM'] - 6.5) * 25) for s in smbhs_3d]
    cols = [s['nc'] / max(nc_max, 1) for s in smbhs_3d]
    ax2.scatter(xs, ys, s=szs, c=cols, cmap='plasma', vmin=0, vmax=1,
                zorder=5, edgecolors='white', linewidths=0.3, alpha=0.9)
    for s, x, y in zip(smbhs_3d, xs, ys):
        if s['nc'] >= 3:
            ax2.text(x+50, y+50, s['name'], color='white', fontsize=6.5, alpha=0.9)

# Gitternetz
for r in [500, 1000, 2000, 3000]:
    ax2.plot(r*np.cos(theta_arr), r*np.sin(theta_arr),
             color='white', lw=0.3, alpha=0.12)
    ax2.text(r*0.707+30, r*0.707+30, f'{r} Mpc', color='white', fontsize=6.5, alpha=0.35)

max_r = R_S * 1.15
ax2.set_xlim(-max_r, max_r); ax2.set_ylim(-max_r, max_r)
ax2.set_xlabel('x [Mpc]  (galaktische Ebene)', color='white', fontsize=10)
ax2.set_ylabel('y [Mpc]  (galaktische Ebene)', color='white', fontsize=10)
ax2.set_title('Querschnitt: Dualitätssphäre + Thales-Schalen + SMBHs', color='white', fontsize=12, pad=10)
ax2.tick_params(colors='white'); ax2.spines[:].set_color('#555')

leg = ax2.legend(loc='lower right', fontsize=7.5, facecolor='#111', edgecolor='gray',
                 labelcolor='white', framealpha=0.7)
plt.tight_layout()
plt.savefig(OUT_CROSS, dpi=150, bbox_inches='tight', facecolor=fig2.get_facecolor())
print(f"Gespeichert: {OUT_CROSS}")
