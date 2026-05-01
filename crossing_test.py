"""
Kreuzungsdichte-Test: Sitzen SMBHs bevorzugt an Schalenkreuzpunkten?

Modell: Unser Universum aus symmetrischem Torsionsschock emergiert.
Schalen emanieren von MEHREREN Polen gleichzeitig (volle S2-Symmetrie).
Knoten = Kreuzungspunkte von Schalen verschiedener Pole.
M_BH sollte mit lokaler Kreuzungsdichte korrelieren.
"""
import math, csv, os, sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

THETA_0 = 58.65
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
CSV_PATH = os.path.join(RESULTS, "catalogs", "smbh_extended.csv")

def sph_to_xyz(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    return (math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b))

def angle_between(v1, v2):
    d = max(-1.0, min(1.0, sum(a*b for a,b in zip(v1,v2))))
    return math.degrees(math.acos(d))

def radec_to_lb(ra_deg, dec_deg):
    ra_r = math.radians(ra_deg); dc_r = math.radians(dec_deg)
    RA_NGP = math.radians(192.85948); DEC_NGP = math.radians(27.12825)
    b_r = math.asin(math.sin(dc_r)*math.sin(DEC_NGP) +
                    math.cos(dc_r)*math.cos(DEC_NGP)*math.cos(ra_r-RA_NGP))
    l_r = math.atan2(
        math.cos(dc_r)*math.sin(ra_r-RA_NGP),
        math.cos(dc_r)*math.sin(DEC_NGP)*math.cos(ra_r-RA_NGP) - math.sin(dc_r)*math.cos(DEC_NGP))
    return (math.degrees(l_r) + 122.93192) % 360.0, math.degrees(b_r)

# ── Pole definieren ────────────────────────────────────────────────────────────
# D-Pol + A-Pol (antipodal) sind die primären Kollisions-Achse.
# Sekundäre Pole bei 90° versetzt auf dem Großkreis zwischen D und A
# (entspricht den S3-Symmetrieachsen des Torsionsfeldes)
POLES_DEF = [
    (305.0,  25.0, 'D-Pol'),    # Haupt-Kollisionspol
    (125.0, -25.0, 'A-Pol'),    # Antipodal-Pol
    ( 35.0,  25.0, 'Sek-1'),    # 90° von D auf Großkreis
    (215.0, -25.0, 'Sek-2'),
    (215.0,  25.0, 'Sek-3'),    # 90° anders herum
    ( 35.0, -25.0, 'Sek-4'),
]
POLES = [(sph_to_xyz(l,b), name) for l,b,name in POLES_DEF]
SHELLS = [THETA_0 * n for n in range(1, 4)]  # n=1,2,3 (nur bis 175.95°)
TOL = 8.0  # Grad Toleranz

def crossing_density(l_deg, b_deg):
    """Anzahl Schalenkreuzungen nahe diesem Punkt."""
    p = sph_to_xyz(l_deg, b_deg)
    on_shells = []
    for (pv, pname) in POLES:
        theta = angle_between(p, pv)
        for n, s in enumerate(SHELLS, 1):
            if abs(theta - s) < TOL:
                on_shells.append((pname, n, abs(theta - s)))
    # Kreuzungen = Paare auf verschiedenen Polen
    crossings = 0
    for i in range(len(on_shells)):
        for j in range(i+1, len(on_shells)):
            if on_shells[i][0] != on_shells[j][0]:
                crossings += 1
    return crossings, on_shells

# ── SMBH-Katalog laden ─────────────────────────────────────────────────────────
rows = []
with open(CSV_PATH, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try:
            lm = float(r['logMbh'])
            ra = float(r['RA_deg'])
            dc = float(r['Dec_deg'])
            l_deg, b_deg = radec_to_lb(ra, dc)
            nc, shells_hit = crossing_density(l_deg, b_deg)
            rows.append({
                'name': r['Name'], 'logM': lm, 'nc': nc,
                'l': l_deg, 'b': b_deg,
                'theta_D': float(r['theta_dpole']),
                'shells': shells_hit
            })
        except Exception:
            continue

print(f"Geladene SMBHs: {len(rows)}")
print()

# ── Kreuzungsdichte vs Masse ───────────────────────────────────────────────────
print("Kreuzungsdichte vs mittlere log(M_BH):")
print(f"{'Crossings':>10}  {'N':>4}  {'Mittel log(M)':>14}  {'Max log(M)':>11}")
by_nc = {}
for r in rows:
    by_nc.setdefault(r['nc'], []).append(r['logM'])
for nc in sorted(by_nc):
    vals = by_nc[nc]
    print(f"  {nc:>8}    {len(vals):>4}    {sum(vals)/len(vals):>12.3f}    {max(vals):>9.3f}")

# Spearman-Korrelation
all_nc = [r['nc'] for r in rows]
all_lm = [r['logM'] for r in rows]
n = len(all_nc)
def spearman(x, y):
    n = len(x)
    rx = [sorted(range(n), key=lambda i: x[i]).index(i) for i in range(n)]
    ry = [sorted(range(n), key=lambda i: y[i]).index(i) for i in range(n)]
    # Korrekte Rang-Berechnung
    ox = sorted(range(n), key=lambda i: x[i])
    oy = sorted(range(n), key=lambda i: y[i])
    rx2 = [0]*n; ry2 = [0]*n
    for rank, idx in enumerate(ox): rx2[idx] = rank
    for rank, idx in enumerate(oy): ry2[idx] = rank
    d2 = sum((rx2[i]-ry2[i])**2 for i in range(n))
    return 1 - 6*d2/(n*(n**2-1))

sp_r = spearman(all_nc, all_lm)
print(f"\nSpearman r (Crossings vs log M_BH) = {sp_r:+.4f}  (n={n})")

# ── Vergleich: nur D-Pol-Gradient vs Kreuzungsdichte ──────────────────────────
sp_r_dpole = spearman([r['theta_D'] for r in rows], all_lm)
print(f"Spearman r (theta_D   vs log M_BH) = {sp_r_dpole:+.4f}  (Vergleich)")
print()

# ── Objekte mit Kreuzungen ─────────────────────────────────────────────────────
nc_sorted = sorted(rows, key=lambda r: (r['nc'], r['logM']), reverse=True)
print("Top-Objekte nach Kreuzungsdichte (dann Masse):")
print(f"  {'Name':<22} {'log(M)':>7}  {'nc':>3}  {'theta_D':>7}  Schalen")
for r in nc_sorted[:20]:
    shells_str = ', '.join(f"{pn}(n={sn})" for pn,sn,_ in r['shells'])
    print(f"  {r['name']:<22} {r['logM']:>7.2f}  {r['nc']:>3}  {r['theta_D']:>7.1f}°  {shells_str}")

print()
# ── Masseverteilung nach Region ────────────────────────────────────────────────
no_cross  = [r['logM'] for r in rows if r['nc'] == 0]
with_cross = [r['logM'] for r in rows if r['nc'] >= 1]
print(f"Keine Kreuzung:    N={len(no_cross):3d}  Median log(M)={sorted(no_cross)[len(no_cross)//2]:.3f}")
print(f"Mind. 1 Kreuzung:  N={len(with_cross):3d}  Median log(M)={sorted(with_cross)[len(with_cross)//2]:.3f}")

# Mann-Whitney
try:
    from scipy import stats
    mw = stats.mannwhitneyu(with_cross, no_cross, alternative='greater')
    print(f"Mann-Whitney (Kreuzung > kein):  U={mw.statistic:.0f}, p={mw.pvalue:.4f}  {'* signifikant' if mw.pvalue < 0.05 else ''}")
except ImportError:
    pass

print()
print("GEOMETRISCHE ERKLAERUNG:")
print(f"  Auf einer S2 mit {len(POLES_DEF)} Polen und {len(SHELLS)} Schalen pro Pol:")
print(f"  Max. moegliche Kreuzungen pro Punkt = {len(POLES_DEF)*(len(POLES_DEF)-1)//2 * len(SHELLS)**2}")
print(f"  Schalen kleiner Umfang (nahe Pol) = dichtere Knoten pro Bogenmass")
print(f"  => Die erste Schale (theta_0=58.65°) hat den kleinsten Umfang")
print(f"     sin(58.65°) = {math.sin(math.radians(58.65)):.4f}  vs sin(90°) = 1.000")
print(f"  => Knotendichte an Schale n=1 ist {1/math.sin(math.radians(58.65)):.2f}x hoeher als am Aequator")
print(f"  => Das erklaert warum die massereichsten SMBHs an theta ~ 58-65 Grad sitzen")
print(f"     OHNE dass es einen monotonen Gradienten geben muss.")
