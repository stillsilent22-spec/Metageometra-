"""
SDSS-Quasar HTM-Knoten-Test
============================
Ladet SDSS DR16 Quasare bei z=0.90-1.04 (Dualitätsphäre z≈0.97)
und testet ob sie bevorzugt an HTM-Winkelknoten sitzen.

Katalog: VII/289 (SDSS DR16Q, Lyke+2020) via VizieR
"""
import math, os, sys, csv, urllib.request, urllib.parse, time
import xml.etree.ElementTree as ET
import statistics as st
import re
import random
from collections import Counter

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(RES, "OT_QSO_htm.txt")
os.makedirs(RES, exist_ok=True)

# ── HTM-Konstanten ──────────────────────────────────────────────────────────
THETA0   = 58.65
R_S      = 4228.3
Z_TARGET = 0.97
Z_WIN    = 0.07   # z = 0.90 … 1.04
TOL_ANG  = 4.0   # verschärft: 8° wäre ~90% Himmelsabdeckung — zu breit

POLES_DEF = [
    (305.0,  25.0, 'D'),
    (125.0, -25.0, 'A'),
    ( 35.0,  25.0, 'S1'),
    (215.0, -25.0, 'S2'),
    (215.0,  25.0, 'S3'),
    ( 35.0, -25.0, 'S4'),
]

# ── Koordinaten ─────────────────────────────────────────────────────────────
def radec_to_lb(ra, dec):
    ra = math.radians(ra); dc = math.radians(dec)
    RN = math.radians(192.85948); DN = math.radians(27.12825)
    b = math.asin(max(-1., min(1.,
        math.sin(dc)*math.sin(DN) + math.cos(dc)*math.cos(DN)*math.cos(ra-RN))))
    y = math.cos(dc)*math.sin(ra-RN)
    x = math.cos(dc)*math.sin(DN)*math.cos(ra-RN) - math.sin(dc)*math.cos(DN)
    l = (math.degrees(math.atan2(y,x)) + 122.93192) % 360.0
    return l, math.degrees(b)

def lb_xyz(l, b):
    l, b = math.radians(l), math.radians(b)
    return (math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b))

def ang_sep_xyz(v1, v2):
    d = sum(a*b for a,b in zip(v1, v2))
    return math.degrees(math.acos(max(-1., min(1., d))))

def crossing_density(l, b):
    sv = lb_xyz(l, b)
    nc = 0
    for pl, pb, _ in POLES_DEF:
        pv = lb_xyz(pl, pb)
        theta = ang_sep_xyz(pv, sv)
        for n in [1, 2, 3]:
            if abs(theta - n*THETA0) < TOL_ANG:
                nc += 1
    return nc

# ── VizieR-Query ─────────────────────────────────────────────────────────────
def vizier_query_allsky(source, col_str, label, constraints='', max_out=100000):
    params = {
        '-source': source,
        '-out': col_str,
        '-out.max': str(max_out),
        '-oc.form': 'dec',
    }
    if constraints:
        params['-c'] = constraints if not constraints.startswith('-') else ''
        # Use constraint fields directly
    url = 'https://vizier.cds.unistra.fr/viz-bin/votable?' + urllib.parse.urlencode(params)
    if constraints:
        url += '&' + constraints
    print(f"  Query {label}...", end=' ', flush=True)
    for attempt in range(3):
        try:
            if attempt > 0:
                time.sleep(5*attempt)
                print(f"  Retry {attempt}...", end=' ', flush=True)
            req = urllib.request.Request(url, headers={'User-Agent': 'Python/HTMQuasar'})
            with urllib.request.urlopen(req, timeout=180) as r:
                raw = r.read().decode('utf-8', 'replace')
            rows = []
            for tm in re.finditer(r'<TABLE\b[^>]*>(.+?)</TABLE>', raw, re.DOTALL):
                tb = tm.group(1)
                fields = re.findall(r'<FIELD[^>]+name="([^"]+)"', tb)
                fields = [f for f in fields if f not in ('recno', '_ivoid')]
                if not fields:
                    continue
                for tr in re.finditer(r'<TR>(.*?)</TR>', tb, re.DOTALL):
                    tds = re.findall(r'<TD>(.*?)</TD>', tr.group(1))
                    if len(tds) < 2:
                        continue
                    rows.append({fields[i]: tds[i].strip() if i < len(tds) else ''
                                 for i in range(len(fields))})
            if rows or attempt == 2:
                print(f"{len(rows)} Zeilen")
                return rows
            print("0 (Retry)...", end=' ', flush=True)
        except Exception as e:
            print(f"FEHLER: {e}")
            if attempt == 2:
                return []
    return []

# ── Laden ─────────────────────────────────────────────────────────────────────
print("=" * 60)
print(f"SDSS DR16 Quasare bei z = {Z_TARGET} ± {Z_WIN}")
print(f"Dualitätssphäre R_S = {R_S} Mpc  (z ≈ 0.97)")
print("=" * 60)

# SDSS DR16Q: VII/289, Felder: RAJ2000, DEJ2000, zspec
# Filter z über VizieR constraint
z_lo = Z_TARGET - Z_WIN
z_hi = Z_TARGET + Z_WIN
constraint = f'zsp={z_lo:.2f}..{z_hi:.2f}'

rows = vizier_query_allsky(
    'VII/289/dr16q',
    'RAJ2000 DEJ2000 zsp',
    f'SDSS DR16Q z={z_lo:.2f}-{z_hi:.2f}',
    constraints=constraint
)

if not rows:
    # Fallback: eBOSS DR14 Quasorkatalog VII/280
    print("  Fallback: eBOSS DR14Q VII/280...")
    rows = vizier_query_allsky(
        'VII/280/dr14qso',
        'RAJ2000 DEJ2000 Z',
        f'eBOSS DR14Q z={z_lo:.2f}-{z_hi:.2f}',
        constraints=f'Z={z_lo:.2f}..{z_hi:.2f}'
    )

print(f"\nGeladene Quasare: {len(rows)}")

if not rows:
    print("FEHLER: Keine Daten erhalten.")
    sys.exit(1)

# ── Crossing-Density für jeden Quasar ────────────────────────────────────────
print("Berechne HTM-Crossing-Density...", flush=True)

ra_key  = next((k for k in rows[0] if 'RA' in k.upper()), None)
dec_key = next((k for k in rows[0] if 'DE' in k.upper() or 'DEC' in k.upper()), None)
z_key   = next((k for k in rows[0] if k.lower() in ('zsp','z','zspec','z_vi')), None)

print(f"  Spalten: RA={ra_key}, Dec={dec_key}, z={z_key}")

qsos = []
for r in rows:
    try:
        ra  = float(r[ra_key])
        dec = float(r[dec_key])
        z   = float(r[z_key]) if z_key else Z_TARGET
        l, b = radec_to_lb(ra, dec)
        nc   = crossing_density(l, b)
        qsos.append({'l': l, 'b': b, 'z': z, 'nc': nc})
    except Exception:
        continue

print(f"  {len(qsos)} Quasare verarbeitet")

# ── Zufalls-MC-Vergleich ─────────────────────────────────────────────────────
# Zufällige Punkte mit gleicher b-Verteilung wie SDSS-Quasare → Erwartungswert nc
random.seed(42)
bs_real = [q['b'] for q in qsos]
n_mc = min(20000, len(qsos))
mc_nc = []
for _ in range(n_mc):
    b_mc = random.choice(bs_real)      # gleiche galaktische-Breiten-Verteilung
    l_mc = random.uniform(0, 360)
    mc_nc.append(crossing_density(l_mc, b_mc))
mc_counts = Counter(mc_nc)
mc_total  = len(mc_nc)

# ── Statistik ─────────────────────────────────────────────────────────────────
nc_counts = Counter(q['nc'] for q in qsos)
total = len(qsos)

lines = []
lines.append("=" * 60)
lines.append(f"SDSS-Quasar HTM-Knoten-Test  (z={z_lo:.2f}–{z_hi:.2f})")
lines.append(f"N = {total}  Quasare")
lines.append("=" * 60)
lines.append("\n── nc-Verteilung ──")
for nc in sorted(nc_counts):
    frac = 100 * nc_counts[nc] / total
    bar  = '█' * int(frac / 0.5)
    lines.append(f"  nc={nc}:  {nc_counts[nc]:6d}  ({frac:5.2f}%)  {bar}")

n_node = sum(nc_counts[nc] for nc in nc_counts if nc >= 1)
f_node = 100 * n_node / total
lines.append(f"\n  nc≥1 (an HTM-Knoten):  {n_node}  ({f_node:.2f}%)")
lines.append(f"  nc=0 (kein Knoten):    {total-n_node}  ({100-f_node:.2f}%)")

# ── Zufallserwartung (analytisch) ────────────────────────────────────────────
# Jeder Pol hat 3 Schalenringe, jeder Ring hat Breite 2*TOL_ANG auf der Sphäre
# Flächenanteil eines Rings: ΔΩ/4π = sin(θ)·Δθ (Δθ in rad)
delta_theta_rad = math.radians(2 * TOL_ANG)
p_hit_one_ring  = 0.0
for pl, pb, _ in POLES_DEF:
    for n in [1, 2, 3]:
        theta = n * THETA0
        p_ring = math.sin(math.radians(theta)) * delta_theta_rad / 2.0
        p_hit_one_ring += p_ring
# p_hit_one_ring = Gesamtflächenanteil aller 18 Ringe
# P(nc≥1) ≈ 1 - (1-p_ring)^18  (approximiert)
p_nc0_rand = 1.0
for pl, pb, _ in POLES_DEF:
    for n in [1, 2, 3]:
        theta = n * THETA0
        p_ring = math.sin(math.radians(theta)) * delta_theta_rad / 2.0
        p_nc0_rand *= (1 - p_ring)
p_nc1_rand = 1 - p_nc0_rand

lines.append(f"\n── Zufallserwartung (MC, N={n_mc}) ──")
lines.append(f"  {'nc':>4}   Beobachtet     MC-Zufall    Faktor")
all_ncs = sorted(set(list(nc_counts.keys()) + list(mc_counts.keys())))
for nc in all_ncs:
    obs_f  = 100 * nc_counts.get(nc, 0) / total
    mc_f   = 100 * mc_counts.get(nc, 0) / mc_total
    fac    = obs_f / mc_f if mc_f > 0 else 0
    flag   = '  <<<' if fac > 1.5 else ('  >>>' if fac < 0.67 else '')
    lines.append(f"  nc={nc}: {obs_f:7.2f}%   {mc_f:7.2f}%   {fac:.2f}x{flag}")
lines.append(f"\n  nc≥1 beobachtet: {f_node:.2f}%   MC: {100*(n_mc-mc_counts.get(0,0))/n_mc:.2f}%")

# ── Richtungstest: lokale SMBH-Knoten in SDSS-Quasaren ───────────────────────
# Bekannte HTM-Knoten-Richtungen aus node_geometry Analyse
# NGC 4278: RA=185.03°, Dec=29.51° → l≈ 142.7°, b=83.0°
# Perseus (NGC1277): RA=49.82°, Dec=41.57° → l≈150.6°, b=-13.4°
key_dirs = [
    ('NGC4278 (nc=6)',     142.7,  83.0),
    ('Perseus/NGC1277',   150.6, -13.4),
    ('D-Pol',             305.0,  25.0),
    ('S3-Pol',            215.0,  25.0),
    ('NGC4889 (nc=10)',    84.2,  86.7),
]

lines.append(f"\n── Quasardichte in bekannten SMBH-Knotenrichtungen (±15°) ──")
for name, ld, bd in key_dirs:
    nearby = [q for q in qsos if
              math.degrees(math.acos(max(-1., min(1.,
                  math.sin(math.radians(q['b']))*math.sin(math.radians(bd)) +
                  math.cos(math.radians(q['b']))*math.cos(math.radians(bd))*
                  math.cos(math.radians(q['l']-ld)))))) < 15.0]
    omega_cap = 2*math.pi*(1 - math.cos(math.radians(15.0)))
    expected = total * omega_cap / (4*math.pi)
    if nearby:
        nc1_near = sum(1 for q in nearby if q['nc'] >= 1)
        lines.append(f"  {name:<22}  N={len(nearby):5d}  erw.={expected:.0f}  "
                     f"nc≥1={nc1_near}({100*nc1_near/len(nearby):.1f}%)")
    else:
        lines.append(f"  {name:<22}  N=0")

lines.append("\n" + "=" * 60)

# ── Ausgabe ───────────────────────────────────────────────────────────────────
for l in lines:
    print(l)

with open(OUT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print(f"\nGespeichert: {OUT}")
