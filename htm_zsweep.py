"""
HTM Redshift-Sweep: nc=3-Ratio als Funktion von z
==================================================

Strategie: EIN großer SDSS-DR16Q-Datenabruf (z=0.10–2.55),
           lokal in ~20 Scheiben aufteilen (Δz≈0.12),
           pro Scheibe nc=3-Ratio vs. MC rechnen.

Vorhersage HTM:
  n=1-Schale  R=4228 Mpc  z≈0.97  → nc=3-Peak erwartet
  n=2-Schale  R=8456 Mpc  z≈2.30  → nc=3-Erhöhung möglich
  Zwischen Schalen         → Ratio ≈ 1.0

Ausgabe: Tabelle + ASCII-Plot in results/OT_HTM_zsweep.txt
"""
import math, os, sys, re, random, time, ssl, urllib.request, urllib.parse
from collections import Counter

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

_SSL = ssl.create_default_context()
_SSL.check_hostname = False
_SSL.verify_mode    = ssl.CERT_NONE

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(RES, "OT_HTM_zsweep.txt")
os.makedirs(RES, exist_ok=True)

# ── HTM-Konstanten ──────────────────────────────────────────────────────────
THETA0  = 58.65
R_S     = 4228.3
H0      = 73.0
TOL_ANG = 4.0

POLES = [
    (305.0,  25.0),
    (125.0, -25.0),
    ( 35.0,  25.0),
    (215.0, -25.0),
    (215.0,  25.0),
    ( 35.0, -25.0),
]

def radec_to_lb(ra, dec):
    ra, dc = math.radians(ra), math.radians(dec)
    RN, DN = math.radians(192.85948), math.radians(27.12825)
    b = math.asin(max(-1., min(1.,
        math.sin(dc)*math.sin(DN) + math.cos(dc)*math.cos(DN)*math.cos(ra-RN))))
    y = math.cos(dc)*math.sin(ra-RN)
    x = math.cos(dc)*math.sin(DN)*math.cos(ra-RN) - math.sin(dc)*math.cos(DN)
    l = (math.degrees(math.atan2(y, x)) + 122.93192) % 360.0
    return l, math.degrees(b)

def lb_xyz(l, b):
    lr, br = math.radians(l), math.radians(b)
    return (math.cos(br)*math.cos(lr), math.cos(br)*math.sin(lr), math.sin(br))

def ang_sep(v1, v2):
    d = max(-1., min(1., sum(a*b for a, b in zip(v1, v2))))
    return math.degrees(math.acos(d))

def crossing_density(l, b):
    sv = lb_xyz(l, b)
    nc = 0
    for pl, pb in POLES:
        pv = lb_xyz(pl, pb)
        th = ang_sep(pv, sv)
        for n in [1, 2, 3]:
            if abs(th - n*THETA0) < TOL_ANG:
                nc += 1
    return nc

def parse_coord(s):
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        pass
    neg = s.startswith('-')
    parts = re.split(r'[\s:]+', s.lstrip('+-'))
    try:
        v = float(parts[0])
        if len(parts) > 1: v += float(parts[1]) / 60
        if len(parts) > 2: v += float(parts[2]) / 3600
        return -v if neg else v
    except Exception:
        return None

# ── VizieR-Abfrage ──────────────────────────────────────────────────────────
def vizier_query(source, cols, constraints='', max_out=100000, label=''):
    params = {'-source': source, '-out': cols,
              '-out.max': str(max_out), '-oc.form': 'dec'}
    url = ('https://vizier.cds.unistra.fr/viz-bin/votable?'
           + urllib.parse.urlencode(params))
    if constraints:
        url += '&' + constraints
    print(f"  [{label}] {source} z-Abfrage ...", end=' ', flush=True)
    for attempt in range(3):
        try:
            if attempt > 0:
                time.sleep(8 * attempt)
                print(f"Retry{attempt}...", end=' ', flush=True)
            req = urllib.request.Request(
                url, headers={'User-Agent': 'Python/HTMSweep 1.0'})
            with urllib.request.urlopen(req, timeout=300, context=_SSL) as r:
                raw = r.read().decode('utf-8', 'replace')
            rows = []
            for tm in re.finditer(r'<TABLE\b[^>]*>(.*?)</TABLE>', raw, re.DOTALL):
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
            print(f"{len(rows)} Zeilen")
            return rows
        except Exception as e:
            print(f"ERR:{e}")
    return []

# ── Hilfsfunktionen ─────────────────────────────────────────────────────────
def mc_baseline(b_samples, n_mc=5000, seed=42):
    random.seed(seed)
    cnt = Counter()
    for _ in range(n_mc):
        b = random.choice(b_samples)
        l = random.uniform(0, 360)
        cnt[crossing_density(l, b)] += 1
    return cnt, n_mc

def nc3_ratio(lbs, n_mc=5000, seed=42):
    """Berechne nc=3-Ratio für eine Menge von (l,b)-Punkten."""
    if len(lbs) < 20:
        return None, None, None
    nc_obs = Counter(crossing_density(l, b) for l, b in lbs)
    N = len(lbs)
    mc_cnt, mc_n = mc_baseline([b for _, b in lbs], n_mc=n_mc, seed=seed)
    results = {}
    for nc in range(6):
        f_o  = nc_obs.get(nc, 0) / N
        f_mc = mc_cnt.get(nc, 0) / mc_n
        results[nc] = (f_o, f_mc, f_o/f_mc if f_mc > 0 else float('nan'))
    return N, nc_obs, results

# ── Hauptprogramm ───────────────────────────────────────────────────────────
print("=" * 68)
print("HTM REDSHIFT-SWEEP — nc=3 Signal-Profil entlang der Sichtlinie")
print("=" * 68)
print(f"\nHTM: θ₀={THETA0}°, TOL_ANG={TOL_ANG}°, R_S={R_S} Mpc")
print(f"Erwarteter Signal-Peak: z≈0.97 (n=1-Schale, R=4228 Mpc)")
print()

# Eine große Abfrage: SDSS DR16Q z=0.10–2.55 (bis 100k Quasare)
rows = vizier_query(
    'VII/289/dr16q',
    'RAJ2000 DEJ2000 Z',
    constraints='Z=0.10..2.55',
    max_out=100000,
    label='SDSS-DR16Q z=0.10-2.55'
)

if not rows:
    print("FEHLER: Keine SDSS-Daten erhalten. Abbruch.")
    sys.exit(1)

# Koordinaten parsen
ra_key  = next((k for k in rows[0] if 'RA'  in k.upper()), None)
dec_key = next((k for k in rows[0] if k.upper() in ('DECJ2000','DEJ2000','DEC')), None)
if dec_key is None:
    dec_key = next((k for k in rows[0] if 'DE' in k.upper()), None)
z_key   = next((k for k in rows[0] if k.upper() in ('Z','ZSPEC','REDSHIFT')), None)
print(f"  Felder: RA='{ra_key}'  Dec='{dec_key}'  z='{z_key}'")

objects = []   # (z, l, b)
skipped = 0
for r in rows:
    try:
        ra  = parse_coord(r.get(ra_key, ''))
        dec = parse_coord(r.get(dec_key, ''))
        z   = float(r.get(z_key, 'nan'))
        if ra is None or dec is None or not math.isfinite(z) or z <= 0:
            skipped += 1
            continue
        l, b = radec_to_lb(ra, dec)
        objects.append((z, l, b))
    except Exception:
        skipped += 1

print(f"  Verwertbar: {len(objects)} Objekte  (übersprungen: {skipped})")

# ── z-Bänder ────────────────────────────────────────────────────────────────
# 21 Bänder von z=0.15 bis z=2.45 (Δz=0.11, Breite=0.14 je Band)
Z_CENTERS = [round(0.12 + i*0.12, 2) for i in range(21)]   # 0.12 … 2.52
HALF_W    = 0.07   # ±0.07 → Bandbreite 0.14

print(f"\n  Analysiere {len(Z_CENTERS)} Redshift-Bänder (Δz=0.14 je Band)...\n")

# HTM-Schalenredshifts aus R=n×R_S (H0=73)
def z_from_r(r_mpc):
    """Grobe Hubble-Näherung z≈v/c = H0*r/c"""
    v = H0 * r_mpc
    c = 299792.458
    # exakter aus kosmologischer Distanzformel schwer, nutze Hubble-Näherung
    # For R~4000-8000 Mpc Hubble-Approx nicht gut — gebrauche gespeicherte Werte
    return v / c

Z_SHELLS = {1: 0.97, 2: 2.30}  # bekannte n=1,2-Schalenrotverschiebungen

results_sweep = []

for z_c in Z_CENTERS:
    z_lo, z_hi = z_c - HALF_W, z_c + HALF_W
    lbs = [(l, b) for z, l, b in objects if z_lo <= z < z_hi]
    N, nc_obs_raw, nc_dict = nc3_ratio(lbs, n_mc=5000, seed=int(z_c*100))
    if N is None:
        results_sweep.append((z_c, 0, None, None, None))
        continue
    f3_obs, f3_mc, r3 = nc_dict.get(3, (0, 0, float('nan')))
    f4_obs, f4_mc, r4 = nc_dict.get(4, (0, 0, float('nan')))
    results_sweep.append((z_c, N, f3_obs*100, f3_mc*100, r3, f4_obs*100, f4_mc*100, r4))
    label = ''
    for n, z_sh in Z_SHELLS.items():
        if abs(z_c - z_sh) < HALF_W + 0.05:
            label = f'  ← n={n}-Schale'
    r3s = f'{r3:.3f}' if r3 is not None and math.isfinite(r3) else '  --  '
    print(f"  z={z_c:.2f}  N={N:5d}  nc3obs={f3_obs*100:5.2f}%  "
          f"nc3mc={f3_mc*100:5.2f}%  Ratio={r3s}{label}")

print()

# ── Tabelle ─────────────────────────────────────────────────────────────────
lines = []
lines.append("=" * 68)
lines.append("HTM REDSHIFT-SWEEP — Ergebnisse")
lines.append(f"Datum: {time.strftime('%Y-%m-%d %H:%M')}")
lines.append(f"θ₀={THETA0}°  TOL_ANG={TOL_ANG}°  R_S={R_S} Mpc")
lines.append("=" * 68)
lines.append("")
lines.append(f"{'z-Mitte':>7} {'N':>6} {'nc3 obs%':>9} {'nc3 mc%':>8} "
             f"{'Ratio':>7} {'nc4 Ratio':>10}")
lines.append("-" * 60)

for row in results_sweep:
    if len(row) < 8 or row[4] is None:
        lines.append(f"  {row[0]:.2f}  {'N<20':>6}")
        continue
    z_c, N, f3o, f3m, r3, f4o, f4m, r4 = row
    note = ''
    for n, z_sh in Z_SHELLS.items():
        if abs(z_c - z_sh) <= HALF_W:
            note = f'  ← n={n}'
    r3s = f'{r3:.3f}' if math.isfinite(r3) else ' ---'
    r4s = f'{r4:.3f}' if math.isfinite(r4) else ' ---'
    lines.append(f"  {z_c:.2f}  {N:6d}  {f3o:8.2f}%  {f3m:7.2f}%  {r3s:>7}  {r4s:>10}{note}")

lines.append("")

# ── ASCII-Plot ──────────────────────────────────────────────────────────────
lines.append("── nc=3-Ratio Profil (ASCII) ──")
lines.append("")

valid = [(row[0], row[4]) for row in results_sweep
         if len(row) >= 5 and row[4] is not None and math.isfinite(row[4]) and row[1] >= 40]

if valid:
    max_r = max(r for _, r in valid)
    min_r = min(r for _, r in valid)
    PLOT_H = 15
    PLOT_W = len(valid)

    y_max = max(max_r * 1.1, 2.5)
    y_min = max(min_r * 0.9, 0.0)

    lines.append(f"  Ratio")
    for row_i in range(PLOT_H, -1, -1):
        y_val = y_min + (y_max - y_min) * row_i / PLOT_H
        bar = f"  {y_val:4.2f} |"
        for z_c, r in valid:
            cell_y = (r - y_min) / (y_max - y_min) * PLOT_H
            diff = abs(cell_y - row_i)
            if diff < 0.55:
                bar += '*'
            else:
                bar += ' '
        lines.append(bar)

    # x-Achse
    lines.append("       " + "-" * (PLOT_W + 2))
    x_tick = "       z:"
    for z_c, _ in valid:
        if round(z_c * 10) % 2 == 0:
            x_tick += f"{z_c:.1f}"
        else:
            x_tick += "   "
    lines.append(x_tick)

    # Schalen markieren
    shell_line = "         "
    for z_c, r in valid:
        for n, z_sh in Z_SHELLS.items():
            if abs(z_c - z_sh) <= HALF_W:
                shell_line += str(n)
                break
        else:
            shell_line += " "
    lines.append(f"  n-Schale:{shell_line[9:]}")
    lines.append("")
    lines.append("  *= nc=3-Ratio-Wert  n=Schalen-Index")
    lines.append("")

# ── Ergebniszusammenfassung ─────────────────────────────────────────────────
lines.append("── Zusammenfassung ──")
lines.append("")

if valid:
    z_peak, r_peak = max(valid, key=lambda x: x[1])
    lines.append(f"  Signal-Peak:  z={z_peak:.2f}  nc=3-Ratio={r_peak:.3f}")
    for n, z_sh in Z_SHELLS.items():
        delta_z = abs(z_peak - z_sh)
        lines.append(f"    Distanz zur n={n}-Schale (z={z_sh}): Δz={delta_z:.2f}")

    lines.append("")
    # Werte an den Schalen
    for n, z_sh in Z_SHELLS.items():
        near = [(z_c, r) for z_c, r in valid if abs(z_c - z_sh) <= HALF_W + 0.01]
        if near:
            z_c, r = min(near, key=lambda x: abs(x[0]-z_sh))
            lines.append(f"  n={n}-Schale z={z_sh}: nc=3-Ratio={r:.3f}")
        else:
            lines.append(f"  n={n}-Schale z={z_sh}: keine Daten im Fenster")

lines.append("")
lines.append("Vollständige Ergebnisse → " + OUT)

# ── Ausgabe ─────────────────────────────────────────────────────────────────
report = '\n'.join(lines)
print(report)

with open(OUT, 'w', encoding='utf-8') as f:
    f.write(report + '\n')

print(f"\nFERTIG. Ergebnisse gespeichert: {OUT}")
