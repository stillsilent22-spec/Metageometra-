"""
HTM-Knoten-Signal: Unabhängige Bestätigung
==========================================
Strategie wie Einsteins Merkur-Test:
  Bekanntes Signal (SDSS nc=3, 1.74×) → teste mit VÖLLIG ANDEREN Daten

DREI Tests:
  A) 2QZ  — 2dF QSO Survey (südl. Himmel, andere Sphärenabdeckung als SDSS)
             VizieR VII/241/2qz, z=0.90-1.04
  B) Null-Test z≈0.50  (Kontrolle: kein HTM-Signal erwartet)
     → SDSS DR16Q bei z=0.46-0.54 (zwischen den Schalen)
  C) Nächste Schale n=2 bei z≈2.3  (R_Thales ≈ 8456 Mpc)
     → SDSS DR16Q bei z=2.25-2.35

Vorhersage:
  A) 2QZ (unabhängiger Himmel) → nc=3 Überschuss wie SDSS → Bestätigung
  B) z≈0.50 (keine HTM-Schale) → nc=3 ≈ MC-Erwartung (Null-Test)
  C) z≈2.30 (n=2-Schale) → nc=3 Überschuss (zweite Schale)

HTM: θ₀=58.65°, R_S=4228.3 Mpc, TOL_ANG=4°
"""
import math, os, sys, re, random, time, ssl, urllib.request, urllib.parse
import xml.etree.ElementTree as ET
from collections import Counter

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# SSL ohne Verifikation (für VizieR/CDS)
_SSL = ssl.create_default_context()
_SSL.check_hostname = False
_SSL.verify_mode    = ssl.CERT_NONE

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(RES, "OT_HTM_indep.txt")
os.makedirs(RES, exist_ok=True)

# ── HTM-Konstanten ─────────────────────────────────────────────────────────
THETA0  = 58.65
R_S     = 4228.3   # Mpc
H0      = 73.0     # km/s/Mpc
TOL_ANG = 4.0      # °

POLES_DEF = [
    (305.0,  25.0, 'D'),
    (125.0, -25.0, 'A'),
    ( 35.0,  25.0, 'S1'),
    (215.0, -25.0, 'S2'),
    (215.0,  25.0, 'S3'),
    ( 35.0, -25.0, 'S4'),
]

# Vorhersagen
PRED_Z = {
    'z0.97_n1': (0.90, 1.04, True,  "n=1-Schale (Dualitätssphäre, Referenzsignal)"),
    'z0.50_null':(0.46, 0.54, False, "z≈0.50 (Null-Kontrolle, keine Schale)"),
    'z2.30_n2': (2.23, 2.37, True,  "n=2-Schale (z≈2.30, R≈8456 Mpc)"),
}

# ── Koordinaten ────────────────────────────────────────────────────────────
def radec_to_lb(ra, dec):
    ra, dc = math.radians(ra), math.radians(dec)
    RN, DN  = math.radians(192.85948), math.radians(27.12825)
    b  = math.asin(max(-1., min(1.,
         math.sin(dc)*math.sin(DN) + math.cos(dc)*math.cos(DN)*math.cos(ra-RN))))
    y  = math.cos(dc)*math.sin(ra-RN)
    x  = math.cos(dc)*math.sin(DN)*math.cos(ra-RN) - math.sin(dc)*math.cos(DN)
    l  = (math.degrees(math.atan2(y,x)) + 122.93192) % 360.0
    return l, math.degrees(b)

def lb_xyz(l, b):
    l, b = math.radians(l), math.radians(b)
    return (math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b))

def ang_sep(v1, v2):
    d = sum(a*b for a,b in zip(v1, v2))
    return math.degrees(math.acos(max(-1., min(1., d))))

def crossing_density(l, b):
    sv = lb_xyz(l, b)
    nc = 0
    for pl, pb, _ in POLES_DEF:
        pv = lb_xyz(pl, pb)
        th = ang_sep(pv, sv)
        for n in [1, 2, 3]:
            if abs(th - n*THETA0) < TOL_ANG:
                nc += 1
    return nc

# ── VizieR VOTable-Abfrage ─────────────────────────────────────────────────
def vizier_query(source, cols, constraints='', max_out=50000, label=''):
    params = {'-source': source, '-out': cols,
              '-out.max': str(max_out), '-oc.form': 'dec'}
    url = 'https://vizier.cds.unistra.fr/viz-bin/votable?' + urllib.parse.urlencode(params)
    if constraints:
        url += '&' + constraints
    print(f"  [{label}] {source} ...", end=' ', flush=True)
    for attempt in range(3):
        try:
            if attempt > 0:
                time.sleep(6*attempt)
                print(f"Retry{attempt}...", end=' ', flush=True)
            req = urllib.request.Request(url,
                headers={'User-Agent': 'Python/HTMVerify 1.0'})
            with urllib.request.urlopen(req, timeout=180, context=_SSL) as r:
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

# ── Monte-Carlo-Baseline ────────────────────────────────────────────────────
def mc_baseline(b_samples, n_mc=20000, seed=42):
    random.seed(seed)
    cnt = Counter()
    for _ in range(n_mc):
        b = random.choice(b_samples)
        l = random.uniform(0, 360)
        cnt[crossing_density(l, b)] += 1
    return cnt, n_mc

# ── Koordinatenparser ─────────────────────────────────────────────────────
def parse_coord(s):
    """Parst Dezimal-Grad-String ODER Sexagesimal (h:m:s oder d:m:s) → float"""
    s = s.strip()
    if not s:
        return None
    # Versuche zuerst direkt als float
    try:
        return float(s)
    except ValueError:
        pass
    # Sexagesimal: "HH MM SS.s" oder "HH:MM:SS.s" oder "-DD:MM:SS.s" etc.
    neg = s.startswith('-')
    s2  = s.lstrip('+-')
    parts = re.split(r'[\s:]+', s2)
    try:
        v = float(parts[0])
        if len(parts) > 1: v += float(parts[1])/60
        if len(parts) > 2: v += float(parts[2])/3600
        if neg: v = -v
        return v
    except Exception:
        return None

# ── Analyse einer Objektliste ───────────────────────────────────────────────
def analyse(rows, ra_key, dec_key, label, n_mc=20000, seed=42, ra_scale=1.0):
    """ra_scale=15.0 fuer RA in Stunden (2QZ liefert h:m:s)"""
    lbs = []
    for r in rows:
        try:
            ra  = parse_coord(r.get(ra_key, ''))
            dec = parse_coord(r.get(dec_key, ''))
            if ra is None or dec is None:
                continue
            ra  *= ra_scale   # Stunden → Grad falls nötig
            lbs.append(radec_to_lb(ra, dec))
        except Exception:
            continue
    if not lbs:
        return None
    nc_obs  = Counter(crossing_density(l, b) for l, b in lbs)
    bs      = [b for _, b in lbs]
    mc_cnt, mc_n = mc_baseline(bs, n_mc=n_mc, seed=seed)
    N = len(lbs)

    results = {}
    for nc in sorted(set(nc_obs.keys()) | set(mc_cnt.keys())):
        f_o  = nc_obs.get(nc, 0) / N
        f_mc = mc_cnt.get(nc, 0) / mc_n
        ratio = f_o/f_mc if f_mc > 0 else float('nan')
        results[nc] = {'obs': nc_obs.get(nc, 0), 'frac_obs': f_o,
                       'frac_mc': f_mc, 'ratio': ratio}
    return {'label': label, 'N': N, 'nc': results,
            'lbs': lbs, 'mc_cnt': mc_cnt, 'mc_n': mc_n}

def print_result(res):
    if res is None:
        print("    -> KEINE DATEN")
        return
    print(f"    N={res['N']}  nc-Verteilung:")
    print(f"    {'nc':>4}  {'obs%':>7}  {'MC%':>7}  {'Ratio':>7}  Signal")
    for nc in sorted(res['nc']):
        d = res['nc'][nc]
        flag = ''
        if d['ratio'] > 1.2 and d['frac_obs']*100 > 2:
            flag = '★ SIGNAL'
        elif d['ratio'] < 0.8 and d['frac_obs']*100 > 2:
            flag = '▼ DEFICIT'
        print(f"    {nc:4d}  {100*d['frac_obs']:7.2f}  {100*d['frac_mc']:7.2f}  "
              f"{d['ratio']:7.3f}  {flag}")

# ═══════════════════════════════════════════════════════════════════════════
print("=" * 68)
print("HTM-Quasar-Signal — UNABHÄNGIGE BESTÄTIGUNG")
print("(analog zu Einsteins Merkur-Test: Effekt mit anderen Daten prüfen)")
print("=" * 68)
print(f"\nHTM: θ₀={THETA0}°, TOL_ANG={TOL_ANG}°")
print(f"Referenz:  SDSS DR16Q z≈0.97 → nc=3 Ratio=1.74× (bereits bekannt)")
print()

all_results = {}

# ── TEST A: 2QZ (südlicher Himmel) ─────────────────────────────────────────
print("─" * 68)
print("TEST A: 2QZ — 2dF QSO Redshift Survey (South, VII/241)")
print("  Südlicher galakt. Pol (SGP), unabhängig von SDSS-Footprint")

rows_2qz = vizier_query(
    'VII/241/2qz', 'RAJ2000 DEJ2000 z1',
    constraints='z1=0.90..1.04',
    max_out=50000, label='2QZ z=0.90-1.04')

if not rows_2qz:
    print("  2QZ nicht verfügbar, versuche eBOSS DR14 (VII/280)...")
    rows_2qz = vizier_query(
        'VII/280/dr14qso', 'RAJ2000 DEJ2000 Z',
        constraints='Z=0.90..1.04',
        max_out=50000, label='eBOSS-DR14 z=0.90-1.04')

if rows_2qz:
    sample = rows_2qz[0]
    ra_k = next((k for k in sample if 'RA' in k.upper()), None)
    de_k = next((k for k in sample if k.upper() in ('DEJ2000','DEC','DE','DECJ2000')), None)
    if de_k is None:
        de_k = next((k for k in sample if 'DE' in k.upper() or 'DEC' in k.upper()), None)
    print(f"  Felder erkannt: RA='{ra_k}', Dec='{de_k}'")
    print(f"  Beispiel: RA={sample.get(ra_k,'?')}  Dec={sample.get(de_k,'?')}")
    # Prüfe ob RA in Stunden (h:m:s) oder Dezimal-Grad
    ra_test = parse_coord(sample.get(ra_k, '0'))
    ra_scale = 15.0 if ra_test is not None and ra_test < 25 else 1.0
    print(f"  RA-Muster: {'Stunden → ×15' if ra_scale == 15.0 else 'Dezimal-Grad'}")
    # Filtere nur QSOs (2QZ enthält auch Stars)
    id_key = next((k for k in sample if k.lower() in ('id1','id','class')), None)
    if id_key:
        before = len(rows_2qz)
        rows_2qz = [r for r in rows_2qz if 'QSO' in r.get(id_key,'').upper()]
        print(f"  QSO-Filter (ID={id_key}): {before} → {len(rows_2qz)}")
    res_A = analyse(rows_2qz, ra_k, de_k, 'TEST-A (2QZ/eBOSS)',
                    seed=101, ra_scale=ra_scale)
    print_result(res_A)
    all_results['A'] = res_A
else:
    print("  KEIN DATENSATZ GEFUNDEN – überspringe Test A")
    all_results['A'] = None

# ── TEST B: Null-Kontrolle z≈0.50 ──────────────────────────────────────────
print("\n" + "─" * 68)
print("TEST B: Null-Kontrolle — SDSS DR16Q bei z≈0.50")
print("  Keine HTM-Schale bei 2140 Mpc → Ratio ≈ 1.0 erwartet")

rows_null = vizier_query(
    'VII/289/dr16q', 'RAJ2000 DEJ2000 zsp',
    constraints='zsp=0.46..0.54',
    max_out=50000, label='SDSS-DR16Q z=0.46-0.54')

if rows_null:
    ra_k = next((k for k in rows_null[0] if 'RA' in k.upper()), None)
    de_k = next((k for k in rows_null[0] if 'DE' in k.upper() or 'DEC' in k.upper()), None)
    res_B = analyse(rows_null, ra_k, de_k, 'TEST-B (z≈0.50, Null)', seed=202)
    print_result(res_B)
    all_results['B'] = res_B
else:
    print("  Keine Daten")
    all_results['B'] = None

# ── TEST C: Zweite HTM-Schale n=2 bei z≈2.30 ───────────────────────────────
print("\n" + "─" * 68)
print("TEST C: Zweite Schale n=2 — SDSS DR16Q bei z≈2.30")
print(f"  R_Thales(n=2) = 2×{R_S:.0f} = {2*R_S:.0f} Mpc → z≈2.30")

rows_n2 = vizier_query(
    'VII/289/dr16q', 'RAJ2000 DEJ2000 zsp',
    constraints='zsp=2.23..2.37',
    max_out=50000, label='SDSS-DR16Q z=2.23-2.37')

if rows_n2:
    ra_k = next((k for k in rows_n2[0] if 'RA' in k.upper()), None)
    de_k = next((k for k in rows_n2[0] if 'DE' in k.upper() or 'DEC' in k.upper()), None)
    res_C = analyse(rows_n2, ra_k, de_k, 'TEST-C (z≈2.30, n=2-Schale)', seed=303)
    print_result(res_C)
    all_results['C'] = res_C
else:
    print("  Keine Daten")
    all_results['C'] = None

# ── REFERENZ: SDSS DR16Q z≈0.97 (für Vergleich) ────────────────────────────
print("\n" + "─" * 68)
print("REFERENZ: SDSS DR16Q bei z≈0.97 (Dualitätssphäre n=1, 100k Quasare)")
rows_ref = vizier_query(
    'VII/289/dr16q', 'RAJ2000 DEJ2000 zsp',
    constraints='zsp=0.90..1.04',
    max_out=100000, label='SDSS-DR16Q z=0.90-1.04 (Referenz)')

if rows_ref:
    ra_k = next((k for k in rows_ref[0] if 'RA' in k.upper()), None)
    de_k = next((k for k in rows_ref[0] if 'DE' in k.upper() or 'DEC' in k.upper()), None)
    res_R = analyse(rows_ref, ra_k, de_k, 'REFERENZ (z≈0.97)', seed=42)
    print_result(res_R)
    all_results['R'] = res_R
else:
    all_results['R'] = None

# ── GESAMTAUSWERTUNG ────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("GESAMTAUSWERTUNG — nc=3 Ratio-Vergleich")
print("=" * 68)
print()
print(f"  {'Test':<28}  {'N':>6}  {'nc=3 Obs%':>10}  {'nc=3 MC%':>9}  {'Ratio':>7}  Erwartung")
print(f"  {'-'*28}  {'-'*6}  {'-'*10}  {'-'*9}  {'-'*7}  {'-'*12}")

def nc3_line(key, label, expect_signal):
    res = all_results.get(key)
    if res is None:
        print(f"  {label:<28}  {'N/A':>6}")
        return
    N   = res['N']
    d3  = res['nc'].get(3, {'frac_obs': 0, 'frac_mc': 0, 'ratio': float('nan')})
    obs = 100*d3['frac_obs']
    mc  = 100*d3['frac_mc']
    ratio = d3['ratio']
    expect = "Signal erwartet" if expect_signal else "Null erwartet"
    flag = ''
    if expect_signal and ratio > 1.2:    flag = ' ✓ Bestätigt'
    elif expect_signal and ratio < 0.9:  flag = ' ✗ Kein Signal'
    elif not expect_signal and ratio < 1.1: flag = ' ✓ Null best.'
    elif not expect_signal and ratio > 1.2: flag = ' ! Unexpected'
    print(f"  {label:<28}  {N:>6}  {obs:>10.3f}  {mc:>9.3f}  {ratio:>7.3f}  {expect}{flag}")

nc3_line('R', 'SDSS z=0.90-1.04 (Ref.)',    True)
nc3_line('A', '2QZ/eBOSS z=0.90-1.04',      True)
nc3_line('B', 'SDSS z=0.46-0.54 (Null)',    False)
nc3_line('C', 'SDSS z=2.23-2.37 (n=2)',     True)

print()
print("Vorhersage-Schema:")
print("  - Ref (z≈0.97):  nc=3 Ratio > 1.5  (bereits bekannt: 1.74×)")
print("  - A (2QZ, unabh.)   nc=3 Ratio > 1.5  → Einstein-Test")
print("  - B (Null, z≈0.50): nc=3 Ratio ≈ 1.0  → Selbstkontrolle")
print("  - C (n=2, z≈2.30):  nc=3 Ratio > 1.5  → Zweite Schale")

# ── nc=1 und nc=2 Vergleich ────────────────────────────────────────────────
print("\n── Vollständiger nc-Vergleich (alle Tests) ──")
for nc in [1, 2, 3, 4]:
    print(f"\n  nc={nc}:")
    for key, lbl in [('R','Ref z=0.97'), ('A','A z=0.97 2QZ'),
                     ('B','B z=0.50 Null'), ('C','C z=2.30 n=2')]:
        res = all_results.get(key)
        if not res:
            continue
        d  = res['nc'].get(nc, {'frac_obs': 0, 'frac_mc': 0, 'ratio': float('nan')})
        print(f"    {lbl:<18}  {100*d['frac_obs']:6.2f}%  MC={100*d['frac_mc']:.2f}%  "
              f"Ratio={d['ratio']:.3f}")

# ── Ausgabe in Datei ────────────────────────────────────────────────────────
import io
buf = io.StringIO()
old_stdout = sys.stdout
sys.stdout = buf

print("=" * 68)
print("HTM-Quasar-Signal: Unabhängige Bestätigung")
print(f"Datum: {time.strftime('%Y-%m-%d %H:%M')}")
print(f"θ₀={THETA0}°  TOL_ANG={TOL_ANG}°  R_S={R_S} Mpc")
print(f"Referenzsignal: SDSS DR16Q z≈0.97, nc=3 Ratio=1.74× (vorher)")
print("=" * 68)
print()
print("Vorhersage (wie Merkur-Perihel-Test):")
print("  HTM-Schale n=1 bei R=4228 Mpc (z≈0.97): nc=3-Exzess erwartet")
print("  HTM-Schale n=2 bei R=8456 Mpc (z≈2.30): nc=3-Exzess erwartet")
print("  Zwischen Schalen (z≈0.50): kein Exzess erwartet")
print()

for key, lbl, expect in [('R','SDSS DR16Q z=0.97 (Referenz)', True),
                          ('A','Test A: 2QZ/eBOSS z=0.97', True),
                          ('B','Test B: SDSS z=0.50 (Null)', False),
                          ('C','Test C: SDSS z=2.30 (n=2-Schale)', True)]:
    res = all_results.get(key)
    print(f"── {lbl} ──")
    if res is None:
        print("   KEIN DATENSATZ")
        continue
    print(f"   N={res['N']}")
    print(f"   {'nc':>4}  {'obs%':>7}  {'MC%':>7}  {'Ratio':>7}")
    for nc in sorted(res['nc'].keys()):
        d = res['nc'][nc]
        print(f"   {nc:4d}  {100*d['frac_obs']:7.3f}  {100*d['frac_mc']:7.3f}  {d['ratio']:7.4f}")
    print()

print()
print("── nc=3 Zusammenfassung ──")
print(f"  {'Test':<30}  {'Ratio':>7}  Erwartung      Ergebnis")
for key, lbl, expect_signal in [
        ('R','SDSS z=0.97 (Referenz)', True),
        ('A','2QZ/eBOSS z=0.97',      True),
        ('B','SDSS z=0.50 (Null)',    False),
        ('C','SDSS z=2.30 (n=2)',     True)]:
    res = all_results.get(key)
    d3 = res['nc'].get(3, {'ratio': float('nan')}) if res else {'ratio': float('nan')}
    r  = d3['ratio']
    if res is None:
        result = 'N/A'
    elif expect_signal and r > 1.2:
        result = '✓ Bestätigt'
    elif expect_signal:
        result = '✗ Kein Signal'
    elif not expect_signal and r < 1.1:
        result = '✓ Null bestätigt'
    else:
        result = '! Unerwartet'
    expect_lbl = 'Exzess>1.2×' if expect_signal else 'Ratio≈1.0'
    print(f"  {lbl:<30}  {r:>7.3f}  {expect_lbl:<14}  {result}")

sys.stdout = old_stdout
txt = buf.getvalue()
with open(OUT, 'w', encoding='utf-8') as fh:
    fh.write(txt)
print(f"\nVollständige Ergebnisse → {OUT}")
print("FERTIG.")
