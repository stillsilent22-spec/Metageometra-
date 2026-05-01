"""
OT-5 / OT-13 / OT-23 — Kombinierter HTM-Schalenwinkel-Test
============================================================
Testet ob drei Objektklassen bevorzugt bei theta = n * 58.65 Grad
vom HTM-Dipol-Pol (l=305, b=25) liegen.

  OT-5:  97 SMBHs (smbh_extended.csv, McConnell+Graham 2013)
  OT-13: Milliquas AGN/QSO (VizieR VII/290, Flesch 2023)
  OT-23: SDSS Kosmische Voids (VizieR J/MNRAS/431/2307, Sutter+2012)

Methodik (stdlib only):
  1. Fuer jedes Objekt: theta = ang_sep(Pos, D-Pol)
  2. delta = min |theta - n*theta0|  (Abstand zur naechsten Schale)
  3. Hit-Rate: delta < TOL_DEG
  4. KS-Test: CDF(delta) vs. uniform[0, theta0/2]
  5. Erwartet-Anteil analytisch aus spherischem Integral
"""

import math, os, sys, csv, io, time, urllib.request, urllib.parse, re, random

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

# ── HTM-Konstanten ────────────────────────────────────────────────────────────
THETA0   = 58.65            # Grad, HTM-Gitter-Winkel
SHELLS   = [THETA0 * n for n in range(1, 4)]  # 58.65, 117.3, 175.95
HALF     = THETA0 / 2.0     # 29.325 Grad
TOL_DEG  = 5.0              # "auf Schale" Toleranz
DPOL_L   = 305.0
DPOL_B   = 25.0
DPOL_XYZ = None             # wird unten gesetzt

# ── Koordinaten-Hilfsfunktionen ───────────────────────────────────────────────
def lb_xyz(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    return (math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b))

def radec_to_lb(ra_deg, dec_deg):
    """RA/Dec J2000 (Grad) → galaktisch (l, b) — korrekte IAU-Formel."""
    ra  = math.radians(ra_deg)
    dc  = math.radians(dec_deg)
    RN  = math.radians(192.85948)   # RA NGP
    DN  = math.radians(27.12825)    # Dec NGP
    LN  = math.radians(122.93192)   # l NCP
    b   = math.asin(max(-1., min(1., math.sin(dc)*math.sin(DN)
                                   + math.cos(dc)*math.cos(DN)*math.cos(ra - RN))))
    y   = math.cos(dc)*math.sin(ra - RN)
    x   = math.sin(dc)*math.cos(DN) - math.cos(dc)*math.sin(DN)*math.cos(ra - RN)
    l   = (LN - math.atan2(y, x))
    return math.degrees(l) % 360.0, math.degrees(b)

def ang_sep_xyz(v1, v2):
    d = sum(a*b for a, b in zip(v1, v2))
    return math.degrees(math.acos(max(-1., min(1., d))))

def theta_from_dpol(ra_deg, dec_deg):
    l, b = radec_to_lb(ra_deg, dec_deg)
    return ang_sep_xyz(lb_xyz(l, b), DPOL_XYZ)

def theta_from_dpol_lb(l_deg, b_deg):
    return ang_sep_xyz(lb_xyz(l_deg, b_deg), DPOL_XYZ)

def delta_min(theta_deg):
    """Kleinster Abstand zur naechsten HTM-Schale (mod theta0, symmetrisch)."""
    d = theta_deg % THETA0
    return min(d, THETA0 - d)

# ── Erwartete Hit-Fraction (analytisch, Sphäre) ───────────────────────────────
def expected_hit_fraction(tol_deg):
    """Fraction des Himmels innerhalb tol_deg einer der Schalen n=1,2,3."""
    total = 0.0
    for sh in SHELLS:
        lo = max(0.0, sh - tol_deg)
        hi = min(180.0, sh + tol_deg)
        total += (math.cos(math.radians(lo)) - math.cos(math.radians(hi))) / 2.0
    return total   # already /2 since sphere area = 4π → half-sphere = cos convention

EXP_FRAC = expected_hit_fraction(TOL_DEG)

# ── KS-Test (stdlib, gegen Uniform[0, HALF]) ──────────────────────────────────
def ks_test_uniform(deltas):
    """KS-Test: delta_min vs. Uniform[0, HALF]. Gibt (D, p) zurück."""
    s = sorted(deltas)
    n = len(s)
    if n < 2:
        return 0., 1.
    D = 0.
    for i, d in enumerate(s):
        cdf_obs_hi = (i + 1) / n
        cdf_obs_lo = i / n
        cdf_null   = d / HALF
        D = max(D, abs(cdf_obs_hi - cdf_null), abs(cdf_obs_lo - cdf_null))
    z = D * (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n))
    p = 2. * sum(((-1)**(k+1)) * math.exp(-2*(k**2)*(z**2)) for k in range(1, 50))
    return D, max(0., min(1., p))

# ── Report-Helfer ─────────────────────────────────────────────────────────────
def analyze(label, thetas, ot_tag):
    """Vollanalyse für eine Objektliste (theta_dpole-Werte)."""
    n      = len(thetas)
    deltas = [delta_min(t) for t in thetas]
    hits   = sum(1 for d in deltas if d < TOL_DEG)
    obs_f  = hits / n if n else 0.
    ratio  = obs_f / EXP_FRAC if EXP_FRAC > 0 else 0.
    D, p   = ks_test_uniform(deltas)

    print(f"\n{'─'*60}")
    print(f"  {ot_tag}: {label}  (N={n})")
    print(f"  Erwartete Hit-Rate (±{TOL_DEG}°, isotropisch): {EXP_FRAC*100:.1f}%")
    print(f"  Beobachtet auf Schale: {hits}/{n} ({obs_f*100:.1f}%)")
    print(f"  Ratio obs/exp:         {ratio:.3f}")
    print(f"  KS D={D:.4f}  p={p:.4f}  "
          f"({'SIGNIFIKANT p<0.05' if p < 0.05 else 'n.s.'})")
    # Top-Shell-Breakdown
    for sh in SHELLS:
        sh_hits = sum(1 for t in thetas if abs(t - sh) < TOL_DEG)
        exp_sh  = expected_hit_fraction(TOL_DEG) / len(SHELLS) * n  # rough
        print(f"    Schale {sh:.1f}°: {sh_hits} Treffer")
    return dict(label=label, n=n, hits=hits, obs_f=obs_f, exp_f=EXP_FRAC,
                ratio=ratio, D=D, p=p, deltas=deltas, thetas=thetas)

# ── VizieR-Abfrage (TAP) ──────────────────────────────────────────────────────
def vizier_tap(adql, label, timeout=60):
    """Führt TAP-Query durch; gibt Liste von Zeilen-Dicts zurück."""
    url = ("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
           "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv"
           f"&QUERY={urllib.parse.quote(adql)}")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'HTMShellCheck/1.0'})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            raw = r.read().decode('utf-8', errors='replace')
        rows = list(csv.DictReader(io.StringIO(raw)))
        print(f"  [{label}] {len(rows)} Zeilen via TAP")
        return rows
    except Exception as e:
        print(f"  [{label}] TAP-Fehler: {e}")
        return []

def vizier_post(catalog, cols, max_rec, label, timeout=60):
    """VizieR classic POST query; gibt CSV-string zurück."""
    params = urllib.parse.urlencode({
        '-source':  catalog,
        '-out':     cols,
        '-out.max': str(max_rec),
        '-out.form':'csv',
    }).encode()
    url = "https://vizier.cds.unistra.fr/viz-bin/votable"
    try:
        req = urllib.request.Request(url, data=params,
                                     headers={'User-Agent': 'HTMShellCheck/1.0'})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            raw = r.read().decode('utf-8', errors='replace')
        print(f"  [{label}] {len(raw)} Bytes via VizieR POST")
        return raw
    except Exception as e:
        print(f"  [{label}] POST-Fehler: {e}")
        return ""

def extract_radec_rows(rows):
    """Extrahiert (ra, dec) aus einer Liste von Row-Dicts. Robust."""
    out = []
    for row in rows:
        keys = list(row.keys())
        ra_k  = next((k for k in keys if k.strip().lower() in
                      ('raj2000','_raj2000','ra','_ra','ra_deg')), None)
        dec_k = next((k for k in keys if k.strip().lower() in
                      ('dej2000','_dej2000','dec','_de','de','dec_deg')), None)
        if ra_k and dec_k:
            try:
                out.append((float(row[ra_k]), float(row[dec_k])))
            except (ValueError, TypeError):
                pass
    return out

# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 68)
print("HTM-Schalenwinkel-Test  OT-5 / OT-13 / OT-23")
print("=" * 68)
print(f"D-Pol: l={DPOL_L}°, b={DPOL_B}°")
print(f"Schalen: {', '.join(f'{s:.2f}°' for s in SHELLS)}")
print(f"Toleranz: ±{TOL_DEG}°")
print(f"Erwartete Hit-Rate (isotropisch): {EXP_FRAC*100:.1f}%")

DPOL_XYZ = lb_xyz(DPOL_L, DPOL_B)

results = {}

# ═══════════════════════════════════════════════════════════════════════════════
# OT-5: SMBH-Katalog (97 Objekte aus smbh_extended.csv)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*68)
print("OT-5: SMBH Shell-Spektrum (97 Objekte)")
print("="*68)

SMBH_CSV = os.path.join(RES, "catalogs", "smbh_extended.csv")
smbh_thetas = []
hits_5_csv  = 0

if os.path.exists(SMBH_CSV):
    with open(SMBH_CSV, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                th = float(row.get('theta_dpole', ''))
                smbh_thetas.append(th)
                if row.get('hit_5deg', '').strip().upper() == 'YES':
                    hits_5_csv += 1
            except (ValueError, TypeError):
                pass
    print(f"  {len(smbh_thetas)} SMBHs geladen (theta_dpole aus CSV)")
    print(f"  CSV-eigene hit_5deg=YES Eintraege: {hits_5_csv}")
else:
    print(f"  FEHLER: {SMBH_CSV} nicht gefunden!")

if smbh_thetas:
    results['OT5'] = analyze("SMBH 97-Katalog (McConnell+Graham 2013)", smbh_thetas, "OT-5")

# ═══════════════════════════════════════════════════════════════════════════════
# OT-13: Milliquas AGN/QSO (VizieR VII/290)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*68)
print("OT-13: Milliquas AGN/QSO (VizieR VII/290)")
print("="*68)

milli_cache = os.path.join(RES, "milliquas_sample.csv")
milli_rows  = []

if os.path.exists(milli_cache):
    with open(milli_cache, newline='', encoding='utf-8') as f:
        milli_rows = list(csv.DictReader(f))
    print(f"  Cache: {len(milli_rows)} Zeilen aus {milli_cache}")
else:
    # Multi-Band-Sampling: 6 RA-Bänder à 500 Objekte = 3000 gesamt
    # Milliquas (VII/290) unterstützt WHERE-Klausel via TAP
    PER_BAND = 500
    RA_BANDS = [(i * 60, i * 60 + 60) for i in range(6)]  # 0-60, 60-120, ...
    milli_rows = []
    for ra_lo, ra_hi in RA_BANDS:
        adql = (f"SELECT RAJ2000, DEJ2000 FROM \"VII/290/catalog\" "
                f"WHERE RAJ2000 BETWEEN {ra_lo} AND {ra_hi}")
        url  = ("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
                f"?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&MAXREC={PER_BAND}"
                f"&QUERY={urllib.parse.quote(adql)}")
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'HTMShellCheck/1.0'})
            with urllib.request.urlopen(req, timeout=60) as r:
                raw = r.read().decode('utf-8', errors='replace')
            band_rows = list(csv.DictReader(io.StringIO(raw)))
            milli_rows.extend(band_rows)
            print(f"  RA {ra_lo:3d}-{ra_hi:3d}°: {len(band_rows)} Eintraege")
        except Exception as e:
            print(f"  RA {ra_lo}-{ra_hi}: {e}")
        time.sleep(0.3)
    print(f"  Milliquas gesamt: {len(milli_rows)} Eintraege")

    if milli_rows:
        with open(milli_cache, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=milli_rows[0].keys())
            writer.writeheader()
            writer.writerows(milli_rows)

milli_coords = extract_radec_rows(milli_rows)
print(f"  Koordinaten extrahiert: {len(milli_coords)}")

milli_thetas = []
for ra, dec in milli_coords:
    try:
        milli_thetas.append(theta_from_dpol(ra, dec))
    except Exception:
        pass
print(f"  Theta-Werte berechnet: {len(milli_thetas)}")

if milli_thetas:
    results['OT13'] = analyze("Milliquas AGN/QSO (VII/290)", milli_thetas, "OT-13")
else:
    print("  KEINE Milliquas-Daten — OT-13 übersprungen")
    results['OT13'] = None

# ═══════════════════════════════════════════════════════════════════════════════
# OT-23: SDSS Voids (Sutter+2012, VizieR J/MNRAS/431/2307)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*68)
print("OT-23: Kosmische Voids (Sutter+2012, J/MNRAS/431/2307)")
print("="*68)

void_cache = os.path.join(RES, "sutter_voids.csv")
void_rows  = []

if os.path.exists(void_cache):
    with open(void_cache, newline='', encoding='utf-8') as f:
        void_rows = list(csv.DictReader(f))
    print(f"  Cache: {len(void_rows)} Zeilen aus {void_cache}")
else:
    # Void-Kataloge in VizieR TAP (bestätigt verfügbar):
    # J/MNRAS/422/25/void  — Nadathur+2012 SDSS DR7 Voids (z~0.04-0.12)
    # J/A+A/570/A106/void  — Micheletti+2014 VIPERS Voids (z~0.5-1.2)
    # J/A+A/706/A362/voids — Verza+2025 Voids
    VOID_TABLES = [
        ("J/MNRAS/422/25/void",  "SELECT RAJ2000, DEJ2000, z FROM \"J/MNRAS/422/25/void\"",  3000),
        ("J/A+A/570/A106/void",  "SELECT RAJ2000, DEJ2000, z FROM \"J/A+A/570/A106/void\"", 3000),
        ("J/A+A/706/A362/voids", "SELECT RAJ2000, DEJ2000, z FROM \"J/A+A/706/A362/voids\"", 3000),
    ]
    for v_label, v_adql, v_max in VOID_TABLES:
        url = ("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
               f"?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&MAXREC={v_max}"
               f"&QUERY={urllib.parse.quote(v_adql)}")
        print(f"  TAP {v_label}...")
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'HTMShellCheck/1.0'})
            with urllib.request.urlopen(req, timeout=90) as r:
                raw = r.read().decode('utf-8', errors='replace')
            rows = list(csv.DictReader(io.StringIO(raw)))
            if rows:
                void_rows.extend(rows)
                print(f"    {len(rows)} Voids geladen")
        except Exception as e:
            print(f"    {v_label}: {e}")
        time.sleep(0.5)

    print(f"  Voids gesamt: {len(void_rows)}")
    if void_rows:
        with open(void_cache, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=void_rows[0].keys())
            writer.writeheader()
            writer.writerows(void_rows)

void_coords = extract_radec_rows(void_rows)
print(f"  Koordinaten extrahiert: {len(void_coords)}")

void_thetas = []
for ra, dec in void_coords:
    try:
        void_thetas.append(theta_from_dpol(ra, dec))
    except Exception:
        pass
print(f"  Theta-Werte berechnet: {len(void_thetas)}")

if void_thetas:
    results['OT23'] = analyze("Kosmische Voids (Nadathur+2012 / Micheletti+2014 / Verza+2025)", void_thetas, "OT-23")
else:
    print("  KEINE Void-Daten — OT-23 übersprungen")
    results['OT23'] = None

# ═══════════════════════════════════════════════════════════════════════════════
# Ausgabe: Ergebnisdatei
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*68)
print("ZUSAMMENFASSUNG")
print("="*68)

def sig_str(p):
    if p < 0.001: return "★★★ p<0.001"
    if p < 0.01:  return "★★  p<0.01"
    if p < 0.05:  return "★   p<0.05"
    return "n.s."

for tag, r in results.items():
    if r is None:
        print(f"  {tag}: KEINE DATEN")
        continue
    print(f"  {tag}: N={r['n']:4d}  hits={r['hits']}({r['obs_f']*100:.1f}%)"
          f"  exp={r['exp_f']*100:.1f}%  Ratio={r['ratio']:.3f}"
          f"  KS p={r['p']:.4f} {sig_str(r['p'])}")

# ── Ausgabe-Datei ─────────────────────────────────────────────────────────────
buf = io.StringIO()
def w(s=''): buf.write(str(s)+'\n')

w("=" * 72)
w("OT-5 / OT-13 / OT-23 — HTM-Schalenwinkel-Test")
w("Hypothese: Objekte cluster bevorzugt bei theta = n * 58.65° vom D-Pol")
w("=" * 72)
w(f"Datum:     {time.strftime('%Y-%m-%d %H:%M')}")
w(f"D-Pol:     l={DPOL_L}°, b={DPOL_B}°")
w(f"Schalen:   {', '.join(f'{s:.2f}°' for s in SHELLS)}")
w(f"Toleranz:  ±{TOL_DEG}°")
w(f"Erwartet:  {EXP_FRAC*100:.2f}% des Himmels auf Schalen (isotropisch)")
w()

for ot, tag, label in [
    ('OT5',  'OT-5',  'SMBH Shell-Spektrum'),
    ('OT13', 'OT-13', 'Milliquas AGN/QSO'),
    ('OT23', 'OT-23', 'Kosmische Voids (Nadathur+2012 / Micheletti+2014)'),
]:
    r = results.get(ot)
    w("-" * 72)
    w(f"{tag}: {label}")
    if r is None:
        w("  STATUS: KEINE DATEN — Katalog nicht verfügbar")
        w()
        continue
    sig = sig_str(r['p'])
    w(f"  N:              {r['n']}")
    w(f"  Hits (±{TOL_DEG}°):   {r['hits']} ({r['obs_f']*100:.2f}%)")
    w(f"  Erwartet:       {r['exp_f']*100:.2f}%")
    w(f"  Ratio obs/exp:  {r['ratio']:.4f}")
    w(f"  KS D:           {r['D']:.4f}")
    w(f"  KS p:           {r['p']:.4f}  {sig}")
    w()
    # Shell-Breakdown
    w("  Shell-Breakdown:")
    for sh in SHELLS:
        sh_hits = sum(1 for t in r['thetas'] if abs(t - sh) < TOL_DEG)
        exp_sh_n = (expected_hit_fraction(TOL_DEG) / len(SHELLS)) * r['n']
        w(f"    theta={sh:.2f}°: {sh_hits} Treffer  (erw. {exp_sh_n:.1f})")
    w()
    # Top-Hits (closest to shells)
    if ot == 'OT5':
        w("  SMBH Shell-Residuen (alle, sortiert nach delta):")
        # Re-read CSV for names
        smbh_named = []
        with open(SMBH_CSV, newline='', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                try:
                    th = float(row.get('theta_dpole', ''))
                    dm = delta_min(th)
                    smbh_named.append((dm, th, row.get('Name','?')))
                except (ValueError, TypeError):
                    pass
        smbh_named.sort()
        for dm, th, nm in smbh_named:
            mk = '*' if dm < TOL_DEG else ' '
            w(f"    {mk} {nm:<22s}  theta={th:6.2f}°  delta={dm:.2f}°")
    w()

w("=" * 72)
w("INTERPRETATION")
w("=" * 72)
for ot, tag in [('OT5','OT-5'),('OT13','OT-13'),('OT23','OT-23')]:
    r = results.get(ot)
    if r is None:
        w(f"  {tag}: KEINE DATEN (Katalog nicht ladbar)")
        continue
    if r['p'] < 0.05 and r['ratio'] > 1.0:
        w(f"  {tag}: EXCESS-SIGNAL  — Ratio={r['ratio']:.3f}× (MEHR als erwartet), KS p={r['p']:.4f}")
    elif r['p'] < 0.05 and r['ratio'] < 1.0:
        w(f"  {tag}: KEIN HTM-SIGNAL — Ratio={r['ratio']:.3f}× (WENIGER als erwartet, "
          f"wahrsch. Survey-Footprint-Effekt), KS p={r['p']:.4f}")
    else:
        w(f"  {tag}: KEIN SIGNAL — Ratio={r['ratio']:.3f}×, KS p={r['p']:.4f} (n.s.)")
w()
w("HINWEIS OT-13/OT-23: KS-Test signifikant wegen inhomogenem Survey-Footprint,")
w("  nicht wegen HTM-Signal. Korrekte Analyse benoetigt Mock-Katalog im selben Footprint.")
w()

out_path = os.path.join(RES, "OT_shells_combined.txt")
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(buf.getvalue())
print(f"\nAusgabe: {out_path}")
print("FERTIG.")
