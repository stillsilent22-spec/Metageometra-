"""
Multi-Katalog Thales-Analyse — Vakuumschalen-Radius
====================================================
Kombiniert alle verfuegbaren Galaxien-Distanz-Kataloge:
  CF3  — CosmicFlows-3  (Tully+2016, VizieR J/AJ/152/50)    ~18k, bis ~500 Mpc
  CF2  — CosmicFlows-2  (Tully+2013, VizieR J/AJ/146/86)     ~8k, bis ~300 Mpc
  2MRS — 2MASS Redshift (Huchra+2012, VizieR J/ApJS/199/26) ~43k, cz/H0
  SMBH — smbh_extended.csv                                   ~61 mit Distanz

Formel: R_Thales = d_Mpc / cos(theta_D)
"""
import math, os, sys, csv, urllib.request, urllib.parse, time
import xml.etree.ElementTree as ET
import statistics as st

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE  = os.path.dirname(os.path.abspath(__file__))
RES   = os.path.join(BASE, "results")
OUT_TXT = os.path.join(RES, "vizier_thales.txt")
OUT_CSV = os.path.join(RES, "catalogs", "vizier_dpol_combined.csv")
os.makedirs(os.path.join(RES, "catalogs"), exist_ok=True)

# ── HTM-Konstanten ──────────────────────────────────────────────────────────
A0  = 1.0964e-10          # m/s^2
C   = 2.998e8             # m/s
MPC = 3.0857e22           # m
r_s = (C**2) / (2*math.pi*A0) / MPC   # ≈ 4228 Mpc
H0  = 73.0                # km/s/Mpc fuer cz->Mpc Konversion
THETA0 = 58.65            # HTM Fraktal-Winkel

DPOLE_L, DPOLE_B = 305.0, 25.0
KEGEL_DEG        = 45.0   # Oeffnungswinkel um D-Pol

# ── Koordinaten-Hilfsfunktionen ──────────────────────────────────────────────
def galactic_to_equatorial(l_deg, b_deg):
    """Galaktisch (J2000-Achse) -> Aequatorial J2000"""
    l = math.radians(l_deg);  b  = math.radians(b_deg)
    RA_NGP  = math.radians(192.85948)
    DEC_NGP = math.radians(27.12825)
    L_NCP   = math.radians(122.93192)
    sin_d   = (math.sin(b)*math.sin(DEC_NGP) +
               math.cos(b)*math.cos(DEC_NGP)*math.cos(L_NCP - l))
    dec     = math.asin(max(-1.0, min(1.0, sin_d)))
    y  = math.cos(b)*math.sin(L_NCP - l)
    x  = (math.sin(b)*math.cos(DEC_NGP) -
          math.cos(b)*math.sin(DEC_NGP)*math.cos(L_NCP - l))
    ra = (RA_NGP + math.atan2(y, x))
    return math.degrees(ra) % 360, math.degrees(dec)

def radec_to_lb(ra_deg, dec_deg):
    ra  = math.radians(ra_deg);  dc  = math.radians(dec_deg)
    RA_NGP  = math.radians(192.85948)
    DEC_NGP = math.radians(27.12825)
    L_NCP   = math.radians(122.93192)
    b   = math.asin(max(-1., min(1., math.sin(dc)*math.sin(DEC_NGP) +
                math.cos(dc)*math.cos(DEC_NGP)*math.cos(ra - RA_NGP))))
    y   = math.cos(dc)*math.sin(ra - RA_NGP)
    x   = (math.sin(dc)*math.cos(DEC_NGP) -
           math.cos(dc)*math.sin(DEC_NGP)*math.cos(ra - RA_NGP))
    l   = (L_NCP - math.atan2(y, x))
    return math.degrees(l) % 360, math.degrees(b)

def ang_sep(l1, b1, l2, b2):
    l1r,b1r = math.radians(l1), math.radians(b1)
    l2r,b2r = math.radians(l2), math.radians(b2)
    c = (math.sin(b1r)*math.sin(b2r) +
         math.cos(b1r)*math.cos(b2r)*math.cos(l1r - l2r))
    return math.degrees(math.acos(max(-1., min(1., c))))

# ── VizieR-Query ─────────────────────────────────────────────────────────────
dpol_ra, dpol_dec = galactic_to_equatorial(DPOLE_L, DPOLE_B)
print(f"D-Pol: l={DPOLE_L}°, b={DPOLE_B}° → RA={dpol_ra:.2f}°, Dec={dpol_dec:.2f}°")
print(f"HTM Vakuumschalen-Radius r_s = {r_s:.1f} Mpc")

# Hilfsfunktion: VizieR abfragen (generisch) – multi-table VOTable aware, mit Retry
def vizier_query(source, cols, label, max_out=30000):
    import re as _re
    center = f"{dpol_ra:.4f} {dpol_dec:+.4f}"  # Space required for negative Dec
    params = {'-source': source, '-c': center, '-c.r': str(KEGEL_DEG),
              '-c.u': 'deg', '-out': cols, '-out.max': str(max_out), '-oc.form': 'dec'}
    url = 'https://vizier.cds.unistra.fr/viz-bin/votable?' + urllib.parse.urlencode(params)
    print(f"  Query {label} [{source}]...", end=' ', flush=True)
    for attempt in range(3):
        try:
            if attempt > 0:
                time.sleep(5 * attempt)
                print(f"  Retry {attempt}...", end=' ', flush=True)
            req = urllib.request.Request(url, headers={'User-Agent': 'Python/HTMResearch'})
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read().decode('utf-8', 'replace')
            # Handle multi-table VOTables: find each TABLE section independently
            rows_all = []
            for table_match in _re.finditer(r'<TABLE\b[^>]*>(.+?)</TABLE>', raw, _re.DOTALL):
                table_body = table_match.group(1)
                fields = _re.findall(r'<FIELD[^>]+name="([^"]+)"', table_body)
                fields = [f for f in fields if f not in ('recno', '_ivoid')]
                if not fields:
                    continue
                for tr_m in _re.finditer(r'<TR>(.*?)</TR>', table_body, _re.DOTALL):
                    tds = _re.findall(r'<TD>(.*?)</TD>', tr_m.group(1))
                    if len(tds) < 2:
                        continue
                    row = {fields[i]: tds[i].strip() if i < len(tds) else '' for i in range(len(fields))}
                    rows_all.append(row)
            if rows_all or attempt == 2:
                print(f"{len(rows_all)} Zeilen")
                return rows_all
            # 0 rows and not last attempt: retry
            print(f"0 (Retry)...", end=' ', flush=True)
        except Exception as e:
            print(f"FEHLER: {e}")
            if attempt == 2:
                return []
    return []

# Koordinaten-Parser (Sexagesimal oder Dezimal)
def parse_coord(s, is_ra=True):
    s = str(s).strip()
    if ':' in s or (s.count(' ') >= 2):
        parts = s.replace(':',' ').split()
        neg = s.startswith('-')
        val = abs(float(parts[0])) + (float(parts[1]) if len(parts)>1 else 0)/60 + (float(parts[2]) if len(parts)>2 else 0)/3600
        if neg: val = -val
        if is_ra: val *= 15
        return val
    return float(s)

# CF3: aus lokalem Cache laden (bereits abgerufen)
center_cf3 = f"{dpol_ra:.4f} {dpol_dec:+.4f}"
params = {
    '-source': 'J/AJ/152/50',
    '-c':      center_cf3,
    '-c.r':    str(KEGEL_DEG),
    '-c.u':    'deg',
    '-out':    'Name,RAJ2000,DEJ2000,Dist',
    '-out.max':'10000',
    '-oc.form':'dec',
}
url = 'https://vizier.cds.unistra.fr/viz-bin/votable?' + urllib.parse.urlencode(params)
print(f"\nQuery: VizieR CosmicFlows-3 in {KEGEL_DEG}°-Kegel um D-Pol...")

try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/HTMResearch'})
    with urllib.request.urlopen(req, timeout=90) as r:
        raw = r.read().decode('utf-8', 'replace')
    print(f"  Antwort: {len(raw):,} Zeichen")
except Exception as e:
    print(f"FEHLER CF3: {e}")
    raw = ''

# ── VOTable parsen (Regex-Methode fuer VizieR) ───────────────────────────────
import re

def parse_vizier_votable(xml_text):
    """Extrahiert FIELD-Namen und TR/TD-Daten aus VizieR-VOTable via Regex."""
    # FIELD-Namen (letzte TABLE mit FIELDs benutzen)
    field_names = re.findall(r'<FIELD[^>]+name="([^"]+)"', xml_text)
    # Leere Hidden-Felder am Ende ignorieren (recno, _ivoid)
    field_names = [f for f in field_names if f not in ('recno',)]

    # Alle TR-Zeilen mit Inhalt
    rows = []
    for tr_match in re.finditer(r'<TR>(.*?)</TR>', xml_text, re.DOTALL):
        tds = re.findall(r'<TD>(.*?)</TD>', tr_match.group(1))
        if len(tds) < 3:  # mindestens Name + RA + Dec
            continue
        row = {}
        for i, fn in enumerate(field_names):
            row[fn] = tds[i].strip() if i < len(tds) else ''
        rows.append(row)
    return rows

rows = parse_vizier_votable(raw)
print(f"  Eintraege geparsed: {len(rows)}")
if rows:
    print(f"  Spalten: {list(rows[0].keys())}")
    print(f"  Beispiel: {rows[0]}")

# ── Thales-Berechnung ─────────────────────────────────────────────────────────
galaxies = []
for row in rows:
    try:
        # CosmicFlows-3 gibt RA/Dec als Sexagesimal oder Dezimal
        name = row.get('Name','?').strip()
        ra_s = row.get('RAJ2000','').strip()
        dec_s= row.get('DEJ2000','').strip()
        dist_s=row.get('Dist','').strip()

        if not dist_s or dist_s in ('','NaN','-','-  '):
            continue

        # RA: "HH MM SS.ss" oder Dezimalgrad
        def parse_angle(s, is_ra=True):
            s = s.strip()
            if ':' in s or (s.count(' ') >= 2):
                parts = s.replace(':',' ').split()
                d = abs(float(parts[0]))
                m = float(parts[1]) if len(parts)>1 else 0
                sc= float(parts[2]) if len(parts)>2 else 0
                val = d + m/60 + sc/3600
                if s.startswith('-'): val = -val
                if is_ra: val *= 15  # Stunden → Grad
                return val
            return float(s)

        ra_deg  = parse_angle(ra_s,  is_ra=True)
        dec_deg = parse_angle(dec_s, is_ra=False)
        d_mpc   = float(dist_s)
        if d_mpc <= 0:
            continue

        l, b    = radec_to_lb(ra_deg, dec_deg)
        theta   = ang_sep(l, b, DPOLE_L, DPOLE_B)

        if theta < 0.5 or theta > KEGEL_DEG:
            continue
        cos_t = math.cos(math.radians(theta))
        if cos_t < 0.05:
            continue

        r_thales = d_mpc / cos_t
        galaxies.append({
            'name': name, 'ra': ra_deg, 'dec': dec_deg,
            'l': l, 'b': b, 'theta_D': theta,
            'd_mpc': d_mpc, 'R_Thales': r_thales, 'source': 'CF3'
        })
    except (ValueError, ZeroDivisionError):
        continue

print(f"  CF3: {len(galaxies)} Galaxien im D-Pol-Kegel")

# ── Deduplizierung ─────────────────────────────────────────────────────────────
seen_keys = set((round(g['l'],1), round(g['b'],1), round(g['d_mpc'],0)) for g in galaxies)

def add_galaxy(g):
    key = (round(g['l'],1), round(g['b'],1), round(g['d_mpc'],0))
    if key not in seen_keys:
        seen_keys.add(key)
        galaxies.append(g)
        return True
    return False

# ── KATALOG 2: CosmicFlows-2 (J/AJ/146/86) ────────────────────────────────────
rows_cf2 = vizier_query('J/AJ/146/86', 'LEDA,RAJ2000,DEJ2000,Dist', 'CF2')
cf2_new = 0
for row in rows_cf2:
    try:
        d = float(row.get('Dist','0'))
        if d <= 0: continue
        name = row.get('LEDA', row.get('Name','CF2?')).strip()
        ra  = parse_coord(row.get('RAJ2000',''), is_ra=True)
        dec = parse_coord(row.get('DEJ2000',''), is_ra=False)
        l, b = radec_to_lb(ra, dec)
        theta = ang_sep(l, b, DPOLE_L, DPOLE_B)
        if theta < 0.5 or theta > KEGEL_DEG: continue
        cos_t = math.cos(math.radians(theta))
        if cos_t < 0.05: continue
        if add_galaxy({'name':name,'ra':ra,'dec':dec,'l':l,'b':b,'theta_D':theta,
                       'd_mpc':d,'R_Thales':d/cos_t,'source':'CF2'}):
            cf2_new += 1
    except (ValueError, TypeError): continue
print(f"  CF2: {cf2_new} neue Objekte hinzugefuegt")
time.sleep(3)  # VizieR rate-limit

# ── KATALOG 3: 2MRS (J/ApJS/199/26) — cz/H0 ──────────────────────────────────
rows_2mrs = vizier_query('J/ApJS/199/26', 'RAJ2000,DEJ2000,cz', '2MRS', max_out=50000)
mrs_new = 0
for row in rows_2mrs:
    try:
        ra  = parse_coord(row.get('RAJ2000','0'), is_ra=True)
        dec = parse_coord(row.get('DEJ2000','0'), is_ra=False)
        cz  = float(row.get('cz','0') or '0')
        if cz < 500: continue   # Eigenbewegung dominiert
        d = cz / H0
        l, b = radec_to_lb(ra, dec)
        theta = ang_sep(l, b, DPOLE_L, DPOLE_B)
        if theta < 0.5 or theta > KEGEL_DEG: continue
        cos_t = math.cos(math.radians(theta))
        if cos_t < 0.05: continue
        name = f"2MRS_{ra:.2f}{dec:+.2f}"
        if add_galaxy({'name':name,'ra':ra,'dec':dec,'l':l,'b':b,'theta_D':theta,
                       'd_mpc':d,'R_Thales':d/cos_t,'source':'2MRS'}):
            mrs_new += 1
    except (ValueError, TypeError): continue
print(f"  2MRS: {mrs_new} neue Objekte (cz/H0, cz>500 km/s)")

# ── KATALOG 4: SMBH-Katalog ────────────────────────────────────────────────────
smbh_csv_path = os.path.join(os.path.join(RES, 'catalogs'), 'smbh_extended.csv')
smbh_new = 0
if os.path.exists(smbh_csv_path):
    with open(smbh_csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                d = float(row.get('distance_mpc', row.get('dist','0')) or '0')
                if d <= 0: continue
                ra  = float(row.get('ra', 0))
                dec = float(row.get('dec', 0))
                l, b = radec_to_lb(ra, dec)
                theta = ang_sep(l, b, DPOLE_L, DPOLE_B)
                if theta < 0.5 or theta > KEGEL_DEG: continue
                cos_t = math.cos(math.radians(theta))
                if cos_t < 0.05: continue
                if add_galaxy({'name':row.get('name','SMBH?'),'ra':ra,'dec':dec,'l':l,'b':b,'theta_D':theta,
                               'd_mpc':d,'R_Thales':d/cos_t,'source':'SMBH'}):
                    smbh_new += 1
            except (ValueError, TypeError): continue
    print(f"  SMBH: {smbh_new} neue Objekte")

galaxies.sort(key=lambda x: x['d_mpc'])
print(f"\nGESAMT: {len(galaxies)} Objekte im D-Pol-Kegel ({KEGEL_DEG}°)")
print(f"  Quellen: CF3={len(galaxies)-cf2_new-mrs_new-smbh_new}, CF2={cf2_new}, 2MRS={mrs_new}, SMBH={smbh_new}")
print()

# ── Ausgabe ───────────────────────────────────────────────────────────────────
lines = []
src_counts = {s: sum(1 for g in galaxies if g.get('source')==s) for s in ['CF3','CF2','2MRS','SMBH']}
lines.append("="*72)
lines.append("MULTI-KATALOG THALES-ANALYSE: VAKUUMSCHALEN-RADIUS")
lines.append(f"D-Pol: l={DPOLE_L}°, b={DPOLE_B}°  |  Kegel: {KEGEL_DEG}°")
lines.append(f"Formel: R_Thales = d_Mpc / cos(theta_D)  |  H0 = {H0} km/s/Mpc")
lines.append(f"HTM Vakuumschalen-Radius r_s = {r_s:.1f} Mpc")
lines.append(f"Kataloge: CF3={src_counts['CF3']}, CF2={src_counts['CF2']}, 2MRS={src_counts['2MRS']}, SMBH={src_counts['SMBH']}")
lines.append("="*72)
lines.append(f"\n{'Name':<24}{'Src':>5}{'d[Mpc]':>9}{'theta':>8}{'R_Thales':>12}")
lines.append("-"*62)

for g in galaxies[:100]:
    lines.append(f"  {g['name']:<22}{g.get('source','?'):>5}{g['d_mpc']:>9.1f}{g['theta_D']:>7.1f}°{g['R_Thales']:>12.1f}")
if len(galaxies) > 100:
    lines.append(f"  ... ({len(galaxies)-100} weitere Objekte)")

if galaxies:
    rv = [g['R_Thales'] for g in galaxies]
    dv = [g['d_mpc']    for g in galaxies]
    med_r = st.median(rv)
    avg_r = st.mean(rv)
    std_r = st.stdev(rv) if len(rv) > 1 else 0

    clust = [r for r in rv if 0.5*med_r < r < 1.5*med_r]

    lines.append("\n" + "="*55)
    lines.append(f"STATISTIK ({len(rv)} Galaxien):")
    lines.append(f"  Median R_Thales   : {med_r:>10.1f} Mpc")
    lines.append(f"  Mittelwert        : {avg_r:>10.1f} Mpc")
    lines.append(f"  Std               : {std_r:>10.1f} Mpc")
    lines.append(f"  Std/Median        : {std_r/med_r:.3f}")
    lines.append(f"  Cluster ±50%      : {len(clust)}/{len(rv)} = {100*len(clust)/len(rv):.0f}%")
    lines.append(f"\n  Distanz-Bereich   : {min(dv):.1f} – {max(dv):.1f} Mpc")
    lines.append(f"  HTM r_s           : {r_s:.1f} Mpc")
    lines.append(f"  Median/r_s        : {med_r/r_s:.4f}  ({100*med_r/r_s:.2f}%)")

    # Distanz-Dekaden
    lines.append("\nNACH DISTANZ-DEKADE (alle Kataloge kombiniert):")
    lines.append(f"  {'Dekade':>12} {'N':>6} {'Median-R':>10} {'Std':>10} {'Cl25%':>8}")
    lines.append("  " + "-"*50)
    dekaden = [(0,5),(5,15),(15,50),(50,100),(100,200),(200,500),(500,1000)]
    for d_lo, d_hi in dekaden:
        g2 = [g for g in galaxies if d_lo < g['d_mpc'] <= d_hi]
        if g2:
            rv2 = [g['R_Thales'] for g in g2]
            med2 = st.median(rv2)
            std2 = st.stdev(rv2) if len(rv2)>1 else 0
            cl25 = sum(1 for r in rv2 if 0.75*med2 < r < 1.25*med2)
            lines.append(f"  d={d_lo:>4}–{d_hi:<5} Mpc  {len(rv2):>6}  {med2:>10.1f}  {std2:>10.1f}  {cl25}/{len(rv2)}")

    # HTM Fraktalschalen-Kandidaten
    lines.append(f"\nHTM FRAKTALSCHALEN (theta0={THETA0}°, r_s={r_s:.0f} Mpc):")
    for n in range(1, 8):
        r_n = r_s * n * THETA0 / 360
        cands = [g for g in galaxies if 0.8*r_n < g['R_Thales'] < 1.2*r_n]
        tag = "  <-- KANDIDAT!" if len(cands) > 30 else ""
        lines.append(f"  n={n}  R_n={r_n:8.0f} Mpc  N={len(cands):5d} Objekte{tag}")

    # Interpretation
    lines.append("\nINTERPRETATION:")
    ratio = std_r / med_r
    if ratio < 0.3:
        lines.append(f"  KOHAERENTE SCHALE gefunden bei R ~ {med_r:.0f} Mpc (Std/Med={ratio:.3f})!")
    elif ratio < 0.7:
        lines.append(f"  Moderate Streuung — moegliche Struktur bei R ~ {med_r:.0f} Mpc (Std/Med={ratio:.3f})")
    else:
        lines.append(f"  Hohe Streuung (Std/Med={ratio:.3f}) — R_Thales skaliert mit d.")
        lines.append(f"  Lokale Kataloge reichen nicht fuer Vakuumschalen-Test.")
    lines.append(f"  Katalog-Limit: {max(dv):.0f} Mpc  |  Benoetigt fuer r_s: z ~ {r_s*H0/3e5:.2f}")
    lines.append(f"  Naechster testbarer Schalen-Kandidat: r_s * ({THETA0}/360) = {r_s*THETA0/360:.0f} Mpc")

    # Interpretation
    lines.append("\nINTERPRETATION:")
    lines.append(f"  Satz des Thales: Wenn alle Galaxien auf einer Schale liegen,")
txt = '\n'.join(lines)
print(txt)

with open(OUT_TXT, 'w', encoding='utf-8') as f:
    f.write(txt)
print(f"\nGespeichert: {OUT_TXT}")

if galaxies:
    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['name','source','ra','dec','l','b','theta_D','d_mpc','R_Thales'],
                          extrasaction='ignore')
        w.writeheader()
        w.writerows(galaxies)
    print(f"CSV: {OUT_CSV} ({len(galaxies)} Objekte)")
