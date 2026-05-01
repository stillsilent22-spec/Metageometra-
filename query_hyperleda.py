"""
VizieR-Query: Galaxien in D-Pol-Richtung fuer Thales-Test
==========================================================
Fragt VizieR (CDS Strassburg) nach CosmicFlows-3-Galaxien mit Distanzen
in einem Kegel um den D-Pol (l=305, b=25).
Berechnet dann R_Thales = d / cos(theta_D) fuer alle Treffer.

Quelle: CosmicFlows-3, Tully et al. 2016, VizieR J/AJ/152/50
"""
import math, os, sys, csv, urllib.request, urllib.parse, xml.etree.ElementTree as ET

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
OUT_TXT  = os.path.join(RESULTS, "hyperleda_thales.txt")
OUT_CSV  = os.path.join(RESULTS, "catalogs", "hyperleda_dpol.csv")

# D-Pol in galaktischen Koordinaten
DPOLE_L, DPOLE_B = 305.0, 25.0

# HTM-Konstante
A0   = 1.0964e-10          # m/s^2
C    = 2.998e8             # m/s
MPC  = 3.0857e22           # m
r_s  = (C**2) / (2*math.pi*A0) / MPC  # Mpc

# ---------------------------------------------------------------
def galactic_to_equatorial(l_deg, b_deg):
    """Galaktisch -> Aequatorial (J2000), benoetigt fuer HyperLeda-Kegel"""
    l = math.radians(l_deg); b = math.radians(b_deg)
    RA_NGP  = math.radians(192.85948)
    DEC_NGP = math.radians(27.12825)
    L_NCP   = math.radians(122.93192)
    sin_d = math.sin(b)*math.sin(DEC_NGP) + math.cos(b)*math.cos(DEC_NGP)*math.cos(L_NCP - l)
    dec   = math.asin(sin_d)
    cos_ra_num = math.cos(b)*math.sin(L_NCP - l)
    cos_ra_den = math.cos(dec)
    ra = RA_NGP + math.atan2(cos_ra_num, cos_ra_den * (math.sin(DEC_NGP)*math.cos(L_NCP-l)*math.cos(b) - math.sin(b)*math.cos(DEC_NGP)) / cos_ra_den)
    # simpler direct formula
    sin_ra_m_RA = math.cos(b)*math.sin(L_NCP - l)
    cos_ra_m_RA = math.sin(b)*math.cos(DEC_NGP) - math.cos(b)*math.sin(DEC_NGP)*math.cos(L_NCP - l)
    ra  = RA_NGP + math.atan2(sin_ra_m_RA, cos_ra_m_RA)
    ra  = math.degrees(ra) % 360
    dec = math.degrees(dec)
    return ra, dec

def radec_to_lb(ra_deg, dec_deg):
    ra_r = math.radians(ra_deg); dc_r = math.radians(dec_deg)
    RA_NGP  = math.radians(192.85948)
    DEC_NGP = math.radians(27.12825)
    L_NCP   = math.radians(122.93192)
    b_r = math.asin(math.sin(dc_r)*math.sin(DEC_NGP) +
                    math.cos(dc_r)*math.cos(DEC_NGP)*math.cos(ra_r - RA_NGP))
    y = math.cos(dc_r)*math.sin(ra_r - RA_NGP)
    x = math.sin(dc_r)*math.cos(DEC_NGP) - math.cos(dc_r)*math.sin(DEC_NGP)*math.cos(ra_r - RA_NGP)
    l_r = L_NCP - math.atan2(y, x)
    return math.degrees(l_r) % 360, math.degrees(b_r)

def angular_sep_deg(l1, b1, l2, b2):
    l1r, b1r = math.radians(l1), math.radians(b1)
    l2r, b2r = math.radians(l2), math.radians(b2)
    cos_theta = (math.sin(b1r)*math.sin(b2r) +
                 math.cos(b1r)*math.cos(b2r)*math.cos(l1r - l2r))
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))

def modz_to_mpc(mod0):
    """Distanzmodulus -> Mpc"""
    return 10**((mod0 - 25) / 5.0)

# ---------------------------------------------------------------
# HyperLeda SQL-Query via Web-API
# Kegel: 45 Grad um D-Pol + 45 Grad Antipol
# Felder: pgc, objname, al2000 (RA), de2000 (Dec), mod0, mabs, logd25
# mod0 = korrigierter Distanzmodulus, mabs = abs. Helligkeit
# ---------------------------------------------------------------

# Konvertiere D-Pol zu Aequatorial fuer die Query
dpol_ra, dpol_dec = galactic_to_equatorial(DPOLE_L, DPOLE_B)
print(f"D-Pol: l={DPOLE_L}, b={DPOLE_B} -> RA={dpol_ra:.2f}, Dec={dpol_dec:.2f}")

# HyperLeda SQL: alle Galaxien mit bekanntem mod0 (Distanzmodulus)
# Wir holen einen breiten Datensatz und filtern lokal
# Limit auf 10000 Eintraege, sortiert nach Helligkeit (viele davon haben mod0)
HYPERLEDA_URL = "http://leda.univ-lyon1.fr/leda/fullsql.html"

# Query: Objekte mit mod0 bekannt, v < 30000 km/s (~z<0.1, d < ~430 Mpc)
# al2000 = RA in Stunden (0-24), de2000 = Dec in Grad
sql_query = """SELECT objname, al2000, de2000, mod0, vgsr, vmaxg, type 
FROM meandata 
WHERE mod0 IS NOT NULL AND mod0 < 38.0
ORDER BY mod0 ASC"""

params = urllib.parse.urlencode({
    'sql': sql_query,
    'nra': 'l',
    'd': 'N'
})

url = f"{HYPERLEDA_URL}?{params}"
print(f"\nSende Anfrage an HyperLeda...")
print(f"URL: {url[:120]}...")

try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/HTM-Research'})
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read().decode('utf-8', errors='replace')
    print(f"Antwort erhalten: {len(raw)} Zeichen")
except Exception as e:
    print(f"HyperLeda-Fehler: {e}")
    print("\nVersuch mit alternativer URL (VizieR/Simbad)...")
    raw = None

# ---------------------------------------------------------------
# Parse HyperLeda HTML-Tabelle
# ---------------------------------------------------------------
galaxies = []

if raw:
    # HyperLeda gibt HTML + eingebettete Datentabelle zurueck
    # Datenzeilen sind durch <pre> oder direkt als Text
    import re
    # Suche nach Datenblock
    pre_match = re.search(r'<pre[^>]*>(.*?)</pre>', raw, re.DOTALL | re.IGNORECASE)
    if pre_match:
        data_text = pre_match.group(1)
    else:
        # Manchmal direkt als Text nach Header
        data_text = raw

    lines = data_text.split('\n')
    parsed = 0
    header_found = False
    col_objname = col_ra = col_dec = col_mod0 = col_v = -1
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        # Header-Zeile erkennen
        if 'objname' in line.lower() and ('al2000' in line.lower() or 'ra' in line.lower()):
            parts = line.split('|')
            if len(parts) < 3:
                parts = line.split()
            for i, p in enumerate(parts):
                p = p.strip().lower()
                if 'objname' in p: col_objname = i
                elif 'al2000' in p or p == 'ra': col_ra = i
                elif 'de2000' in p or p == 'dec': col_dec = i
                elif 'mod0' in p: col_mod0 = i
                elif 'vgsr' in p or p == 'v': col_v = i
            header_found = True
            print(f"Header gefunden: objname={col_objname}, ra={col_ra}, dec={col_dec}, mod0={col_mod0}")
            continue
        
        if not header_found:
            continue
        if line.startswith('---') or line.startswith('==='):
            continue
            
        # Datenzeile
        parts = line.split('|')
        if len(parts) < 4:
            parts = line.split()
        if len(parts) < 4:
            continue
        
        try:
            objname = parts[col_objname].strip() if col_objname >= 0 else parts[0].strip()
            ra_h    = float(parts[col_ra].strip())  if col_ra  >= 0 else float(parts[1])
            dec_d   = float(parts[col_dec].strip()) if col_dec >= 0 else float(parts[2])
            mod0    = float(parts[col_mod0].strip()) if col_mod0 >= 0 else float(parts[3])
        except (ValueError, IndexError):
            continue
        
        ra_deg = ra_h * 15.0  # Stunden -> Grad
        l, b   = radec_to_lb(ra_deg, dec_d)
        theta  = angular_sep_deg(l, b, DPOLE_L, DPOLE_B)
        d_mpc  = modz_to_mpc(mod0)
        
        galaxies.append({
            'name': objname, 'ra': ra_deg, 'dec': dec_d,
            'l': l, 'b': b, 'theta_D': theta,
            'd_mpc': d_mpc, 'mod0': mod0
        })
        parsed += 1

    print(f"Geparst: {parsed} Galaxien")

# Fallback: falls HyperLeda nicht erreichbar, lokale smbh_extended verwenden
if not galaxies:
    print("\nFallback: Verwende smbh_extended.csv...")
    csv_path = os.path.join(RESULTS, "catalogs", "smbh_extended.csv")
    if os.path.exists(csv_path):
        with open(csv_path, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    d = float(row.get('distance_mpc', row.get('dist', '')))
                    ra = float(row.get('ra', 0))
                    dec = float(row.get('dec', 0))
                    l, b = radec_to_lb(ra, dec)
                    theta = angular_sep_deg(l, b, DPOLE_L, DPOLE_B)
                    galaxies.append({
                        'name': row.get('name', '?'),
                        'ra': ra, 'dec': dec, 'l': l, 'b': b,
                        'theta_D': theta, 'd_mpc': d,
                        'mod0': None
                    })
                except (ValueError, TypeError):
                    pass

# ---------------------------------------------------------------
# Thales-Berechnung + Kegel-Filter
# ---------------------------------------------------------------
KEGEL = 45.0  # Grad Oeffnungswinkel um D-Pol

dpol_objects = []
for g in galaxies:
    theta = g['theta_D']
    d     = g['d_mpc']
    if theta > KEGEL or theta < 1.0:  # theta < 1 Grad: Formel instabil
        continue
    if d <= 0:
        continue
    cos_t = math.cos(math.radians(theta))
    if cos_t < 0.01:
        continue
    r_thales = d / cos_t
    g['R_Thales'] = r_thales
    dpol_objects.append(g)

# Sortiere nach Distanz
dpol_objects.sort(key=lambda x: x['d_mpc'])

print(f"\nObjekte im D-Pol-Kegel ({KEGEL}°): {len(dpol_objects)}")

# ---------------------------------------------------------------
# Ausgabe
# ---------------------------------------------------------------
lines_out = []
lines_out.append("="*72)
lines_out.append("HYPERLEDA THALES-ANALYSE: Vakuumschalen-Radius")
lines_out.append(f"D-Pol: l={DPOLE_L}°, b={DPOLE_B}°  |  Kegel: {KEGEL}°")
lines_out.append(f"Formel: R_Thales = d_Mpc / cos(theta_D)")
lines_out.append(f"HTM r_s = {r_s:.1f} Mpc")
lines_out.append("="*72)
lines_out.append(f"\n{'Name':<22} {'d[Mpc]':>9} {'theta':>7} {'R_Thales':>11}")
lines_out.append("-"*55)

r_values = []
for g in dpol_objects:
    r  = g['R_Thales']
    r_values.append(r)
    lines_out.append(f"  {g['name']:<20} {g['d_mpc']:>9.1f} {g['theta_D']:>6.1f}° {r:>11.1f}")

if r_values:
    import statistics as st
    med = st.median(r_values)
    avg = st.mean(r_values)
    std = st.stdev(r_values) if len(r_values) > 1 else 0
    
    # Cluster um Median ±50%
    cluster = [r for r in r_values if 0.5*med < r < 1.5*med]
    
    lines_out.append("-"*55)
    lines_out.append(f"\nSTATISTIK ({len(r_values)} Objekte):")
    lines_out.append(f"  Median R_Thales : {med:>10.1f} Mpc")
    lines_out.append(f"  Mittel          : {avg:>10.1f} Mpc")
    lines_out.append(f"  Std             : {std:>10.1f} Mpc")
    lines_out.append(f"  Std/Median      : {std/med if med>0 else 0:>10.3f}")
    lines_out.append(f"  Cluster ±50%    : {len(cluster)}/{len(r_values)} Objekte")
    lines_out.append(f"\n  HTM r_s         : {r_s:>10.1f} Mpc")
    lines_out.append(f"  r_s / Median    : {r_s/med if med>0 else 0:>10.2f}")
    
    # Abstandsdekaden-Analyse
    lines_out.append("\nNACH DISTANZ-DEKADE:")
    dekaden = [(0,10),(10,50),(50,200),(200,500),(500,5000)]
    for d_min, d_max in dekaden:
        gruppe = [g for g in dpol_objects if d_min < g['d_mpc'] <= d_max]
        if gruppe:
            rv = [g['R_Thales'] for g in gruppe]
            lines_out.append(f"  d={d_min}-{d_max} Mpc: N={len(rv):3d}  "
                             f"Median-R={st.median(rv):8.1f}  "
                             f"Std={st.stdev(rv) if len(rv)>1 else 0:8.1f}")

# Naechste 10 Objekte (lokale Schalen-Struktur)
lines_out.append(f"\nNAECHSTE 10 IM KEGEL:")
for g in dpol_objects[:10]:
    lines_out.append(f"  {g['name']:<20} d={g['d_mpc']:.1f} Mpc  theta={g['theta_D']:.1f}°  R={g['R_Thales']:.1f} Mpc")

# ---------------------------------------------------------------
output_str = '\n'.join(lines_out)
print(output_str)

os.makedirs(os.path.join(RESULTS, "catalogs"), exist_ok=True)
with open(OUT_TXT, 'w', encoding='utf-8') as f:
    f.write(output_str)

# CSV speichern
if dpol_objects:
    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['name','l','b','theta_D','d_mpc','mod0','R_Thales'])
        writer.writeheader()
        writer.writerows(dpol_objects)
    print(f"\nCSV gespeichert: {OUT_CSV}")

print(f"Bericht gespeichert: {OUT_TXT}")
