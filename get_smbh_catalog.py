"""
Fetch larger SMBH catalog via VizieR vizquery (not TAP) and SIMBAD.
Targets:
 1. Graham & Scott 2013 (J/ApJ/764/151) - all dynamical SMBH
 2. Kormendy & Ho 2013 compilation
 3. Fill coordinates via SIMBAD

Saves to results/catalogs/smbh_large.csv
"""
import os, math, csv, time, io, urllib.request, urllib.parse, sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "catalogs")
os.makedirs(RESULTS, exist_ok=True)

D_POLE_L, D_POLE_B = 305.0, 25.0
THETA_0 = 58.65

def eq_to_gal(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    ra_ngp = math.radians(192.859508)
    dec_ngp = math.radians(27.128336)
    l_ncp = math.radians(122.932)
    sin_b = (math.sin(dec)*math.sin(dec_ngp) +
             math.cos(dec)*math.cos(dec_ngp)*math.cos(ra-ra_ngp))
    b = math.asin(max(-1., min(1., sin_b)))
    x = math.cos(dec)*math.sin(ra-ra_ngp)
    y = (math.sin(dec)*math.cos(dec_ngp) -
         math.cos(dec)*math.sin(dec_ngp)*math.cos(ra-ra_ngp))
    return math.degrees(l_ncp - math.atan2(x, y)) % 360, math.degrees(b)

def theta_dpole(ra, dec):
    l, b = eq_to_gal(ra, dec)
    l1r, b1r = math.radians(D_POLE_L), math.radians(D_POLE_B)
    l2r, b2r = math.radians(l), math.radians(b)
    cos_c = math.sin(b1r)*math.sin(b2r) + math.cos(b1r)*math.cos(b2r)*math.cos(l1r-l2r)
    return math.degrees(math.acos(max(-1., min(1., cos_c))))

def fetch(url, data=None, timeout=30, method='GET'):
    req = urllib.request.Request(url,
        data=data,
        headers={'User-Agent': 'Python/metageometra-catalog-builder',
                 'Content-Type': 'application/x-www-form-urlencoded'},
        method=method)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')

def simbad_coord(name):
    """Return (ra_deg, dec_deg) for galaxy name via SIMBAD, or None."""
    script = f'output console=off script=off\nformat object "%COO(d6;A)|%COO(d6;D)"\nquery id {name}\n'
    data = urllib.parse.urlencode({'script': script}).encode()
    url = 'http://simbad.cds.unistra.fr/simbad/sim-script'
    try:
        raw = fetch(url, data=data, timeout=15, method='POST')
        for line in raw.splitlines():
            line = line.strip()
            if '|' in line and not line.startswith(':'):
                parts = [p.strip() for p in line.split('|')]
                if len(parts) >= 2:
                    try:
                        return float(parts[0]), float(parts[1])
                    except ValueError:
                        pass
    except Exception:
        pass
    return None

# ── VizieR: Kormendy & Ho 2013 Compilation (J/ARA+A/51/511) ─────────────────
def try_tap(q, label, timeout=45):
    url = ('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
           '?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=' + urllib.parse.quote(q))
    try:
        raw = fetch(url, timeout=timeout)
        lns = raw.strip().splitlines()
        if len(lns) > 1 and 'error' not in raw[:300].lower():
            print(f'  {label}: {len(lns)-1} rows | {lns[0][:80]}')
            return raw
        else:
            print(f'  {label}: no rows or error — {raw[:200]}')
    except Exception as e:
        print(f'  {label} FAILED: {e}')
    return None

print('=== VizieR TAP catalog probe ===')
for q, lbl in [
    ('SELECT TOP 5 * FROM "J/ARA+A/51/511/table2"', 'Kormendy+Ho2013'),
    ('SELECT TOP 5 * FROM "J/ARA+A/51/511/table1"', 'Kormendy+Ho2013 t1'),
    ('SELECT TOP 5 * FROM "J/ApJ/764/151/table1"', 'Graham+Scott2013'),
    ('SELECT TOP 5 * FROM "J/ApJ/764/151/table2"', 'Graham+Scott2013 t2'),
    ('SELECT TOP 5 * FROM "J/ApJS/228/14/table3"', 'vdBosch2016 t3'),
    ('SELECT TOP 5 * FROM "J/ApJS/228/14/table5"', 'vdBosch2016 t5'),
    ('SELECT TOP 5 * FROM "J/ApJS/228/14/mb"', 'vdBosch2016 mb'),
    ('SELECT TOP 5 * FROM "J/ApJS/228/14/table1b"', 'vdBosch2016 t1b'),
    ('SELECT TOP 5 * FROM "J/ApJS/228/14/bhcat"', 'vdBosch2016 bhcat'),
    ('SELECT TOP 5 * FROM "J/MNRAS/460/3119/table2"', 'Sahu2019 A'),
    ('SELECT TOP 5 * FROM "J/ApJ/887/10/table2"', 'Sahu2019 B'),
    ('SELECT TOP 5 * FROM "J/ApJ/831/134/table3"', 'Thomas2016'),
]:
    result = try_tap(q, lbl)
    if result and len(result.strip().splitlines()) > 1:
        print(f'  *** FOUND: {lbl} ***')
        outf = os.path.join(RESULTS, f'tap_{lbl.replace(" ","_").replace("+","")}.csv')
        with open(outf, 'w', encoding='utf-8') as f:
            f.write(result)
        print(f'  Saved: {outf}')
