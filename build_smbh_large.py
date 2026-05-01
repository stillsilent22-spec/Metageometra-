"""
Fetch larger SMBH catalogs from VizieR and build extended position catalog.
Downloads Thomas+2016 (J/ApJ/831/134) and probes other accessible catalogs.
"""
import os, math, csv, io, urllib.request, urllib.parse, sys, time

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

def tap_fetch(q, timeout=60):
    url = ('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
           '?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=' + urllib.parse.quote(q))
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/metageometra'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')

def simbad_coord(name):
    script = (f'output console=off script=off\n'
              f'format object "%COO(d6;A)|%COO(d6;D)"\n'
              f'query id {name}\n')
    data = urllib.parse.urlencode({'script': script}).encode()
    url = 'http://simbad.cds.unistra.fr/simbad/sim-script'
    req = urllib.request.Request(url, data=data,
            headers={'User-Agent': 'Python/metageometra',
                     'Content-Type': 'application/x-www-form-urlencoded'})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            raw = r.read().decode('utf-8', 'replace')
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

print('=== Step 1: Download Thomas+2016 full catalog ===')
try:
    q_thomas_full = 'SELECT recno,Name,SimbadName,_RA,_DE,logBHMass,E_logBHMass,e_logBHMass,logsigma,Dist FROM "J/ApJ/831/134/table3" MAXREC=5000'
    thomas_raw = tap_fetch(q_thomas_full)
    lines_t = thomas_raw.strip().splitlines()
    print(f'Thomas+2016 full: {len(lines_t)-1} objects')
    out_t = os.path.join(RESULTS, 'thomas2016_full.csv')
    with open(out_t, 'w', encoding='utf-8') as f:
        f.write(thomas_raw)
    print(f'Saved: {out_t}')
    print(f'Headers: {lines_t[0]}')
    if len(lines_t) > 1:
        print(f'Row 1:   {lines_t[1]}')
except Exception as e:
    print(f'Thomas+2016 download error: {e}')
    thomas_raw = None

# ── Step 2: Probe more VizieR catalogs for SMBH positions ─────────────────────
print('\n=== Step 2: Probe additional VizieR SMBH catalogs ===')
probe_queries = [
    ('SELECT TOP 5 * FROM "J/MNRAS/460/3119/table1"', 'Savorgnan2016 t1'),
    ('SELECT TOP 5 * FROM "J/MNRAS/460/3119/table2"', 'Savorgnan2016 t2'),
    ('SELECT TOP 5 * FROM "J/MNRAS/460/3119/table3"', 'Savorgnan2016 t3'),
    ('SELECT TOP 5 * FROM "J/MNRAS/460/3119/bhcat"', 'Savorgnan2016 bhcat'),
    ('SELECT TOP 5 * FROM "J/ApJ/831/134/table1"', 'Thomas2016 t1'),
    ('SELECT TOP 5 * FROM "J/ApJ/831/134/table2"', 'Thomas2016 t2'),
    ('SELECT TOP 5 * FROM "J/ApJ/831/134/table4"', 'Thomas2016 t4'),
    ('SELECT TOP 5 * FROM "J/ApJ/831/134/table5"', 'Thomas2016 t5'),
    ('SELECT TOP 5 * FROM "J/ApJ/831/134/table6"', 'Thomas2016 t6'),
    ('SELECT TOP 5 * FROM "J/MNRAS/451/3413/table1"', 'Saglia2016 t1'),
    ('SELECT TOP 5 * FROM "J/MNRAS/451/3413/table2"', 'Saglia2016 t2'),
    ('SELECT TOP 5 * FROM "J/ApJ/813/82/table2"', 'McConnell2015'),
    ('SELECT TOP 5 * FROM "J/A+A/578/A110/table1"', 'Krajnovic2015'),
    ('SELECT TOP 5 * FROM "J/MNRAS/486/4726/table1"', 'Nyland2020'),
    ('SELECT TOP 5 * FROM "J/ApJS/245/25/table3"', 'Nguyen2019 t3'),
    ('SELECT TOP 5 * FROM "J/ApJS/245/25/bhcat"', 'Nguyen2019 bhcat'),
]

found_catalogs = []
for q, lbl in probe_queries:
    try:
        raw = tap_fetch(q, timeout=45)
        lns = raw.strip().splitlines()
        if len(lns) > 1 and 'error' not in raw[:300].lower():
            print(f'  FOUND {lbl}: {len(lns)-1} rows | {lns[0][:80]}')
            found_catalogs.append((lbl, q.replace('TOP 5', 'MAXREC=1000'), lns[0]))
        else:
            print(f'  skip {lbl}: too short')
    except Exception as e:
        print(f'  fail {lbl}: {e}')

# ── Step 3: Download all found catalogs ────────────────────────────────────────
print('\n=== Step 3: Download full versions of found catalogs ===')
for lbl, q_full, hdr in found_catalogs:
    try:
        raw = tap_fetch(q_full, timeout=120)
        lns = raw.strip().splitlines()
        fname = f'smbh_{lbl.replace(" ","_").replace("+","_").lower()}.csv'
        out = os.path.join(RESULTS, fname)
        with open(out, 'w', encoding='utf-8') as f:
            f.write(raw)
        print(f'  Downloaded {lbl}: {len(lns)-1} rows -> {fname}')
    except Exception as e:
        print(f'  Download failed {lbl}: {e}')

# ── Step 4: Build combined high-count catalog with coordinates ─────────────────
print('\n=== Step 4: Build combined SMBH catalog ===')

# Load Thomas+2016 as base (already has _RA, _DE)
all_objects = []  # (name, ra, dec, logM, source)
seen_names = set()

def norm_name(n): return str(n).strip().upper().replace(' ', '').replace('-', '')

if thomas_raw:
    reader = csv.DictReader(io.StringIO(thomas_raw))
    for row in reader:
        name = (row.get('SimbadName') or row.get('Name') or '').strip()
        if not name or name == ' ':
            continue
        ra_s = row.get('_RA', '').strip()
        de_s = row.get('_DE', '').strip()
        logm = row.get('logBHMass', '').strip()
        if not ra_s or not de_s or not logm:
            continue
        try:
            ra = float(ra_s)
            dec = float(de_s)
            logM = float(logm)
        except ValueError:
            continue
        nn = norm_name(name)
        if nn not in seen_names:
            seen_names.add(nn)
            all_objects.append({'Name': name, 'RA_deg': ra, 'Dec_deg': dec,
                                'logMbh': logM, 'source': 'Thomas2016'})

print(f'  Thomas2016: {len(all_objects)} objects with coords')

# Load our existing smbh_extended.csv
extended_csv = os.path.join(RESULTS, 'smbh_extended.csv')
base_before = len(all_objects)
with open(extended_csv, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row.get('Name', '').strip()
        nn = norm_name(name)
        if nn in seen_names:
            continue
        try:
            ra = float(row['RA_deg'])
            dec = float(row['Dec_deg'])
            logM = float(row['logMbh'])
        except (ValueError, KeyError):
            continue
        seen_names.add(nn)
        all_objects.append({'Name': name, 'RA_deg': ra, 'Dec_deg': dec,
                            'logMbh': logM, 'source': row.get('source', 'existing')})
print(f'  After merge with smbh_extended: +{len(all_objects)-base_before} = {len(all_objects)} total')

# Compute HTM angles for all objects
for obj in all_objects:
    td = theta_dpole(obj['RA_deg'], obj['Dec_deg'])
    obj['theta_dpole'] = td
    obj['theta_apole'] = 180.0 - td

# Save combined catalog
out_combined = os.path.join(RESULTS, 'smbh_large.csv')
fieldnames = ['Name', 'RA_deg', 'Dec_deg', 'logMbh', 'source', 'theta_dpole', 'theta_apole']
with open(out_combined, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
    w.writeheader()
    w.writerows(all_objects)

print(f'\nFinal catalog: {len(all_objects)} SMBH objects')
print(f'Saved to: {out_combined}')
