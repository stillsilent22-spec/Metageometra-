"""
Download SMBH catalogs from VizieR TAP using correct ADQL TOP syntax.
"""
import os, math, csv, io, urllib.request, urllib.parse, sys

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

def tap_fetch(q, timeout=90):
    url = ('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
           '?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=' + urllib.parse.quote(q))
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/metageometra'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')

print('=== Downloading SMBH catalogs via VizieR TAP (SELECT TOP syntax) ===')

all_objects = []
seen = set()

def norm(n): return str(n).strip().upper().replace(' ','').replace('-','').replace('_','')

# ── 1. Thomas+2016 (J/ApJ/831/134) — SMBH measurements ──────────────────────
print('\n--- Thomas+2016 (J/ApJ/831/134/table3) ---')
for tbl in ['table3', 'table2', 'table4']:
    try:
        q = f'SELECT TOP 1000 * FROM "J/ApJ/831/134/{tbl}"'
        raw = tap_fetch(q)
        lines = raw.strip().splitlines()
        print(f'  {tbl}: {len(lines)-1} rows | Headers: {lines[0][:80]}')
        if len(lines) > 1:
            reader = csv.DictReader(io.StringIO(raw))
            n_added = 0
            for row in reader:
                name = (row.get('SimbadName') or row.get('Name') or '').strip()
                ra_s = (row.get('_RA') or row.get('RAJ2000') or '').strip()
                de_s = (row.get('_DE') or row.get('DEJ2000') or '').strip()
                logm = (row.get('logBHMass') or '').strip()
                if not ra_s or not de_s or not logm or not name or name == ' ':
                    continue
                try:
                    ra, dec, logM = float(ra_s), float(de_s), float(logm)
                    nn = norm(name)
                    if nn not in seen:
                        seen.add(nn)
                        all_objects.append({'Name': name, 'RA_deg': ra, 'Dec_deg': dec,
                                            'logMbh': logM, 'source': 'Thomas2016'})
                        n_added += 1
                except ValueError:
                    pass
            print(f'  Added {n_added} objects with coords from {tbl}')
    except Exception as e:
        print(f'  {tbl}: {e}')

print(f'After Thomas2016: {len(all_objects)} objects')

# ── 2. Probe more accessible catalogs ────────────────────────────────────────
print('\n--- Probing other catalogs ---')
extra_probes = [
    ('J/ApJ/831/134/table1', 'Thomas2016-t1'),
    ('J/ApJ/831/134/table5', 'Thomas2016-t5'),
    ('J/ApJ/831/134/table6', 'Thomas2016-t6'),
    ('J/MNRAS/451/3413/table2', 'Saglia2016-t2'),
    ('J/MNRAS/451/3413/table3', 'Saglia2016-t3'),
    ('J/ApJ/813/82/table2', 'McConnell2015-t2'),
    ('J/A+A/578/A110/table1', 'Krajnovic2015-t1'),
    ('J/MNRAS/460/3119/table2', 'Savorgnan2016-t2'),
    ('J/MNRAS/460/3119/table1', 'Savorgnan2016-t1'),
    ('J/ApJS/245/25/table3', 'Nguyen2019-t3'),
]
for tbl_id, label in extra_probes:
    try:
        q = f'SELECT TOP 5 * FROM "{tbl_id}"'
        raw = tap_fetch(q, timeout=30)
        lines = raw.strip().splitlines()
        if len(lines) > 1 and 'error' not in raw[:200].lower():
            print(f'  FOUND {label}: cols={lines[0][:80]}')
        else:
            print(f'  {label}: empty/error')
    except Exception as e:
        print(f'  {label}: {e}')

# ── 3. Merge with existing smbh_extended.csv ─────────────────────────────────
print('\n--- Merging with smbh_extended.csv ---')
ext_path = os.path.join(RESULTS, 'smbh_extended.csv')
n_before = len(all_objects)
with open(ext_path, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        name = row.get('Name','').strip()
        nn = norm(name)
        if nn in seen:
            continue
        try:
            ra, dec, logM = float(row['RA_deg']), float(row['Dec_deg']), float(row['logMbh'])
            seen.add(nn)
            all_objects.append({'Name': name, 'RA_deg': ra, 'Dec_deg': dec,
                                'logMbh': logM, 'source': row.get('source','')})
        except (ValueError, KeyError):
            pass
print(f'+{len(all_objects)-n_before} from smbh_extended.csv → total: {len(all_objects)}')

# ── 4. Compute HTM angles ────────────────────────────────────────────────────
for obj in all_objects:
    td = theta_dpole(obj['RA_deg'], obj['Dec_deg'])
    obj['theta_dpole'] = td
    obj['theta_apole'] = 180.0 - td

# ── 5. Save combined catalog ─────────────────────────────────────────────────
out = os.path.join(RESULTS, 'smbh_large.csv')
cols = ['Name','RA_deg','Dec_deg','logMbh','source','theta_dpole','theta_apole']
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
    w.writeheader()
    w.writerows(all_objects)

print(f'\nFinal catalog: {len(all_objects)} SMBH objects')
print(f'Saved: {out}')
