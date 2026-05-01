"""
Download van den Bosch 2016 full SMBH catalog from VizieR TAP.
Saves to results/catalogs/vdb2016_full.csv
"""
import urllib.request, urllib.parse, os, csv, io

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "catalogs")
os.makedirs(RESULTS, exist_ok=True)

from urllib.error import HTTPError

def try_vizier_query(q, label):
    tap_url = ('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
               '?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=' + urllib.parse.quote(q))
    req = urllib.request.Request(tap_url, headers={'User-Agent': 'Python/metageometra'})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read().decode('utf-8')
        lines = data.strip().splitlines()
        if len(lines) > 1 and 'error' not in data[:200].lower():
            print(f"  {label}: {len(lines)-1} rows | Headers: {lines[0][:80]}")
            return data
        else:
            print(f"  {label}: response too short or error: {data[:200]}")
    except HTTPError as e:
        print(f"  {label}: HTTP {e.code}")
    except Exception as e:
        print(f"  {label}: {e}")
    return None

print("Probing VizieR for SMBH catalogs...")
data = None

# Try various van den Bosch 2016 table names
for q, lbl in [
    ('SELECT * FROM "J/ApJS/228/14/table1" MAXREC=5', 'vdB2016 table1'),
    ('SELECT * FROM "J/ApJS/228/14/table2" MAXREC=5', 'vdB2016 table2'),
    ('SELECT * FROM "J/ApJS/228/14/table3" MAXREC=5', 'vdB2016 table3'),
    ('SELECT * FROM "J/ApJS/228/14/table4" MAXREC=5', 'vdB2016 table4'),
    ('SELECT * FROM "J/ApJS/228/14/galaxies" MAXREC=5', 'vdB2016 galaxies'),
    # McConnell & Ma 2013
    ('SELECT * FROM "J/ApJ/764/184/table1" MAXREC=5', 'McConnell2013 table1'),
    ('SELECT * FROM "J/ApJ/764/184/table2" MAXREC=5', 'McConnell2013 table2'),
    # Kormendy & Ho 2013
    ('SELECT * FROM "J/ARA+A/51/511/table2" MAXREC=5', 'KH2013 table2'),
    # Sahu 2019
    ('SELECT * FROM "J/ApJ/887/10/table1" MAXREC=5', 'Sahu2019 table1'),
]:
    data = try_vizier_query(q, lbl)
    if data:
        print(f"SUCCESS: {lbl}")
        break

if not data:
    print("  All probes failed, trying SIMBAD-based approach")
lines = data.strip().splitlines() if data else []

out = os.path.join(RESULTS, "vdb2016_full.csv")
with open(out, 'w', encoding='utf-8') as f:
    f.write(data)
print(f"  Saved to: {out}")
