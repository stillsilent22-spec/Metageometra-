"""Download CF3 and SPARC catalogs via multiple fallback methods."""
import urllib.request
import urllib.parse
import os, sys

def try_url(url, dest, label, decode='utf-8'):
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = r.read()
        with open(dest, 'wb') as f:
            f.write(data)
        lines = data.decode(decode, 'replace').strip().split('\n')
        data_lines = [l for l in lines if l.strip() and not l.startswith('#') and not l.startswith('-')]
        print(f"  {label}: {len(lines)} Zeilen total, {len(data_lines)} Datenzeilen")
        for l in data_lines[:3]:
            print(f"    {repr(l[:120])}")
        return len(data_lines) > 0
    except Exception as e:
        print(f"  FAIL {label}: {e}")
        return False

os.makedirs('results/catalogs', exist_ok=True)

print("=== CF3 (Cosmicflows-3) ===")
# Method 1: VizieR TAP with correct table name
query1 = 'SELECT TOP 200 RAJ2000,DEJ2000,cz,e_cz,Dist,e_Dist,l,b FROM "J/AJ/152/50/table1"'
params1 = urllib.parse.urlencode({
    'REQUEST': 'doQuery',
    'LANG': 'ADQL',
    'FORMAT': 'csv',
    'QUERY': query1
})
url_tap = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync?' + params1
try_url(url_tap, 'results/catalogs/cf3_tap.csv', 'CF3 TAP ADQL')

# Method 2: VizieR catalog query form
url_viz2 = ('https://vizier.cds.unistra.fr/viz-bin/votable?'
            '-source=J/AJ/152/50/table1&-out.all&-out.max=200')
try_url(url_viz2, 'results/catalogs/cf3_votable.xml', 'CF3 VOTable')

# Method 3: VizieR plain form
url_viz3 = ('https://vizier.cds.unistra.fr/viz-bin/asu-tsv?'
            '-source=J%2FAJ%2F152%2F50%2Ftable1'
            '&-out.all&-out.max=200&-out.type=sorted')
try_url(url_viz3, 'results/catalogs/cf3_asu2.tsv', 'CF3 ASU-TSV v2')

print("\n=== SPARC Rotationskurven ===")
# Method 1: VizieR RotCur table
url_sparc1 = ('https://vizier.cds.unistra.fr/viz-bin/asu-tsv?'
              '-source=J%2FAJ%2F152%2F157%2Frotcur'
              '&-out=Galaxy,Rad,Vobs,e_Vobs,Vgas,Vdisk,Vbul'
              '&-out.max=50000')
try_url(url_sparc1, 'results/catalogs/sparc_rotcur.tsv', 'SPARC RotCur TSV')

# Method 2: VOTable
url_sparc2 = ('https://vizier.cds.unistra.fr/viz-bin/votable?'
              '-source=J/AJ/152/157/rotcur&-out.all&-out.max=50000')
try_url(url_sparc2, 'results/catalogs/sparc_rotcur.xml', 'SPARC RotCur VOTable')

# Method 3: Direct from Stacy McGaugh's page
url_sparc3 = 'http://astroweb.cwru.edu/SPARC/SPARC_Lelli2016c.mrt'
try_url(url_sparc3, 'results/catalogs/sparc_mrt.txt', 'SPARC MRT McGaugh')

# Method 4: Zenodo (Lelli 2016)
url_sparc4 = 'https://zenodo.org/record/5562847/files/Lelli2016c_SPARC_Rotcur.mrt'
try_url(url_sparc4, 'results/catalogs/sparc_zenodo.mrt', 'SPARC Zenodo MRT')

print("\nDone.")
