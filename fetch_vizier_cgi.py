"""
Download van den Bosch 2016 (J/ApJS/228/14) via VizieR CGI (not TAP).
"""
import os, urllib.request, urllib.parse, sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "catalogs")
os.makedirs(RESULTS, exist_ok=True)

def try_vizier_cgi(catid, table, max_rows=1000, label=''):
    # VizieR CGI interface (not TAP)
    params = {
        '-source': f'{catid}/{table}',
        '-out.max': str(max_rows),
        '-out.form': 'CSV',
        '-out': '*',
        '-sort': '_r',
    }
    url = 'https://vizier.cds.unistra.fr/cgi-bin/VizieR-6?' + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/metageometra'})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read().decode('utf-8', 'replace')
        lines = data.strip().splitlines()
        # Remove comment lines starting with #
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        print(f'  {label or catid+"/"+table}: {len(data_lines)} data lines')
        if data_lines:
            print(f'  Headers: {data_lines[0][:100]}')
        return data
    except Exception as e:
        print(f'  {label or catid}/{table}: FAILED {e}')
        return None

def try_vizier_cgi_v3(catid, table, max_rows=500):
    # Different VizieR endpoint (v3)
    params = {
        '-source': f'{catid}/{table}',
        '-out.max': str(max_rows),
        '-out.form': 'CSV',
    }
    url = 'https://vizier.cds.unistra.fr/cgi-bin/VizieR-3?' + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/metageometra'})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read().decode('utf-8', 'replace')
        lines = [l for l in data.strip().splitlines() if l.strip() and not l.startswith('#')]
        print(f'  v3 {catid}/{table}: {len(lines)} lines | {lines[0][:80] if lines else "empty"}')
        return data
    except Exception as e:
        print(f'  v3 FAILED: {e}')
        return None

print('=== Probing VizieR CGI endpoints ===')
for catid, table, label in [
    ('J/ApJS/228/14', 'table1', 'vdBosch2016 t1'),
    ('J/ApJS/228/14', 'table2', 'vdBosch2016 t2'),
    ('J/ApJS/228/14', 'table3', 'vdBosch2016 t3'),
    ('J/ApJS/228/14', 'bhcat', 'vdBosch2016 bhcat'),
    ('J/ApJ/831/134', 'table3', 'Thomas2016 t3'),
    ('J/ARA+A/51/511', 'table2', 'KormendyHo2013'),
    ('J/ApJ/764/184', 'table1', 'McConnellMa2013'),
]:
    data = try_vizier_cgi(catid, table, max_rows=500, label=label)
    if data:
        lines = [l for l in data.strip().splitlines() if l.strip() and not l.startswith('#')]
        if len(lines) > 2:
            fname = f'viz_{label.replace(" ","_").replace("/","_")}.csv'
            out = os.path.join(RESULTS, fname)
            # Save only data lines (strip comments)
            with open(out, 'w', encoding='utf-8') as f:
                for line in data.strip().splitlines():
                    if not line.startswith('#'):
                        f.write(line + '\n')
            print(f'  Saved {len(lines)-1} data rows to {fname}')
        else:
            print(f'  Too few rows: {len(lines)}')
