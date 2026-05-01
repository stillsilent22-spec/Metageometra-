"""Prepare CF3 full-sky theta_D catalog for OT-40 permutation test."""
import math, csv, os

BASE = os.path.dirname(os.path.abspath(__file__))

# D-pole: l=305°, b=+25°
L_DP = math.radians(305.0); B_DP = math.radians(25.0)
dp_x = math.cos(B_DP) * math.cos(L_DP)
dp_y = math.cos(B_DP) * math.sin(L_DP)
dp_z = math.sin(B_DP)

IN_TSV = os.path.join(BASE, "results", "catalogs", "cf3_all.tsv")
OUT_CSV = os.path.join(BASE, "results", "catalogs", "cf3_fullsky_thetaD.csv")

rows = []
with open(IN_TSV, encoding='utf-8') as f:
    lines = f.readlines()

# Find the header line (first non-comment, non-dash line)
headers = None
skip = 0
for i, l in enumerate(lines):
    if l.startswith('#') or l.startswith('-') or l.strip() == '':
        continue
    if headers is None:
        headers = l.strip().split('\t')
        skip = i + 1  # skip header + units + dashes = header+2 lines
        print(f"Header at line {i}: {headers[:12]}")
        break

# Skip also units and dashes lines
data_start = skip
while data_start < len(lines):
    l = lines[data_start].strip()
    if l.startswith('-') or l == '' or all(c in '-\t ' for c in l):
        data_start += 1
        continue
    # Check if it's a units line (first column should be numeric for data)
    first_col = l.split('\t')[0].strip()
    if first_col.replace(' ', '') in ('', 'deg', 'mag', 'Mpc', 'km/s', '[Lsun]', 'TMsun'):
        data_start += 1
        continue
    break

try:
    i_glon = headers.index('GLON')
    i_glat = headers.index('GLAT')
    i_dist = headers.index('<Dist>')
except ValueError as e:
    print(f"Column not found: {e}")
    print("Available columns:", headers)
    raise

count = 0
for line in lines[data_start:]:
    parts = line.strip().split('\t')
    if len(parts) <= max(i_glon, i_glat):
        continue
    try:
        glon_deg = float(parts[i_glon].strip())
        glat_deg = float(parts[i_glat].strip())
        dist_mpc = float(parts[i_dist].strip()) if parts[i_dist].strip() else 0.
    except ValueError:
        continue

    # Convert to unit vector
    l_rad = math.radians(glon_deg); b_rad = math.radians(glat_deg)
    gx = math.cos(b_rad) * math.cos(l_rad)
    gy = math.cos(b_rad) * math.sin(l_rad)
    gz = math.sin(b_rad)

    # Angular distance from D-pole
    dot = min(1., max(-1., gx*dp_x + gy*dp_y + gz*dp_z))
    theta_D = math.degrees(math.acos(dot))

    rows.append({'glon': glon_deg, 'glat': glat_deg, 'dist_mpc': dist_mpc, 'theta_D': theta_D})
    count += 1

print(f"Parsed {count} CF3 groups")

with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['glon', 'glat', 'dist_mpc', 'theta_D'])
    writer.writeheader()
    writer.writerows(rows)

print(f"Saved: {OUT_CSV}")

import numpy as np
thetas = [r['theta_D'] for r in rows]
thetas = np.array(thetas)
print(f"theta_D range: {thetas.min():.1f} - {thetas.max():.1f} deg")
hist, bins = np.histogram(thetas, bins=[0, 20, 40, 60, 80, 100, 120, 140, 160, 180])
print(f"Distribution: {list(hist)}")
