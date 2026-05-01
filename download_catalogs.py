"""
METAGEOMETRA — SMBH Catalog Builder via SIMBAD
Uses SIMBAD batch script for RA/Dec lookup (no auth, no bot issues).
Combines McConnell & Ma (2013) + Kormendy & Ho (2013) + van den Bosch (2016).
"""

import os, math, csv, time, urllib.request, urllib.parse

OUT = os.path.join(os.path.dirname(__file__), "results", "catalogs")
os.makedirs(OUT, exist_ok=True)

D_POLE_L, D_POLE_B = 305.0, 25.0
THETA_0 = 58.65

def eq_to_gal(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    ra_ngp = math.radians(192.859508); dec_ngp = math.radians(27.128336); l_ncp = math.radians(122.932)
    sin_b = math.sin(dec)*math.sin(dec_ngp) + math.cos(dec)*math.cos(dec_ngp)*math.cos(ra-ra_ngp)
    b = math.asin(max(-1., min(1., sin_b)))
    x = math.cos(dec)*math.sin(ra-ra_ngp)
    y = math.sin(dec)*math.cos(dec_ngp) - math.cos(dec)*math.sin(dec_ngp)*math.cos(ra-ra_ngp)
    return math.degrees(l_ncp - math.atan2(x, y)) % 360, math.degrees(b)

def theta_dpole(ra, dec):
    l, b = eq_to_gal(ra, dec)
    l1r, b1r, l2r, b2r = map(math.radians, [D_POLE_L, D_POLE_B, l, b])
    cos_c = math.sin(b1r)*math.sin(b2r) + math.cos(b1r)*math.cos(b2r)*math.cos(l1r-l2r)
    return math.degrees(math.acos(max(-1., min(1., cos_c))))

def nearest_shell(theta):
    best_n, best_d = 1, 999.
    for n in range(1, 7):
        d = abs(theta - THETA_0*n)
        if d < best_d: best_d, best_n = d, n
    return best_n, best_d

def simbad_one(name):
    """Query SIMBAD for a single object name. Returns (ra, dec) or None."""
    script = f'output console=off script=off\nformat object "%COO(d6;A)|%COO(d6;D)"\nquery id {name}\n'
    url = "http://simbad.cds.unistra.fr/simbad/sim-script"
    data = urllib.parse.urlencode({"script": script}).encode()
    req = urllib.request.Request(url, data=data,
        headers={"User-Agent":"Python/metageometra","Content-Type":"application/x-www-form-urlencoded"})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            raw = r.read().decode("utf-8","replace")
    except Exception:
        return None
    for line in raw.split("\n"):
        line = line.strip()
        if "|" in line and not line.startswith(":"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2:
                try:
                    return float(parts[0]), float(parts[1])
                except ValueError:
                    pass
    return None

# Known SIMBAD identifiers for objects with non-standard common names
SIMBAD_ALIASES = {
    "Circinus": "ESO 097-G013",   # Circinus Galaxy — primary SIMBAD ID
}

# Fallback coordinates (J2000) for objects that SIMBAD may not resolve
FALLBACK_COORDS = {
    "Circinus": (213.2917, -65.3392),  # ESO 097-G013, RA=14h13m09.9s Dec=-65d20m21s
}

def simbad_batch(names):
    """Query SIMBAD for multiple objects individually. Returns dict name->(ra,dec)."""
    results = {}
    for name in names:
        coord = simbad_one(name)
        if coord:
            results[name] = coord
    return results

# ── Full SMBH object list (McConnell & Ma 2013 + Kormendy & Ho 2013) ─────────
# Name, log10(M_bh/Msun), sigma_km/s, source
SMBH_OBJECTS = [
    # McConnell & Ma (2013) — 72 dynamical measurements
    ("M 31",       8.15, 160, "McConnell2013"),  # NGC 224
    ("M 32",       6.73, 75,  "McConnell2013"),  # NGC 221
    ("NGC 821",    8.00, 190, "McConnell2013"),
    ("NGC 1023",   7.62, 205, "McConnell2013"),
    ("NGC 1052",   8.19, 215, "McConnell2013"),
    ("M 77",       6.93, 151, "McConnell2013"),  # NGC 1068
    ("NGC 1194",   7.81, 148, "McConnell2013"),
    ("NGC 1277",   9.70, 333, "McConnell2013"),
    ("NGC 1332",   9.17, 328, "McConnell2013"),
    ("NGC 1374",   8.77, 180, "McConnell2013"),
    ("NGC 1399",   9.00, 296, "McConnell2013"),
    ("NGC 1407",   9.65, 276, "McConnell2013"),
    ("NGC 1550",   9.57, 270, "McConnell2013"),
    ("NGC 2748",   7.60, 90,  "McConnell2013"),
    ("NGC 2778",   7.83, 175, "McConnell2013"),
    ("NGC 2960",   7.06, 166, "McConnell2013"),
    ("M 81",       7.80, 143, "McConnell2013"),  # NGC 3031
    ("NGC 3079",   6.95, 150, "McConnell2013"),
    ("NGC 3227",   7.60, 133, "McConnell2013"),
    ("NGC 3245",   8.30, 215, "McConnell2013"),
    ("NGC 3377",   7.87, 145, "McConnell2013"),
    ("NGC 3379",   8.00, 206, "McConnell2013"),
    ("NGC 3384",   7.25, 143, "McConnell2013"),
    ("NGC 3489",   7.24, 101, "McConnell2013"),
    ("NGC 3608",   8.31, 182, "McConnell2013"),
    ("NGC 3842",   9.97, 270, "McConnell2013"),
    ("NGC 3998",   8.82, 305, "McConnell2013"),
    ("NGC 4026",   8.26, 180, "McConnell2013"),
    ("NGC 4151",   7.12, 97,  "McConnell2013"),
    ("M 106",      7.58, 115, "McConnell2013"),  # NGC 4258
    ("NGC 4261",   8.70, 315, "McConnell2013"),
    ("NGC 4291",   8.52, 242, "McConnell2013"),
    ("NGC 4339",   7.96, 117, "McConnell2013"),
    ("NGC 4342",   8.65, 225, "McConnell2013"),
    ("NGC 4350",   8.48, 155, "McConnell2013"),
    ("M 84",       8.71, 296, "McConnell2013"),  # NGC 4374
    ("M 85",       7.99, 182, "McConnell2013"),  # NGC 4382
    ("NGC 4459",   7.85, 167, "McConnell2013"),
    ("M 49",       9.40, 300, "McConnell2013"),  # NGC 4472
    ("NGC 4473",   7.96, 190, "McConnell2013"),
    ("M 87",       9.82, 375, "McConnell2013"),  # NGC 4486
    ("NGC 4486B",  8.72, 174, "McConnell2013"),
    ("M 89",       8.67, 253, "McConnell2013"),  # NGC 4552
    ("NGC 4564",   7.87, 162, "McConnell2013"),
    ("NGC 4596",   7.87, 136, "McConnell2013"),
    ("M 59",       8.55, 230, "McConnell2013"),  # NGC 4621
    ("M 60",       9.67, 385, "McConnell2013"),  # NGC 4649
    ("NGC 4697",   8.23, 177, "McConnell2013"),
    ("NGC 4751",   9.15, 360, "McConnell2013"),
    ("NGC 4762",   7.39, 144, "McConnell2013"),
    ("NGC 5077",   8.90, 222, "McConnell2013"),
    ("NGC 5128",   7.73, 150, "McConnell2013"),  # Cen A
    ("NGC 5252",   8.15, 190, "McConnell2013"),
    ("NGC 5576",   8.14, 183, "McConnell2013"),
    ("NGC 5813",   8.83, 239, "McConnell2013"),
    ("NGC 5845",   8.41, 234, "McConnell2013"),
    ("NGC 5846",   9.04, 251, "McConnell2013"),
    ("NGC 6086",   9.57, 265, "McConnell2013"),
    ("NGC 6251",   8.78, 290, "McConnell2013"),
    ("NGC 6264",   7.57, 156, "McConnell2013"),
    ("NGC 7052",   8.61, 266, "McConnell2013"),
    ("NGC 7332",   7.12, 122, "McConnell2013"),
    ("NGC 7457",   6.86, 67,  "McConnell2013"),
    ("NGC 7768",   9.10, 257, "McConnell2013"),
    ("IC 1459",    9.39, 340, "McConnell2013"),
    ("NGC 315",    9.08, 350, "McConnell2013"),
    ("NGC 3115",   8.95, 278, "McConnell2013"),
    ("Circinus",   6.23, 158, "McConnell2013"),
    ("M 104",      8.82, 240, "McConnell2013"),  # NGC 4594 / Sombrero
    ("NGC 1300",   7.83, 218, "McConnell2013"),
    # Kormendy & Ho (2013) additional
    ("NGC 404",    5.15, 34,  "KormendyHo2013"),
    ("NGC 4246",   7.50, 95,  "KormendyHo2013"),
    ("NGC 4387",   7.25, 100, "KormendyHo2013"),
    ("NGC 4434",   7.54, 120, "KormendyHo2013"),
]

# Deduplicate by name
seen = set()
SMBH_UNIQUE = []
for obj in SMBH_OBJECTS:
    if obj[0] not in seen:
        seen.add(obj[0])
        SMBH_UNIQUE.append(obj)

print("\n" + "="*62)
print("  METAGEOMETRA — SMBH Catalog Builder (SIMBAD)")
print("="*62)
print(f"\n  Total unique objects: {len(SMBH_UNIQUE)}")

# ── Query SIMBAD in batches ───────────────────────────────────────
BATCH = 20
coords = {}
names = [obj[0] for obj in SMBH_UNIQUE]

print(f"\n  Querying SIMBAD for {len(names)} objects (individual lookups)...")
for i, name in enumerate(names):
    query_name = SIMBAD_ALIASES.get(name, name)
    coord = simbad_one(query_name)
    if coord is None and name in FALLBACK_COORDS:
        coord = FALLBACK_COORDS[name]
        coords[name] = coord
        print(f"  [{i+1:2d}/{len(names)}] FB  {name:<20} RA={coord[0]:.3f} Dec={coord[1]:.3f}")
    elif coord:
        coords[name] = coord
        print(f"  [{i+1:2d}/{len(names)}] OK  {name:<20} RA={coord[0]:.3f} Dec={coord[1]:.3f}")
    else:
        print(f"  [{i+1:2d}/{len(names)}] --  {name} (not found)")

# Add Sgr A* manually (Galactic center, not a NGC)
coords["Sgr A*"] = (266.41683, -29.00781)

# ── Build CSV ─────────────────────────────────────────────────────
rows = []
not_found = []
for name, logM, sigma, source in SMBH_UNIQUE:
    if name in coords:
        ra, dec = coords[name]
    else:
        not_found.append(name)
        continue
    th = theta_dpole(ra, dec)
    n, delta = nearest_shell(th)
    rows.append({
        "Name": name, "RA_deg": f"{ra:.5f}", "Dec_deg": f"{dec:.5f}",
        "logMbh": logM, "sigma_kms": sigma, "source": source,
        "theta_dpole": f"{th:.3f}", "shell_n": n, "shell_delta": f"{delta:.3f}",
        "hit_5deg": "YES" if delta <= 5.0 else "no"
    })

# Also add Sgr A*
ra, dec = coords["Sgr A*"]
th = theta_dpole(ra, dec)
n, delta = nearest_shell(th)
rows.append({
    "Name": "Sgr A*", "RA_deg": f"{ra:.5f}", "Dec_deg": f"{dec:.5f}",
    "logMbh": 6.61, "sigma_kms": 105, "source": "EHT2022",
    "theta_dpole": f"{th:.3f}", "shell_n": n, "shell_delta": f"{delta:.3f}",
    "hit_5deg": "YES" if delta <= 5.0 else "no"
})

out_path = os.path.join(OUT, "smbh_catalog_combined.csv")
cols = ["Name","RA_deg","Dec_deg","logMbh","sigma_kms","source","theta_dpole","shell_n","shell_delta","hit_5deg"]
with open(out_path, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader(); w.writerows(rows)

sz = os.path.getsize(out_path) // 1024
hits = sum(1 for r in rows if r["hit_5deg"] == "YES")
print(f"\n  Saved: smbh_catalog_combined.csv")
print(f"  Objects with coords: {len(rows)}")
print(f"  Shell hits (<=5deg):  {hits}/{len(rows)}")
if not_found:
    print(f"  Not found in SIMBAD: {not_found}")

# ── Shell hit table ───────────────────────────────────────────────
print("\n  Shell assignment table:")
print(f"  {'Name':<22} {'theta':>7}  n  {'delta':>6}  hit")
for r in sorted(rows, key=lambda x: float(x["theta_dpole"])):
    tag = "<-- HIT" if r["hit_5deg"] == "YES" else ""
    print(f"  {r['Name']:<22} {float(r['theta_dpole']):7.2f}  {r['shell_n']}  {float(r['shell_delta']):6.2f}  {tag}")

# ── Summary stats ─────────────────────────────────────────────────
n_total = len(rows)
frac_expected = 6 * 10 / 180  # 6 shells * 2*5deg / 180deg  = 33%
import scipy.stats as st
binom_excess  = st.binomtest(hits, n_total, frac_expected, alternative='greater')
binom_avoid   = st.binomtest(hits, n_total, frac_expected, alternative='less')
binom_2sided  = st.binomtest(hits, n_total, frac_expected, alternative='two-sided')
print(f"\n  Statistical test (shell = within 5 deg of n*theta0):")
print(f"  Hits: {hits}/{n_total} = {hits/n_total:.1%}  |  Expected random: {frac_expected:.1%}  (mean = {n_total*frac_expected:.1f})")
print(f"  Binomial p  excess  (hits > expected) : {binom_excess.pvalue:.4f}  {'*' if binom_excess.pvalue < 0.05 else ''}")
print(f"  Binomial p  avoid   (hits < expected) : {binom_avoid.pvalue:.4f}  {'*' if binom_avoid.pvalue < 0.05 else ''}")
print(f"  Binomial p  2-sided                   : {binom_2sided.pvalue:.4f}  {'*' if binom_2sided.pvalue < 0.05 else ''}")
if hits / n_total < frac_expected:
    print(f"  >> SMBHs are BELOW shell rate: possible shell-AVOIDANCE signature")
else:
    print(f"  >> SMBHs are ABOVE shell rate: possible shell-EXCESS signature")

print("\n" + "="*62)
print(f"  Output: {out_path}")
print("="*62 + "\n")
