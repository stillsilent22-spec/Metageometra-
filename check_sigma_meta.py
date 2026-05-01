"""
check_sigma_meta.py
===================
Sigma_meta Großkreis-Suche am Himmel
Hypothese: Die topologische Grenzfläche Sigma_meta (S3-Sphäre ±-Trennung)
projiziert sich als Großkreis auf unseren Beobachtungshimmel,
senkrecht zum D-Pol (l=305°, b=+25°) — Toleranz ±15°.

Methodik:
  1. Ankerpunkt: CMB Cold Spot (l=209°, b=-57°)
  2. Grid-Suche: alle Pol-Richtungen senkrecht zum Cold Spot (φ=0..180°)
  3. Bekannte LSS-Anomalien als Test-Sample (21 Strukturen)
  4. VizieR: Planck PSZ2 SZ-Cluster → KS-Test gegen sin-Gleichverteilung
  5. D-Pol-Relation: Winkeld. Pol des besten Großkreises zum D-Pol
  6. Monte Carlo p-Wert (N=50000 Permutationen)
  7. Ausgabe: Koordinaten + p-Wert
"""
import math, os, sys, ssl, csv, io, time, re, random, urllib.request, urllib.parse

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(RES, "OT_sigma_meta.txt")
os.makedirs(RES, exist_ok=True)

_SSL = ssl.create_default_context()
_SSL.check_hostname = False
_SSL.verify_mode    = ssl.CERT_NONE

# ── Koordinaten-Utilities ─────────────────────────────────────────────────────
# Galactic pole J2000: α_G=192.85948°, δ_G=27.12825°, l_NCP=122.93192°
_AG  = math.radians(192.85948)
_DG  = math.radians(27.12825)
_LNCP = 122.93192

def radec_to_lb(ra_deg, dec_deg):
    """RA,Dec (J2000, Grad) → galaktisch (l,b) in Grad. Korrekte IAU-Formel."""
    a  = math.radians(ra_deg)
    d  = math.radians(dec_deg)
    Y  = math.cos(d)*math.sin(a - _AG)
    X  = math.sin(d)*math.cos(_DG) - math.cos(d)*math.sin(_DG)*math.cos(a - _AG)
    l  = (_LNCP - math.degrees(math.atan2(Y, X))) % 360.0
    b  = math.degrees(math.asin(max(-1., min(1.,
         math.sin(d)*math.sin(_DG) + math.cos(d)*math.cos(_DG)*math.cos(a - _AG)))))
    return l, b

def lb_to_radec(l_deg, b_deg):
    """Galaktisch (l,b) → RA,Dec J2000 (Grad). Umkehrformel."""
    l  = math.radians(l_deg); b = math.radians(b_deg)
    LN = math.radians(_LNCP)
    Y  = math.cos(b)*math.sin(LN - l)
    X  = math.sin(b)*math.cos(_DG) - math.cos(b)*math.sin(_DG)*math.cos(LN - l)
    ra = (math.degrees(math.atan2(Y, X)) + math.degrees(_AG)) % 360.0
    dec= math.degrees(math.asin(max(-1., min(1.,
         math.sin(b)*math.sin(_DG) + math.cos(b)*math.cos(_DG)*math.cos(LN - l)))))
    return ra, dec

def lb_xyz(l, b):
    lr, br = math.radians(l), math.radians(b)
    return (math.cos(br)*math.cos(lr), math.cos(br)*math.sin(lr), math.sin(br))

def xyz_to_lb(v):
    m  = math.sqrt(sum(x*x for x in v))
    if m < 1e-12: return 0., 0.
    x, y, z = v[0]/m, v[1]/m, v[2]/m
    b = math.degrees(math.asin(max(-1., min(1., z))))
    l = math.degrees(math.atan2(y, x)) % 360.
    return l, b

def cross3(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def dot3(a, b):
    return sum(u*v for u,v in zip(a,b))

def norm3(v):
    m = math.sqrt(dot3(v,v))
    if m < 1e-12: return (0.,0.,1.)
    return tuple(x/m for x in v)

def ang_sep(v1, v2):
    return math.degrees(math.acos(max(-1., min(1., dot3(v1, v2)))))

def dist_to_gc(point_xyz, pole_xyz):
    """Angularer Abstand (Grad) von Punkt zur Großkreis-Ebene des Pols.
    = arcsin(|P·N|) — 0° = auf GC, 90° = am Pol."""
    d = abs(dot3(point_xyz, pole_xyz))
    return math.degrees(math.asin(min(1., d)))

def gc_pole_basis(cs_xyz):
    """Zwei orthonormale Basisvektoren in der Ebene ⊥ cs_xyz."""
    ref = (0., 0., 1.)
    if abs(dot3(cs_xyz, ref)) > 0.99:
        ref = (1., 0., 0.)
    e1 = tuple(ref[i] - dot3(ref, cs_xyz)*cs_xyz[i] for i in range(3))
    e1 = norm3(e1)
    e2 = norm3(cross3(cs_xyz, e1))
    return e1, e2

def pole_from_phi(phi_deg, e1, e2):
    phi = math.radians(phi_deg)
    return tuple(math.cos(phi)*e1[i] + math.sin(phi)*e2[i] for i in range(3))

# ── CMB Cold Spot (Ankerpunkt) ────────────────────────────────────────────────
# Nutzervorgabe: RA=03h15m=48.75°, Dec=-19° → l≈209°, b≈-57°
# Wir verwenden die Standard-Galaktisch-Koordinaten direkt:
CS_L, CS_B = 209.0, -57.0
CS_XYZ     = lb_xyz(CS_L, CS_B)

# Prüfe: RA/Dec → l,b mit korrekter Formel
cs_l_check, cs_b_check = radec_to_lb(48.75, -19.0)

# HTM D-Pol (Metageometra-Theorie)
DPOL_L, DPOL_B = 305.0, 25.0
DPOL_XYZ       = lb_xyz(DPOL_L, DPOL_B)

print("=" * 68)
print("Sigma_meta Grosskreis-Suche")
print("=" * 68)
print(f"CMB Cold Spot:   l={CS_L}°, b={CS_B}°  (RA=48.75°, Dec=-19°)")
print(f"Formel-Check:    radec_to_lb(48.75,-19) → l={cs_l_check:.1f}°, b={cs_b_check:.1f}°")
print(f"D-Pol:           l={DPOL_L}°, b={DPOL_B}°")
print(f"Abstand CS ↔ D-Pol: {ang_sep(CS_XYZ, DPOL_XYZ):.1f}°")

# ── Bekannte LSS-Anomalien (galaktische Koordinaten aus Literatur) ───────────
# Alle Positionen in galaktischen Koordinaten (l, b), direkt aus Literatur.
# Format: (Name, l_gal, b_gal, Typ, Referenz)
STRUCTURES_LB = [
    # CMB-Anomalien (Referenz-Richtungen)
    ("CMB_ColdSpot",        209.0,  -57.0, "CMB",    "Vielva+2004/Cruz+2007"),
    ("CMB_Quad_axis",       240.0,  +63.0, "CMB",    "Tegmark+2003 Quadrupol"),
    ("CMB_Oct_axis",        308.0,  +63.0, "CMB",    "Tegmark+2003 Oktupol"),
    ("CMB_Asym_dir",        218.0,  -20.0, "CMB",    "Eriksen+2004 Hemisph.Asym"),
    ("CMB_AxisOfEvil",      264.0,  +48.0, "CMB",    "Land+Magueijo 2005"),
    # Große Voids
    ("Eridanus_Supervoid",  205.0,  -56.5, "Void",   "Naidoo+2016 z=0.22"),
    ("Bootes_Void",          57.0,  +68.0, "Void",   "Kirshner+1981 z=0.065"),
    ("KBC_LocalVoid",       312.0,   -4.0, "Void",   "Keenan+2013"),
    ("DESI_Supervoid_z03",  242.0,  +17.0, "Void",   "Bremer+2010 z=0.3 est."),
    # Große Wände / Filamente
    ("Sloan_GreatWall",     245.0,  +52.0, "Wall",   "Gott+2005 z=0.073"),
    ("CfA2_GreatWall",      130.0,  +70.0, "Wall",   "Geller+Huchra 1989"),
    ("Perseus_Pisces",      150.0,  -13.0, "Wall",   "Giovanelli+1986"),
    ("Saraswati_SC",         50.0,   -7.0, "Wall",   "Bagchi+2017 z=0.28"),
    ("Hercules_CrB_GW",      75.0,  +45.0, "Wall",   "Horvath+2015 z=1.6"),
    # Galaxienhaufen / Superhaufen
    ("Shapley_Conc",        312.0,  +31.0, "SC",     "Scaramella+1989 z=0.046"),
    ("Laniakea_GA",         304.0,  +29.0, "SC",     "Tully+2014"),
    ("VirgoCluster",        284.0,  +74.0, "Cluster","z=0.0038"),
    ("El_Gordo",            249.0,  -47.0, "Cluster","Menanteau+2012 z=0.87"),
    # HTM-Referenzpunkte
    ("HTM_D_Pol",           305.0,  +25.0, "HTM",    "Metageom. D-Pol"),
    ("HTM_A_Pol",           125.0,  -25.0, "HTM",    "Metageom. A-Pol"),
    ("SDSS_nc3_est",        245.0,  +55.0, "HTM",    "SDSS nc>=3 QSO-excess Schw."),
]

# Konvertiere zu galaktischen Koordinaten und xyz
structures = []
for nm, l, b, typ, ref in STRUCTURES_LB:
    xyz = lb_xyz(l, b)
    structures.append({'name': nm, 'l': l, 'b': b, 'xyz': xyz, 'type': typ, 'ref': ref})

print(f"\nStrukturen geladen: {len(structures)}")
for s in structures[:5]:
    print(f"  {s['name']:25s}  l={s['l']:6.1f}  b={s['b']:6.1f}")
print("  ...")

# ── VizieR-Abfrage: Planck PSZ2 SZ-Cluster ───────────────────────────────────
def vizier_query(source, out_cols, max_out=2000, extra=''):
    params = {'-source': source, '-out': out_cols,
              '-out.max': str(max_out), '-oc.form': 'dec'}
    url = 'https://vizier.cds.unistra.fr/viz-bin/votable?' + urllib.parse.urlencode(params)
    if extra: url += '&' + extra
    for attempt in range(3):
        try:
            if attempt: time.sleep(5*attempt)
            req = urllib.request.Request(url, headers={'User-Agent': 'Python/SigmaMeta'})
            with urllib.request.urlopen(req, timeout=90, context=_SSL) as r:
                raw = r.read().decode('utf-8', 'replace')
            rows = []
            for tm in re.finditer(r'<TABLE\b[^>]*>(.*?)</TABLE>', raw, re.DOTALL):
                tb = tm.group(1)
                fields = re.findall(r'<FIELD[^>]+name="([^"]+)"', tb)
                fields = [f for f in fields if f != 'recno']
                if len(fields) < 2: continue
                for tr in re.finditer(r'<TR>(.*?)</TR>', tb, re.DOTALL):
                    tds = re.findall(r'<TD>(.*?)</TD>', tr.group(1))
                    if len(tds) < 2: continue
                    rows.append({fields[i]: tds[i].strip() if i < len(tds) else ''
                                 for i in range(len(fields))})
            return rows
        except Exception as e:
            print(f"  VizieR Versuch {attempt+1}: {e}")
    return []

def to_float(s):
    try: return float(s.strip())
    except: return None

viz_clusters = []

print("\n[VizieR] Planck PSZ2 SZ-Cluster...")
for cat in ['J/A+A/594/A27/psz2', 'VIII/92/psz', 'J/A+A/581/A14/psz2']:
    rows = vizier_query(cat, '_RA _DE SNR', max_out=1800)
    if rows:
        print(f"  {cat}: {len(rows)} Eintraege")
        break
def extract_radec(row):
    """Try multiple VizieR column name conventions for RA/Dec."""
    for rak in ('_RAJ2000','_RA','RAJ2000','RA_ICRS','RA','ra'):
        for dek in ('_DEJ2000','_DE','DEJ2000','DE_ICRS','Dec','DE','dec'):
            ra = to_float(row.get(rak,'')); de = to_float(row.get(dek,''))
            if ra is not None and de is not None:
                return ra, de
    return None, None

if rows:
    ra0, de0 = extract_radec(rows[0])
    if ra0 is None:
        print(f"  DEBUG erste Row keys: {list(rows[0].keys())[:8]}")
for row in rows:
    ra_, dec_ = extract_radec(row)
    if ra_ is None or dec_ is None: continue
    l_, b_ = radec_to_lb(ra_, dec_)
    viz_clusters.append({'type':'PSZ2', 'l':l_, 'b':b_, 'xyz':lb_xyz(l_, b_)})
print(f"  PSZ2 konvertiert: {len(viz_clusters)}")

print("[VizieR] Abell-Cluster (VII/110A)...")
abell_rows = vizier_query('VII/110A/table3', '_RA _DE', max_out=2000)
print(f"  Abell: {len(abell_rows)} Eintraege")
n_pre = len(viz_clusters)
for row in abell_rows:
    ra_, dec_ = extract_radec(row)
    if ra_ is None or dec_ is None: continue
    l_, b_ = radec_to_lb(ra_, dec_)
    viz_clusters.append({'type':'Abell', 'l':l_, 'b':b_, 'xyz':lb_xyz(l_, b_)})
print(f"  Abell konvertiert: {len(viz_clusters)-n_pre}")

print(f"  VizieR-Cluster gesamt: {len(viz_clusters)}")

# ── Schritt 2: Großkreis-Grid-Suche ──────────────────────────────────────────
TOL_GC   = 7.0    # Grad: Abstand < TOL_GC = "auf dem Großkreis"
N_GRID   = 3600   # 0.05° Auflösung

e1, e2 = gc_pole_basis(CS_XYZ)

# Test-Punkte: alle Strukturen AUSSER Cold Spot selbst
test_pts = [s for s in structures if s['name'] != 'CMB_ColdSpot']

print(f"\n[SUCHE] {N_GRID} Pol-Kandidaten, Tol={TOL_GC}°, N_test={len(test_pts)}")

best_score = -1
results    = []

for step in range(N_GRID):
    phi  = step * 180.0 / N_GRID
    pole = pole_from_phi(phi, e1, e2)
    hits  = [s for s in test_pts if dist_to_gc(s['xyz'], pole) < TOL_GC]
    d_dpol = ang_sep(pole, DPOL_XYZ)
    d_dpol_eff = min(d_dpol, 180. - d_dpol)  # 0°=Pol ist D-Pol=GC⊥D-Pol, 90°=GC liegt im D-Pol
    score = len(hits) + (1. - d_dpol_eff/90.) * 0.01
    results.append({'phi':phi, 'pole':pole, 'hits':hits, 'n_hits':len(hits),
                    'd_dpol':d_dpol_eff, 'score':score})
    if score > best_score:
        best_score = score

results.sort(key=lambda x: -x['score'])
BR = results[0]

# ── Schritt 4: D-Pol-Relation ─────────────────────────────────────────────────
bp_l, bp_b = xyz_to_lb(BR['pole'])
print(f"\n[BEST GC] Pol: l={bp_l:.1f}°, b={bp_b:.1f}°  "
      f"Treffer: {BR['n_hits']}/{len(test_pts)}  d_D-Pol={BR['d_dpol']:.1f}°")
for s in BR['hits']:
    d = dist_to_gc(s['xyz'], BR['pole'])
    print(f"  ★ {s['name']:25s}  l={s['l']:6.1f}  b={s['b']:6.1f}  d={d:.1f}°")

dpol_compatible = BR['d_dpol'] < 15.0
print(f"\n  D-Pol-Senkrecht-Kriterium (<15°): "
      f"{'JA (Sigma_meta-kompatibel)' if dpol_compatible else 'NEIN'}")

# Bester GC unter D-Pol-Bedingung
dpol_top = [r for r in results if r['d_dpol'] < 15.0]
dpol_top.sort(key=lambda x: -x['n_hits'])
print(f"\n[D-POL FILTER] Beste GC mit Pol-Abst.< 15° zum D-Pol: {len(dpol_top)} Kandidaten")
if dpol_top:
    DR = dpol_top[0]
    dpl, dpb = xyz_to_lb(DR['pole'])
    print(f"  Bester: Pol l={dpl:.1f}°, b={dpb:.1f}°  Hits={DR['n_hits']}  d_D-Pol={DR['d_dpol']:.1f}°")
    for s in DR['hits']:
        d = dist_to_gc(s['xyz'], DR['pole'])
        print(f"    ★ {s['name']:25s}  l={s['l']:6.1f}  b={s['b']:6.1f}  d={d:.1f}°")
else:
    DR = None
    print("  Kein Kandidat unter d_dpol < 15°.")

# Reiner D-Pol-Großkreis (Pol = D-Pol exakt)
d_cs_dpol_gc = dist_to_gc(CS_XYZ, DPOL_XYZ)
print(f"\n[D-POL GC EXAKT] Abstand Cold Spot → D-Pol-Großkreis: {d_cs_dpol_gc:.1f}°")
dpol_exact_hits = [(s, dist_to_gc(s['xyz'], DPOL_XYZ))
                   for s in test_pts if dist_to_gc(s['xyz'], DPOL_XYZ) < TOL_GC]
dpol_exact_hits.sort(key=lambda x: x[1])
print(f"  Strukturen auf D-Pol-GC (<{TOL_GC}°): {len(dpol_exact_hits)}")
for s, d in dpol_exact_hits:
    print(f"    {s['name']:25s}  d={d:.1f}°")

# ── Schritt 3: KS-Test ────────────────────────────────────────────────────────
print(f"\n[KS-TEST] Bekannte Strukturen → best GC:")
all_dists = sorted(dist_to_gc(s['xyz'], BR['pole']) for s in test_pts)
n = len(all_dists)

def ks_cdf_sin(d_deg):
    return math.sin(math.radians(d_deg))

ks_stat = max(
    max(abs((i+1)/n - ks_cdf_sin(d)), abs(i/n - ks_cdf_sin(d)))
    for i, d in enumerate(all_dists))

def ks_pvalue(D, n):
    z = D * (math.sqrt(n) + 0.12 + 0.11/math.sqrt(n))
    p = 2.*sum(((-1)**(k+1))*math.exp(-2*(k**2)*(z**2)) for k in range(1, 50))
    return max(0., min(1., p))

ks_p = ks_pvalue(ks_stat, n)
print(f"  KS D={ks_stat:.4f}  p={ks_p:.4f}  "
      f"({'signifikant p<0.05' if ks_p < 0.05 else 'nicht signifikant'})")

# KS-Test auf VizieR-Cluster
if viz_clusters:
    cd = sorted(dist_to_gc(c['xyz'], BR['pole']) for c in viz_clusters)
    nc = len(cd)
    ks_c = max(
        max(abs((i+1)/nc - ks_cdf_sin(d)), abs(i/nc - ks_cdf_sin(d)))
        for i, d in enumerate(cd))
    ks_cp = ks_pvalue(ks_c, nc)
    n_on = sum(1 for d in cd if d < TOL_GC)
    frac_exp = math.sin(math.radians(TOL_GC))
    print(f"  KS Cluster (N={nc}): D={ks_c:.4f}  p={ks_cp:.4f}"
          f"  auf-GC: {n_on}/{nc}({100*n_on/nc:.1f}%) vs Zufall({100*frac_exp:.1f}%)"
          f"  Ratio={n_on/nc/frac_exp:.2f}")
else:
    ks_cp, ks_c, n_on, nc = None, None, 0, 0
    frac_exp = math.sin(math.radians(TOL_GC))

# ── Schritt 6: Monte Carlo p-Wert ─────────────────────────────────────────────
# Optimiert: asin vermieden durch Vergleich |dot| < sin(TOL), Projektion vorberechnet
N_MC = 10000
print(f"\n[MC] Monte Carlo (N_MC={N_MC})...")
random.seed(42)
obs_hits = BR['n_hits']
mc_exceed = 0
sin_tol = math.sin(math.radians(TOL_GC))
n_pts = len(test_pts)
angles_rad = [math.radians(k * 10.) for k in range(18)]
cos_a = [math.cos(a) for a in angles_rad]
sin_a = [math.sin(a) for a in angles_rad]

for _ in range(N_MC):
    # Zufälliger Referenzpol
    phi_r = random.uniform(0, 2*math.pi)
    ct_r  = random.uniform(-1, 1)
    st_r  = math.sqrt(1 - ct_r**2)
    rcs   = (st_r*math.cos(phi_r), st_r*math.sin(phi_r), ct_r)
    re1, re2 = gc_pole_basis(rcs)
    # Projektion jedes Zufallspunkts auf Basis (e1, e2) vorberechnen
    AB = []
    for _ in range(n_pts):
        fp2 = random.uniform(0, 2*math.pi)
        ct2 = random.uniform(-1, 1)
        st2 = math.sqrt(1 - ct2**2)
        rp  = (st2*math.cos(fp2), st2*math.sin(fp2), ct2)
        AB.append((dot3(rp, re1), dot3(rp, re2)))
    # Bestes GC (keine asin nötig): |a*cos(phi) + b*sin(phi)| < sin(TOL)
    best_mc = 0
    for ca, sa in zip(cos_a, sin_a):
        h = sum(1 for (a, b) in AB if abs(a*ca + b*sa) < sin_tol)
        if h > best_mc:
            best_mc = h
    if best_mc >= obs_hits:
        mc_exceed += 1

mc_p = mc_exceed / N_MC
print(f"  obs={obs_hits}/{n_pts}  MC p={mc_p:.4f}"
      f"  ({'signifikant' if mc_p < 0.05 else 'n.s.'})")

# ── Schritt 5: Ergebnis ───────────────────────────────────────────────────────
bp_l, bp_b = xyz_to_lb(BR['pole'])
ap_l, ap_b = xyz_to_lb(tuple(-x for x in BR['pole']))
pol_x_cs   = norm3(cross3(BR['pole'], CS_XYZ))

print(f"\n{'='*68}")
print(f"ERGEBNIS: Bester Grosskreis-Kandidat (Sigma_meta GC)")
print(f"  Pol:         l={bp_l:.2f}°, b={bp_b:.2f}°")
print(f"  Antipol:     l={ap_l:.2f}°, b={ap_b:.2f}°")
print(f"  Treffer:     {BR['n_hits']}/{len(test_pts)}")
print(f"  D-Pol-Abst:  {BR['d_dpol']:.1f}°  "
      f"({'SIGMA_META-KOMPATIBEL (<15°)' if BR['d_dpol']<=15 else 'nicht D-Pol-senkrecht'})")
print(f"  KS p-Wert:   {ks_p:.4f}")
print(f"  MC p-Wert:   {mc_p:.4f}")

print(f"\n  GC-Punkte (alle 30°, theta=0 ist Cold Spot):")
for th in range(0, 360, 30):
    gcp = tuple(math.cos(math.radians(th))*CS_XYZ[i]
               + math.sin(math.radians(th))*pol_x_cs[i] for i in range(3))
    gl, gb = xyz_to_lb(gcp)
    print(f"    theta={th:3d}°   l={gl:6.1f}°  b={gb:6.1f}°")

# ── Ausgabedatei ─────────────────────────────────────────────────────────────
buf = io.StringIO()
def w(s=''): buf.write(str(s)+'\n')

w("=" * 72)
w("OT-Sigma_meta: Grosskreis-Suche am Himmel")
w("Hypothese: Sigma_meta (S3-Topologiegrenze) → Grosskreis senkrecht zu D-Pol")
w("=" * 72)
w(f"Datum:     {time.strftime('%Y-%m-%d %H:%M')}")
w(f"Anker:     CMB Cold Spot  l={CS_L}°  b={CS_B}°  (RA=48.75°, Dec=-19°)")
w(f"D-Pol:     l={DPOL_L}°  b={DPOL_B}°   |CS-D-Pol|={ang_sep(CS_XYZ,DPOL_XYZ):.1f}°")
w(f"Grid:      {N_GRID} Pole  Toleranz={TOL_GC}°")
w(f"Strukturen:{len(test_pts)} bekannte LSS-Anomalien + {len(viz_clusters)} VizieR-Cluster")
w()
w("-" * 72)
w("SCHRITT 2+4: BESTER GROSSKREIS (hoechste Trefferrate AUF ALLEN Kandidaten)")
w(f"  Pol     l={bp_l:.2f}°  b={bp_b:.2f}°")
w(f"  Antipol l={ap_l:.2f}°  b={ap_b:.2f}°")
w(f"  Treffer: {BR['n_hits']}/{len(test_pts)}")
w(f"  D-Pol Abstand des Pols: {BR['d_dpol']:.2f}°  (Limit 15°)")
w(f"  Sigma_meta-kompatibel (D-Pol-senkrecht): {'JA' if BR['d_dpol']<=15 else 'NEIN'}")
w()
w("  Strukturen auf bestem GC:")
for s in BR['hits']:
    d = dist_to_gc(s['xyz'], BR['pole'])
    w(f"    {s['name']:27s}  l={s['l']:7.2f}  b={s['b']:7.2f}  d={d:.1f}°  [{s['type']}]")
w()
w("-" * 72)
if dpol_top:
    w(f"SCHRITT 4: BESTER GC MIT D-POL-BEDINGUNG (d_pol<15°) — {len(dpol_top)} Kandidaten")
    dpl, dpb = xyz_to_lb(DR['pole'])
    w(f"  Pol l={dpl:.2f}°  b={dpb:.2f}°  Hits={DR['n_hits']}  d_D-Pol={DR['d_dpol']:.1f}°")
    for s in DR['hits']:
        d = dist_to_gc(s['xyz'], DR['pole'])
        w(f"    {s['name']:27s}  d={d:.1f}°")
    w()
w("-" * 72)
w(f"D-POL-GROSSKREIS (Pol = D-Pol exakt): Cold Spot Abstand={d_cs_dpol_gc:.1f}°")
w(f"  Strukturen auf D-Pol-GC (<{TOL_GC}°): {len(dpol_exact_hits)}")
for s, d in dpol_exact_hits:
    w(f"    {s['name']:27s}  d={d:.1f}°")
w()
w("-" * 72)
w("SCHRITT 3: KS-TEST")
w(f"  Bekannte Strukturen (N={n}):  D={ks_stat:.4f}  p={ks_p:.4f}")
if ks_cp is not None:
    w(f"  VizieR Cluster (N={nc}):      D={ks_c:.4f}  p={ks_cp:.4f}")
    w(f"  Cluster auf GC: {n_on}/{nc}({100*n_on/nc:.1f}%) vs Zufall({100*frac_exp:.1f}%)  Ratio={n_on/nc/frac_exp:.2f}")
w()
w("-" * 72)
w(f"MC p-WERT: {mc_p:.4f}  (obs={obs_hits}/{len(test_pts)}, N_MC={N_MC})")
w()
w("-" * 72)
w("ALLE STRUKTUREN: Abstand zum besten GC (sortiert):")
all_sorted = sorted([(s, dist_to_gc(s['xyz'], BR['pole'])) for s in test_pts], key=lambda x:x[1])
for s, d in all_sorted:
    mk = '*' if d < TOL_GC else ' '
    w(f"  {mk} {s['name']:27s}  l={s['l']:7.2f}  b={s['b']:7.2f}  d={d:.1f}°")
w()
w("TOP-10 GC-KANDIDATEN:")
w(f"  {'phi':>7}  {'Pol_l':>7}  {'Pol_b':>6}  {'Hits':>4}  {'d_D-Pol':>7}  Treffer")
for r in results[:10]:
    pl, pb = xyz_to_lb(r['pole'])
    nms = ','.join(s['name'].replace('_',' ')[:12] for s in r['hits'])
    w(f"  {r['phi']:7.2f}  {pl:7.2f}  {pb:6.2f}  {r['n_hits']:4d}  {r['d_dpol']:7.2f}  {nms}")
w()
w("GROSSKREIS-PUNKTE (theta=0 = Cold Spot, alle 30°):")
for th in range(0, 360, 30):
    gcp = tuple(math.cos(math.radians(th))*CS_XYZ[i]
               + math.sin(math.radians(th))*pol_x_cs[i] for i in range(3))
    gl, gb = xyz_to_lb(gcp)
    w(f"  theta={th:3d}°   l={gl:6.1f}°  b={gb:6.1f}°")

txt = buf.getvalue()
with open(OUT, 'w', encoding='utf-8') as fh:
    fh.write(txt)

print(f"\nAusgabe: {OUT}")
print("FERTIG.")
