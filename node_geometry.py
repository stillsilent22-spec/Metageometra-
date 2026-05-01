"""
HTM Knotenpunkt Geometrie-Analyse
===================================
Geht ALLE möglichen Kandidaten im kombinierten Galaxien-Katalog durch.
Für jeden Kandidaten nahe einem HTM-Knoten:
  1. Abstand zur nächsten Thales-Sphäre (|R_Thales - R_n| / R_n)
  2. Nächste bekannte SMBHs auf Knoten (nc >= 1) suchen
  3. Mit mind. 3-4 SMBHs: Kollinearitäts-Test (Großkreis auf Himmelssphäre)
  4. 3D-Filament-Test (wenn Distanzen bekannt)
  5. Geometrische Interpretation

HTM-Knoten = Punkt, wo Schalen-Ringe MEHRERER Pole sich schneiden.
Crossing-Density nc = Anzahl der Paare (Pol_i, n_i) / (Pol_j, n_j) mit i≠j,
die beide innerhalb TOL_ANG vom Punkt liegen.
"""
import math, csv, os, sys, itertools, statistics

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE     = os.path.dirname(os.path.abspath(__file__))
RES      = os.path.join(BASE, "results")
OUT_TXT  = os.path.join(RES, "OT_NG_v2.txt")
SMBH_CSV = os.path.join(RES, "catalogs", "smbh_extended.csv")
GAL_CSV  = os.path.join(RES, "catalogs", "vizier_dpol_combined.csv")

# ── Konstanten ──────────────────────────────────────────────────────────────────
THETA0      = 58.65    # HTM Fraktal-Winkel [°]
R_S         = 4228.3   # HTM Vakuum-Sphäre [Mpc]
H0          = 73.0     # km/s/Mpc
TOL_ANG     = 8.0      # Toleranz Schalen-Winkel [°]
TOL_GC      = 7.0      # Max. RMS für Großkreis-Alignment [°]
TOL_3D      = 20.0     # Max. RMS für 3D-Filament [Mpc]
NEIGHBOR_DEG = 25.0    # Suchhalbmesser für Nachbar-SMBHs [°]
MAX_COMBO   = 6        # Max. Kombinations-Größe für Kollinearitäts-Test

# HTM-Pole (galaktische Koordinaten)
POLES_DEF = [
    (305.0,  25.0, 'D'),    # Haupt-Kollisionspol
    (125.0, -25.0, 'A'),    # Antipodal
    ( 35.0,  25.0, 'S1'),   # Sekundär 90° auf D-A-Großkreis
    (215.0, -25.0, 'S2'),
    (215.0,  25.0, 'S3'),
    ( 35.0, -25.0, 'S4'),
]
SHELLS_N   = [1, 2, 3]
SHELLS_DEG = [THETA0 * n for n in SHELLS_N]

# ── Geometrie-Helfer ─────────────────────────────────────────────────────────────
def lb_xyz(l, b):
    l, b = math.radians(l), math.radians(b)
    return math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)

def xyz_lb(x, y, z):
    r = math.sqrt(x*x + y*y + z*z)
    if r < 1e-15: return 0.0, 0.0
    b = math.degrees(math.asin(max(-1.0, min(1.0, z/r))))
    l = math.degrees(math.atan2(y, x)) % 360.0
    return l, b

def ang_sep(v1, v2):
    d = max(-1.0, min(1.0, v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]))
    return math.degrees(math.acos(d))

def cross3(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def dot3(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm3(v):
    m = math.sqrt(dot3(v, v))
    return (v[0]/m, v[1]/m, v[2]/m) if m > 1e-15 else (0.0, 0.0, 0.0)

def ang_lb(l1, b1, l2, b2):
    return ang_sep(lb_xyz(l1, b1), lb_xyz(l2, b2))

def radec_to_lb(ra_deg, dec_deg):
    ra = math.radians(ra_deg); dc = math.radians(dec_deg)
    RA_NGP = math.radians(192.85948); DEC_NGP = math.radians(27.12825)
    b = math.asin(math.sin(dc)*math.sin(DEC_NGP) +
                  math.cos(dc)*math.cos(DEC_NGP)*math.cos(ra - RA_NGP))
    y = math.cos(dc)*math.sin(ra - RA_NGP)
    x = (math.cos(dc)*math.sin(DEC_NGP)*math.cos(ra - RA_NGP)
         - math.sin(dc)*math.cos(DEC_NGP))
    l = (math.degrees(math.atan2(y, x)) + 122.93192) % 360.0
    return l, math.degrees(b)

# Vorberechnete Pol-Vektoren
POLES_XYZ = [(lb_xyz(l, b), nm) for l, b, nm in POLES_DEF]

def get_hits(l, b, tol=TOL_ANG):
    """Alle (Pol, n, Δ°, θ_gesamt) wo dieser Punkt nahe einer HTM-Schale liegt."""
    p = lb_xyz(l, b)
    hits = []
    for pv, pn in POLES_XYZ:
        theta = ang_sep(p, pv)
        for n, sdeg in zip(SHELLS_N, SHELLS_DEG):
            delta = abs(theta - sdeg)
            if delta < tol:
                hits.append((pn, n, round(delta, 2)))
    return hits

def crossing_density(l, b, tol=TOL_ANG):
    """Anzahl Pol-Schalen-Kreuzungen + Details."""
    hits = get_hits(l, b, tol)
    nc = sum(1 for i, h1 in enumerate(hits)
             for h2 in hits[i+1:] if h1[0] != h2[0])
    return nc, hits

def best_thales_shell(r_thales):
    """Nächste HTM Thales-Sphäre. Returns (n, R_n, frac_abstand)."""
    best = None
    for n in range(1, 15):
        Rn = R_S * n * THETA0 / 360.0
        frac = abs(r_thales - Rn) / Rn
        if best is None or frac < best[2]:
            best = (n, round(Rn, 1), round(frac, 4))
    return best

# ── Großkreis-Anpassung ──────────────────────────────────────────────────────────
def great_circle_fit(lb_list):
    """
    Passt Großkreis an N Himmelspunkte (l,b) an.
    Returns: (normal_lb, residuals_deg, rms_deg)
    Methode: Mittel aller Paar-Normalvektoren (robust für N≥3)
    """
    vecs = [lb_xyz(l, b) for l, b in lb_list]
    pairs = list(itertools.combinations(vecs, 2))
    if not pairs:
        return None, [], 999.9
    normals = []
    for va, vb in pairs:
        n_ab = cross3(va, vb)
        m = math.sqrt(dot3(n_ab, n_ab))
        if m > 1e-10:
            normals.append((n_ab[0]/m, n_ab[1]/m, n_ab[2]/m))
    if not normals:
        return None, [], 999.9
    # Alle Normalen in dieselbe Hemisphäre
    ref = normals[0]
    aligned = [(n if dot3(n, ref) >= 0 else (-n[0], -n[1], -n[2])) for n in normals]
    avg_n = norm3((
        sum(n[0] for n in aligned) / len(aligned),
        sum(n[1] for n in aligned) / len(aligned),
        sum(n[2] for n in aligned) / len(aligned)
    ))
    residuals = []
    for v in vecs:
        sin_r = dot3(v, avg_n)
        residuals.append(abs(math.degrees(math.asin(max(-1.0, min(1.0, sin_r))))))
    rms = math.sqrt(sum(r*r for r in residuals) / len(residuals))
    nl, nb = xyz_lb(*avg_n)
    return (nl, nb), residuals, rms

# ── 3D Geradenfit ──────────────────────────────────────────────────────────────
def line_fit_3d(points_3d):
    """
    Passt Gerade an 3D-Punkte (x,y,z in Mpc) an.
    Returns: (direction_lb, residuals_mpc, rms_mpc)
    """
    n = len(points_3d)
    if n < 2:
        return None, [], 999.9
    cx = sum(p[0] for p in points_3d) / n
    cy = sum(p[1] for p in points_3d) / n
    cz = sum(p[2] for p in points_3d) / n
    ctr = [(p[0]-cx, p[1]-cy, p[2]-cz) for p in points_3d]
    # Anfangsrichtung: erstes zu zweitem Punkt
    d = norm3((ctr[1][0]-ctr[0][0], ctr[1][1]-ctr[0][1], ctr[1][2]-ctr[0][2]))
    for _ in range(60):
        s = [0.0, 0.0, 0.0]
        for p in ctr:
            sgn = 1.0 if dot3(p, d) >= 0 else -1.0
            s[0] += sgn*p[0]; s[1] += sgn*p[1]; s[2] += sgn*p[2]
        d_new = norm3(tuple(s))
        if math.sqrt(sum((a-b)**2 for a, b in zip(d, d_new))) < 1e-9:
            break
        d = d_new
    residuals = []
    for p in ctr:
        proj = dot3(p, d)
        perp = (p[0]-proj*d[0], p[1]-proj*d[1], p[2]-proj*d[2])
        residuals.append(math.sqrt(dot3(perp, perp)))
    rms = math.sqrt(sum(r*r for r in residuals) / len(residuals))
    dl, db = xyz_lb(*d)
    return (dl, db), residuals, round(rms, 2)

# ── Daten laden ──────────────────────────────────────────────────────────────────
print("Lade SMBH-Katalog...", flush=True)
smbhs = []
with open(SMBH_CSV, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try:
            lm  = float(r['logMbh'])
            ra  = float(r['RA_deg'])
            dc  = float(r['Dec_deg'])
            l, b = radec_to_lb(ra, dc)
            d_mpc = None
            try:
                v = r.get('dist_mpc', '').strip()
                if v: d_mpc = float(v)
            except ValueError:
                pass
            theta_d = ang_lb(l, b, 305.0, 25.0)
            r_thales = None
            if d_mpc and d_mpc > 0:
                ct = math.cos(math.radians(theta_d))
                if ct > 0.05:
                    r_thales = d_mpc / ct
            nc, hits = crossing_density(l, b)
            ths = best_thales_shell(r_thales) if r_thales else None
            xyz3 = None
            if d_mpc and d_mpc > 0:
                ux, uy, uz = lb_xyz(l, b)
                xyz3 = (ux*d_mpc, uy*d_mpc, uz*d_mpc)
            smbhs.append({
                'name': r['Name'].strip(), 'l': l, 'b': b,
                'logM': lm, 'd_mpc': d_mpc, 'xyz3': xyz3,
                'theta_d': theta_d, 'r_thales': r_thales,
                'nc': nc, 'hits': hits, 'thales': ths
            })
        except Exception:
            continue

print(f"  {len(smbhs)} SMBHs, davon {sum(1 for s in smbhs if s['d_mpc'])} mit Distanz", flush=True)

print("Lade Galaxien-Katalog...", flush=True)
galaxies = []
if os.path.exists(GAL_CSV):
    with open(GAL_CSV, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            try:
                l  = float(r['l']); b = float(r['b'])
                d  = float(r['d_mpc']); rt = float(r['R_Thales'])
                nc, hits = crossing_density(l, b)
                ths = best_thales_shell(rt)
                galaxies.append({
                    'name': r['name'], 'src': r.get('source', '?'),
                    'l': l, 'b': b, 'd_mpc': d, 'R_Thales': rt,
                    'nc': nc, 'hits': hits, 'thales': ths,
                    'theta_d': ang_lb(l, b, 305.0, 25.0)
                })
            except Exception:
                continue
print(f"  {len(galaxies)} Galaxien geladen", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
lines = []
lines.append("=" * 74)
lines.append("HTM KNOTENPUNKT GEOMETRIE-ANALYSE")
lines.append(f"THETA0={THETA0}°  r_s={R_S} Mpc  Pole={[p[2] for p in POLES_DEF]}")
lines.append(f"Tol: ±{TOL_ANG}° (Schalen) | GK-RMS≤{TOL_GC}° | 3D-RMS≤{TOL_3D} Mpc | Nachbar≤{NEIGHBOR_DEG}°")
lines.append("=" * 74)

# ── ABSCHNITT 1: SMBHs an Knoten ─────────────────────────────────────────────
lines.append("\n══ 1. SMBHs AN HTM-KNOTEN ══")
node_smbhs = sorted([s for s in smbhs if s['nc'] >= 1], key=lambda x: (-x['nc'], -x['logM']))
all_smbhs_nc0 = sorted([s for s in smbhs if s['nc'] == 0], key=lambda x: -x['logM'])

lines.append(f"\nAlle {len(smbhs)} SMBHs nach Kreuzungsdichte:")
lines.append(f"  nc=0 (auf keiner Schale):   {sum(1 for s in smbhs if s['nc']==0):3d} SMBHs")
lines.append(f"  nc=1 (1 Kreuzung):           {sum(1 for s in smbhs if s['nc']==1):3d} SMBHs")
lines.append(f"  nc=2 (2 Kreuzungen = Knoten):{sum(1 for s in smbhs if s['nc']>=2):3d} SMBHs")

lines.append(f"\nSMBHs mit nc ≥ 1 (auf mind. einer Schalen-Kreuzung):")
lines.append(f"{'Name':<24} {'logM':>6}  {'nc':>3}  {'θ_D':>6}  {'d[Mpc]':>8}  Schalen-Hits")
lines.append("-" * 78)
for s in node_smbhs:
    hits_str = '  '.join(f"{pn}(n={n},δ={d}°)" for pn, n, d in s['hits'])
    d_str = f"{s['d_mpc']:.1f}" if s['d_mpc'] else "  ?"
    ths_str = ""
    if s['thales']:
        ths_str = f"  R_Th={s['r_thales']:.0f}→n={s['thales'][0]}±{100*s['thales'][2]:.0f}%"
    lines.append(f"  {s['name']:<22} {s['logM']:>6.2f}  {s['nc']:>3}  {s['theta_d']:>5.1f}°  "
                 f"{d_str:>8}  {hits_str}{ths_str}")

# ── ABSCHNITT 2: Großkreis-Kollinearität ──────────────────────────────────────
lines.append("\n\n══ 2. GROßKREIS-KOLLINEARITÄT (alle SMBHs, nc ≥ 1) ══")
cands_gc = node_smbhs[:]        # alle mit nc ≥ 1

gc_found = []
n_tested = 0
MIN_SEP = 2.0  # min. Winkeltrennung zwischen je zwei Objekten [°]
for sz in range(3, min(MAX_COMBO + 1, len(cands_gc) + 1)):
    combos = list(itertools.combinations(cands_gc, sz))
    if n_tested + len(combos) > 50000:
        break
    n_tested += len(combos)
    for combo in combos:
        # Triviale Alignments filtern: mind. MIN_SEP zwischen je zwei Objekten
        too_close = any(ang_lb(combo[i]['l'], combo[i]['b'],
                                combo[j]['l'], combo[j]['b']) < MIN_SEP
                        for i in range(len(combo))
                        for j in range(i+1, len(combo)))
        if too_close:
            continue
        lbs = [(s['l'], s['b']) for s in combo]
        n_circ, residuals, rms = great_circle_fit(lbs)
        if n_circ and rms < TOL_GC:
            # Winkelseparationen zwischen benachbarten Objekten
            seps = [round(ang_lb(combo[i]['l'], combo[i]['b'],
                                 combo[i+1]['l'], combo[i+1]['b']), 1)
                    for i in range(len(combo) - 1)]
            # Passiert der GK-Normalvektor in der Nähe eines Pols?
            pole_dists = {pn: ang_lb(n_circ[0], n_circ[1], pl, pb)
                          for pl, pb, pn in POLES_DEF}
            near_poles = [(pn, round(d, 1)) for pn, d in pole_dists.items()
                          if d < 20 or d > 160]
            # Sind Separationen Vielfache von THETA0?
            theta_matches = [s for s in seps if abs(s % THETA0) < 5 or abs(s % THETA0 - THETA0) < 5]
            gc_found.append({
                'combo': [s['name'] for s in combo],
                'logMs': [s['logM'] for s in combo],
                'ncs':   [s['nc'] for s in combo],
                'nc_sum': sum(s['nc'] for s in combo),
                'lbs':   lbs,
                'sz': sz,
                'normal_lb': n_circ,
                'residuals': [round(r, 2) for r in residuals],
                'rms': round(rms, 2),
                'seps': seps,
                'near_poles': near_poles,
                'theta_matches': theta_matches
            })

gc_found.sort(key=lambda x: (x['rms'], -x['nc_sum'], -x['sz']))
lines.append(f"Getestete Kombinationen: {n_tested}  |  GK-Alignments (RMS≤{TOL_GC}°): {len(gc_found)}")

if gc_found:
    lines.append(f"\nTop-{min(25, len(gc_found))} Großkreis-Alignments:")
    for i, gc in enumerate(gc_found[:25]):
        seps_str = " → ".join(f"{s}°" for s in gc['seps'])
        ths_str = f"  *** GK-NORMAL nahe {gc['near_poles']}!" if gc['near_poles'] else ""
        tm_str  = f"  [Δ≈{gc['theta_matches']}° ≈ n×{THETA0}°]" if gc['theta_matches'] else ""
        lines.append(f"\n  [{i+1}] {' | '.join(gc['combo'])}  (n={gc['sz']})")
        lines.append(f"       GK-Normal: l={gc['normal_lb'][0]:.1f}°, b={gc['normal_lb'][1]:.1f}°  "
                     f"RMS={gc['rms']:.2f}°{ths_str}")
        lines.append(f"       Winkelsep: {seps_str}{tm_str}")
        lines.append(f"       logM:  {[f'{m:.1f}' for m in gc['logMs']]}  "
                     f"nc: {gc['ncs']}  Residuen: {gc['residuals']}°")
else:
    lines.append("  Kein Großkreis-Alignment mit RMS ≤ 7° gefunden.")
    lines.append("  → Versuch mit nc ≥ 0 (alle SMBHs):")
    cands_all = sorted(smbhs, key=lambda x: -x['logM'])[:30]
    for sz in range(3, 5):
        for combo in itertools.combinations(cands_all, sz):
            lbs = [(s['l'], s['b']) for s in combo]
            _, residuals, rms = great_circle_fit(lbs)
            if rms < TOL_GC:
                seps = [round(ang_lb(combo[i]['l'], combo[i]['b'],
                                     combo[i+1]['l'], combo[i+1]['b']), 1)
                        for i in range(len(combo) - 1)]
                lines.append(f"    GK (nc≥0): {[s['name'] for s in combo]}  RMS={rms:.2f}°  seps={seps}°")

# ── ABSCHNITT 2b: nc≥2 hochdichte Knoten – vollständige Combo-Suche ─────────
lines.append("\n\n══ 2b. 4+-OBJEKT ALIGNMENTS (nur nc≥2 SMBHs, vollständig) ══")
high_nc = [s for s in node_smbhs if s['nc'] >= 2]
lines.append(f"  Hochdichte SMBHs (nc≥2): {len(high_nc)}  → {[s['name'] for s in high_nc]}")
gc_high = []
for sz in range(4, min(7, len(high_nc) + 1)):
    for combo in itertools.combinations(high_nc, sz):
        too_close = any(ang_lb(combo[i]['l'], combo[i]['b'],
                               combo[j]['l'], combo[j]['b']) < MIN_SEP
                        for i in range(len(combo))
                        for j in range(i+1, len(combo)))
        if too_close:
            continue
        lbs = [(s['l'], s['b']) for s in combo]
        n_circ, residuals, rms = great_circle_fit(lbs)
        if n_circ and rms < TOL_GC:
            seps = [round(ang_lb(combo[i]['l'], combo[i]['b'],
                                 combo[i+1]['l'], combo[i+1]['b']), 1)
                    for i in range(len(combo) - 1)]
            pole_dists = {pn: ang_lb(n_circ[0], n_circ[1], pl, pb)
                          for pl, pb, pn in POLES_DEF}
            near_poles = [(pn, round(d, 1)) for pn, d in pole_dists.items()
                          if d < 20 or d > 160]
            theta_matches = sorted(set(
                round(s, 1) for s in seps
                if any(abs(s - k * THETA0) < 3.0 for k in [1, 2, 3])))
            nc_sum = sum(s['nc'] for s in combo)
            gc_high.append({
                'combo': [s['name'] for s in combo], 'sz': sz,
                'normal_lb': n_circ,
                'residuals': [round(r, 2) for r in residuals],
                'rms': round(rms, 2), 'seps': seps,
                'near_poles': near_poles, 'theta_matches': theta_matches,
                'nc_sum': nc_sum, 'ncs': [s['nc'] for s in combo],
                'logMs': [s['logM'] for s in combo],
            })
gc_high.sort(key=lambda x: (x['rms'], -x['nc_sum'], -x['sz']))
lines.append(f"  Gefundene 4+-Alignments (RMS≤{TOL_GC}°): {len(gc_high)}")
if gc_high:
    lines.append(f"\n  Top-{min(20, len(gc_high))} Hochdichte-Alignments (n≥4):")
    for i, gc in enumerate(gc_high[:20]):
        seps_str = " → ".join(f"{s}°" for s in gc['seps'])
        ths_str  = f"  *** GK-NORMAL nahe {gc['near_poles']}!" if gc['near_poles'] else ""
        tm_str   = f"  [Δ≈{gc['theta_matches']}° ≈ n×{THETA0}°]" if gc['theta_matches'] else ""
        lines.append(f"\n  [{i+1}] {' | '.join(gc['combo'])}  (n={gc['sz']})")
        lines.append(f"       GK-Normal: l={gc['normal_lb'][0]:.1f}°, b={gc['normal_lb'][1]:.1f}°  "
                     f"RMS={gc['rms']:.2f}°{ths_str}")
        lines.append(f"       Winkelsep: {seps_str}{tm_str}")
        lines.append(f"       logM: {[f'{m:.1f}' for m in gc['logMs']]}  nc: {gc['ncs']}  "
                     f"Residuen: {gc['residuals']}°")
else:
    lines.append("  Kein 4+-Objekt GK-Alignment unter nc≥2 SMBHs gefunden.")

# ── ABSCHNITT 3: Thales-Sphären-Kandidaten ───────────────────────────────────
lines.append("\n\n══ 3. THALES-SPHÄREN-KANDIDATEN (aus Galaxien-Katalog) ══")
if galaxies:
    for n in range(1, 8):
        Rn = R_S * n * THETA0 / 360.0
        cands_10 = [g for g in galaxies if abs(g['R_Thales'] - Rn) / Rn < 0.10]
        cands_nc  = [g for g in cands_10 if g['nc'] >= 1]
        lines.append(f"\n  n={n}  R_n={Rn:.0f} Mpc:  {len(cands_10)} Galaxien (±10%),  "
                     f"{len(cands_nc)} davon an HTM-Knoten (nc≥1)")
        if cands_nc:
            for g in sorted(cands_nc, key=lambda x: (-x['nc'], x['thales'][2]))[:8]:
                hits_str = '  '.join(f"{pn}(n={nn},δ={d}°)" for pn, nn, d in g['hits'])
                lines.append(f"    {g['name']:<20} {g['src']:>5}  d={g['d_mpc']:>5.0f}Mpc  "
                             f"R_Th={g['R_Thales']:>6.0f}Mpc  δ={100*g['thales'][2]:.1f}%  "
                             f"nc={g['nc']}  [{hits_str}]")
            # Nächste bekannte SMBHs zu diesen Kandidaten
            nearby = {}
            for gc_g in cands_nc:
                for s in smbhs:
                    sep = ang_lb(gc_g['l'], gc_g['b'], s['l'], s['b'])
                    if sep < NEIGHBOR_DEG:
                        if s['name'] not in nearby or nearby[s['name']]['sep'] > sep:
                            nearby[s['name']] = {'smbh': s, 'sep': round(sep, 1)}
            if nearby:
                lines.append(f"    ➜ Nächste SMBHs innerhalb {NEIGHBOR_DEG}° dieser Thales-Kandidaten:")
                for nm, info in sorted(nearby.items(), key=lambda x: x[1]['sep']):
                    s = info['smbh']
                    lines.append(f"      {s['name']:<22} logM={s['logM']:.2f}  nc={s['nc']}  "
                                 f"θ_D={s['theta_d']:.1f}°  sep={info['sep']}°")
                # Großkreis dieser SMBHs
                nb_list = [info['smbh'] for info in nearby.values() if info['smbh']['nc'] >= 1]
                if len(nb_list) >= 3:
                    lbs_nb = [(s['l'], s['b']) for s in nb_list[:10]]
                    for sz in range(min(MAX_COMBO, len(lbs_nb)), 2, -1):
                        found_gc = False
                        for combo_nb in itertools.combinations(range(len(lbs_nb)), sz):
                            lbs_sub = [lbs_nb[i] for i in combo_nb]
                            n_circ, residuals_nb, rms_nb = great_circle_fit(lbs_sub)
                            if rms_nb < TOL_GC:
                                names_nb = [nb_list[i]['name'] for i in combo_nb]
                                pole_d = {pn: ang_lb(n_circ[0], n_circ[1], pl, pb)
                                          for pl, pb, pn in POLES_DEF}
                                near_p = [(pn, round(d, 1)) for pn, d in pole_d.items()
                                          if d < 20 or d > 160]
                                lines.append(f"    *** GK für benachbarte SMBHs bei n={n}: "
                                             f"{names_nb}  RMS={rms_nb:.2f}°")
                                if near_p:
                                    lines.append(f"        GK-Normal nahe Pol: {near_p}")
                                found_gc = True
                        if found_gc:
                            break
else:
    lines.append("  Galaxien-Katalog nicht gefunden.")

# ── ABSCHNITT 4: 3D-Filament-Test ────────────────────────────────────────────
lines.append("\n\n══ 4. 3D-FILAMENT-TEST (SMBHs mit Distanz + nc ≥ 1) ══")
# 3D-Duplikate entfernen: max. einen Vertreter innerhalb 1 Mpc / innerhalb 1°
all_3d_raw = [s for s in smbhs if s['xyz3'] is not None and s['nc'] >= 1]
smbhs_3d = []
for s in all_3d_raw:
    too_close_3d = False
    for already in smbhs_3d:
        dx = s['xyz3'][0] - already['xyz3'][0]
        dy = s['xyz3'][1] - already['xyz3'][1]
        dz = s['xyz3'][2] - already['xyz3'][2]
        sep3d = math.sqrt(dx*dx + dy*dy + dz*dz)
        if sep3d < 2.0:   # innerhalb 2 Mpc = selber Cluster-Knoten
            too_close_3d = True; break
    if not too_close_3d:
        smbhs_3d.append(s)
lines.append(f"SMBHs mit nc≥1 und Distanz: {len(smbhs_3d)}")

if len(smbhs_3d) >= 3:
    best_filaments = []
    n_3d_tested = 0
    for sz in range(3, min(MAX_COMBO + 1, len(smbhs_3d) + 1)):
        combos = list(itertools.combinations(range(len(smbhs_3d)), sz))
        if n_3d_tested + len(combos) > 20000:
            break
        n_3d_tested += len(combos)
        for idxs in combos:
            pts = [smbhs_3d[i]['xyz3'] for i in idxs]
            # Nur sinnvoll wenn mind. 3 nicht kollineare Punkte (4+ für echten Test)
            if sz == 3:
                # Für 3 Punkte: prüfe ob sie wirklich einen Winkel bilden
                v01 = (pts[1][0]-pts[0][0], pts[1][1]-pts[0][1], pts[1][2]-pts[0][2])
                v02 = (pts[2][0]-pts[0][0], pts[2][1]-pts[0][1], pts[2][2]-pts[0][2])
                cp = cross3(v01, v02)
                if math.sqrt(dot3(cp, cp)) < 0.1 * math.sqrt(dot3(v01,v01)) * math.sqrt(dot3(v02,v02)):
                    continue  # casi collinear by construction
            dl, res_3d, rms_3d = line_fit_3d(pts)
            if dl and rms_3d < TOL_3D:
                names = [smbhs_3d[i]['name'] for i in idxs]
                dir_d = ang_lb(dl[0], dl[1], 305.0, 25.0)
                dir_d = min(dir_d, 180.0 - dir_d)   # 0-90° (undirected)
                best_filaments.append({
                    'names': names, 'sz': sz,
                    'direction': dl, 'rms': rms_3d,
                    'residuals': [round(r, 1) for r in res_3d],
                    'dir_dpol': round(dir_d, 1)
                })
    best_filaments.sort(key=lambda x: (x['rms'], -x['sz']))
    lines.append(f"Getestete 3D-Kombinationen: {n_3d_tested}  "
                 f"|  Filamente (RMS≤{TOL_3D} Mpc): {len(best_filaments)}")
    if best_filaments:
        lines.append(f"\nTop-{min(15, len(best_filaments))} 3D-Filamente:")
        for f in best_filaments[:15]:
            dpol_note = " *** ALONG D-POL AXIS" if f['dir_dpol'] < 15 else ""
            lines.append(f"  [{', '.join(f['names'])}]")
            lines.append(f"    Richtung l={f['direction'][0]:.1f}°, b={f['direction'][1]:.1f}°  "
                        f"RMS={f['rms']:.1f} Mpc  Δ(D-Pol-Achse)={f['dir_dpol']}°{dpol_note}")
            lines.append(f"    Residuen: {f['residuals']} Mpc")
    else:
        lines.append("  Kein 3D-Filament (RMS ≤ 20 Mpc) gefunden.")
        # Relax threshold
        best_rms = 999.9
        best_f = None
        for sz in range(3, min(5, len(smbhs_3d)+1)):
            for idxs in itertools.combinations(range(len(smbhs_3d)), sz):
                pts = [smbhs_3d[i]['xyz3'] for i in idxs]
                dl, res_3d, rms_3d = line_fit_3d(pts)
                if dl and rms_3d < best_rms:
                    best_rms = rms_3d
                    best_f = {'names': [smbhs_3d[i]['name'] for i in idxs],
                              'direction': dl, 'rms': rms_3d,
                              'residuals': [round(r, 1) for r in res_3d],
                              'dir_dpol': round(min(ang_lb(dl[0], dl[1], 305.0, 25.0),
                                                    180-ang_lb(dl[0], dl[1], 305.0, 25.0)), 1)}
        if best_f:
            lines.append(f"  Bestes 3D-Filament (relaxed): {best_f['names']}  RMS={best_f['rms']:.1f} Mpc  "
                         f"Δ(D-Pol)={best_f['dir_dpol']}°")

# ── ABSCHNITT 5: Super-Knoten (angular + Thales kombiniert) ──────────────────
lines.append("\n\n══ 5. KOMBINIERTE SUPER-KNOTEN ══")
lines.append("(nc ≥ 1 angular UND Thales-Sphäre <15%)")
if galaxies:
    super_nodes = [(g, g['thales']) for g in galaxies
                   if g['nc'] >= 1 and g['thales'] and g['thales'][2] < 0.15]
    super_nodes.sort(key=lambda x: (x[1][2], -x[0]['nc']))
    lines.append(f"  {len(super_nodes)} Galaxien an kombinierten Knoten")
    if super_nodes:
        lines.append(f"\n  {'Name':<20} {'Src':>5}  {'d':>6}  {'R_n':>6}  {'δ%':>5}  "
                     f"{'nc':>3}  Schalen")
        lines.append("  " + "-" * 65)
        for g, ths in super_nodes[:30]:
            hits_str = ' '.join(f"{pn}n{nn}" for pn, nn, d in g['hits'])
            lines.append(f"  {g['name']:<20} {g['src']:>5}  {g['d_mpc']:>6.0f}  "
                        f"{ths[1]:>6.0f}  {100*ths[2]:>5.1f}%  {g['nc']:>3}  [{hits_str}]")
        # GK für super-node SMBHs in der Nähe
        nearby_sn = {}
        for g, _ in super_nodes[:30]:
            for s in smbhs:
                sep = ang_lb(g['l'], g['b'], s['l'], s['b'])
                if sep < 15.0:
                    if s['name'] not in nearby_sn or nearby_sn[s['name']]['sep'] > sep:
                        nearby_sn[s['name']] = {'smbh': s, 'sep': round(sep, 1)}
        if nearby_sn:
            lines.append(f"\n  SMBHs in der Nähe der Super-Knoten (≤15°):")
            for nm, info in sorted(nearby_sn.items(), key=lambda x: x[1]['sep'])[:15]:
                s = info['smbh']
                lines.append(f"    {s['name']:<22} logM={s['logM']:.2f}  nc={s['nc']}  "
                             f"sep={info['sep']}°")

# ── ABSCHNITT 6: Geometrische Interpretation ─────────────────────────────────
lines.append("\n\n══ 6. GEOMETRISCHE INTERPRETATION ══")
lines.append(f"\n  Gesamt SMBHs:          {len(smbhs)}")
lines.append(f"  davon an Knoten nc≥1:  {len(node_smbhs)}")
lines.append(f"  davon an Knoten nc≥2:  {sum(1 for s in smbhs if s['nc'] >= 2)}")
lines.append(f"\n  Großkreis-Alignments: {len(gc_found)} (RMS ≤ {TOL_GC}°)")

if gc_found:
    best = gc_found[0]
    lines.append(f"\n  STÄRKSTES ALIGNMENT ({best['sz']} SMBHs):")
    lines.append(f"    Objekte: {' ↔ '.join(best['combo'])}")
    lines.append(f"    GK-Normal: l={best['normal_lb'][0]:.1f}°, b={best['normal_lb'][1]:.1f}°")
    lines.append(f"    RMS = {best['rms']:.2f}°  (≪ 7° = kollinear auf Großkreis)")
    lines.append(f"    Winkelabstände zwischen Objekten: {best['seps']} °")
    if best['near_poles']:
        lines.append(f"    *** Großkreis-Normal liegt nahe einem HTM-Pol: {best['near_poles']}")
        lines.append(f"    → Dieser Großkreis IST eine HTM-Schalen-Ebene!")
    if best['theta_matches']:
        lines.append(f"    → Winkelabstand ≈ n×{THETA0}° — HTM-Gitter-Spacing bestätigt!")
    # Separations commentary
    seps = best['seps']
    for sep in seps:
        ratio = sep / THETA0
        n_near = round(ratio)
        if abs(ratio - n_near) < 0.1 and n_near > 0:
            lines.append(f"    Sep {sep}° ≈ {n_near}×{THETA0}° — HTM-Resonanz!")

    # Mass trend along the alignment
    ms = best['logMs']
    if len(ms) >= 3:
        if ms[0] == max(ms) or ms[-1] == max(ms):
            lines.append(f"    Masseprofil: [{', '.join(f'{m:.1f}' for m in ms)}] — höchste Masse am Rand")
        elif ms[len(ms)//2] == max(ms):
            lines.append(f"    Masseprofil: [{', '.join(f'{m:.1f}' for m in ms)}] — höchste Masse in Mitte")
else:
    lines.append("\n  Kein klares Großkreis-Alignment (RMS ≤ 7°) unter nc≥1 SMBHs.")
    lines.append("  Mögliche Gründe:")
    lines.append("    1. SMBH-Katalog zu klein/lückenhaft (<100 Einträge)")
    lines.append("    2. Distanz-Messungen zu ungenau für 3D-Test")
    lines.append("    3. Echte isotrope Verteilung (Null-Ergebnis)")
    # Show closest miss
    if gc_found:
        pass  # already shown
    else:
        # Compute RMS for all 3-combinations of node smbhs and show best
        best_miss = None
        for combo in itertools.combinations(node_smbhs[:15], 3):
            lbs = [(s['l'], s['b']) for s in combo]
            _, _, rms = great_circle_fit(lbs)
            if best_miss is None or rms < best_miss[1]:
                best_miss = ([s['name'] for s in combo], rms)
        if best_miss:
            lines.append(f"\n  Nächster Miss (3 SMBHs): {best_miss[0]}  GK-RMS={best_miss[1]:.2f}°")

lines.append("\n" + "=" * 74)

# ── Ausgabe ────────────────────────────────────────────────────────────────────
txt = "\n".join(lines)
print(txt)
with open(OUT_TXT, "w", encoding="utf-8") as f:
    f.write(txt)
print(f"\nGespeichert: {OUT_TXT}")
