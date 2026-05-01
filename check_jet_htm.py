"""
SMBH Jet-Spinachsen vs HTM-Pole — Alignment-Test
=================================================
Hypothese (Tidal Torque Theory + HTM):
  Galaxien / SMBHs an HTM-Knoten (nc>=2) sollten Jet-Achsen
  bevorzugt in Richtung der nächsten HTM-Pole zeigen.

Methode:
  1) SMBH-Katalog (97 Quellen) → nc per Quelle
  2) Jet-PA aus Literatur (hardcoded) + MOJAVE J/AJ/147/143
  3) Für jede Quelle: projiziere nächsten HTM-Pol auf Himmel → erw. PA
  4) Alignment-Winkel δ = min(|jet_PA - erw_PA|, 180-|...|) mod 90
  5) nc≥2 vs nc=0: ist δ kleiner für hohe-nc Quellen?

HTM-Rahmen: θ₀=58.65°, TOL_ANG=4°, 6 Pole
"""
import math, os, sys, csv, ssl, random, time, re, io
import urllib.request, urllib.parse

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(RES, "OT_jet_htm.txt")
os.makedirs(RES, exist_ok=True)

_SSL = ssl.create_default_context()
_SSL.check_hostname = False
_SSL.verify_mode    = ssl.CERT_NONE

# ── HTM Konstanten ────────────────────────────────────────────────────────────
THETA0  = 58.65
TOL_ANG = 4.0

POLES_DEF = [
    (305.0,  25.0, 'D'),
    (125.0, -25.0, 'A'),
    ( 35.0,  25.0, 'S1'),
    (215.0, -25.0, 'S2'),
    (215.0,  25.0, 'S3'),
    ( 35.0, -25.0, 'S4'),
]

# ── Koordinaten ───────────────────────────────────────────────────────────────
def radec_to_lb(ra, dec):
    ra_r = math.radians(ra); dc_r = math.radians(dec)
    RN   = math.radians(192.85948); DN = math.radians(27.12825)
    b  = math.asin(max(-1., min(1.,
         math.sin(dc_r)*math.sin(DN) +
         math.cos(dc_r)*math.cos(DN)*math.cos(ra_r - RN))))
    y  = math.cos(dc_r)*math.sin(ra_r - RN)
    x  = math.cos(dc_r)*math.sin(DN)*math.cos(ra_r - RN) - math.sin(dc_r)*math.cos(DN)
    l  = (math.degrees(math.atan2(y, x)) + 122.93192) % 360.0
    return l, math.degrees(b)

def lb_to_radec(l, b):
    """Galaktisch → J2000 RA, Dec (Grad)"""
    lr, br = math.radians(l), math.radians(b)
    RN  = math.radians(192.85948); DN = math.radians(27.12825)
    LN  = math.radians(122.93192)
    sin_d = (math.sin(br)*math.sin(DN) +
             math.cos(br)*math.cos(DN)*math.cos(LN - lr))
    dec_r  = math.asin(max(-1., min(1., sin_d)))
    y      = math.cos(br)*math.sin(LN - lr)
    x      = (math.sin(br)*math.cos(DN) -
               math.cos(br)*math.sin(DN)*math.cos(LN - lr))
    ra_r   = math.atan2(y, x) + RN
    return math.degrees(ra_r) % 360, math.degrees(dec_r)

def lb_xyz(l, b):
    lr, br = math.radians(l), math.radians(b)
    return (math.cos(br)*math.cos(lr), math.cos(br)*math.sin(lr), math.sin(br))

def ang_sep_xyz(v1, v2):
    d = max(-1., min(1., sum(a*b for a,b in zip(v1, v2))))
    return math.degrees(math.acos(d))

def crossing_density(l, b):
    sv = lb_xyz(l, b)
    nc, which = 0, []
    for pl, pb, name in POLES_DEF:
        pv  = lb_xyz(pl, pb)
        th  = ang_sep_xyz(pv, sv)
        for n in [1, 2, 3]:
            if abs(th - n*THETA0) < TOL_ANG:
                nc += 1
                which.append((name, n, th))
    return nc, which

def pole_projected_pa(ra_src, dec_src, ra_pole, dec_pole):
    """
    Positionswinkel (°, N->E, äquatorial, mod 180 = Achse)
    der Richtung Quelle → Pol auf der Himmelssphäre.
    Benutzt Großkreis-Formel (Vincenty-style).
    """
    d = math.pi / 180
    dra   = (ra_pole - ra_src) * d
    dec_s = dec_src  * d
    dec_p = dec_pole * d
    y = math.sin(dra) * math.cos(dec_p)
    x = (math.cos(dec_s)*math.sin(dec_p) -
         math.sin(dec_s)*math.cos(dec_p)*math.cos(dra))
    pa = math.degrees(math.atan2(y, x)) % 180
    return pa

def alignment_angle(jet_pa, expected_pa):
    """Minimaler Winkel zw. zwei Achsen (je mod 180), Ergebnis 0-90°."""
    diff = abs((jet_pa % 180) - (expected_pa % 180))
    if diff > 90:
        diff = 180 - diff
    return diff

def crossproduct_expected_pa(ra_src, dec_src, which):
    """
    Fragmentierungs-Modell: Jet-Achse = P1 × P2 (Kreuzprodukt der Pol-Vektoren).
    Physik: Zwei Schalen kollisionieren am nc>=2 Knoten. Der Drehimpulsvektor
    der Kollision steht senkrecht auf beiden Pol-Richtungen — die 'freie Richtung'.

    Für antipodal-Paare (S1+S2, D+A, S3+S4 mit P2=-P1) ist der Kreuzprodukt = 0
    → keine eindeutige Richtung. Bei nc>=3 wird das Paar mit dem größten
    |P1×P2| gewählt (maximum Orthogonalität).

    Gibt (erw_PA_deg, polname_string) oder (None, None) zurück.
    """
    # Sammle alle beteiligten eindeutigen Pol-Richtungsvektoren
    seen, pole_vecs = set(), []
    for pname, n, theta in which:
        if pname in seen:
            continue
        seen.add(pname)
        for pl_l, pl_b, p in POLES_DEF:
            if p == pname:
                pole_vecs.append((lb_xyz(pl_l, pl_b), pname))
                break
    if len(pole_vecs) < 2:
        return None, None
    # Wähle das Paar mit dem größten |P1×P2| (am wenigsten kollinear/antipodal)
    best_mag = 0.0
    best_cx = best_cy = best_cz = 0.0
    best_n1 = best_n2 = ''
    for i in range(len(pole_vecs)):
        for j in range(i+1, len(pole_vecs)):
            v1, n1 = pole_vecs[i]
            v2, n2 = pole_vecs[j]
            cx = v1[1]*v2[2] - v1[2]*v2[1]
            cy = v1[2]*v2[0] - v1[0]*v2[2]
            cz = v1[0]*v2[1] - v1[1]*v2[0]
            mag = math.sqrt(cx*cx + cy*cy + cz*cz)
            if mag > best_mag:
                best_mag = mag
                best_cx, best_cy, best_cz = cx, cy, cz
                best_n1, best_n2 = n1, n2
    if best_mag < 0.1:  # Alle Paare fast antipodal (sin(angle) < 0.1 ≈ <6°)
        return None, None
    best_cx /= best_mag; best_cy /= best_mag; best_cz /= best_mag
    # Projiziere auf Himmelstangentialebene bei (ra_src, dec_src)
    ra_r  = math.radians(ra_src)
    dec_r = math.radians(dec_src)
    n_hat = (-math.sin(dec_r)*math.cos(ra_r),
             -math.sin(dec_r)*math.sin(ra_r),
              math.cos(dec_r))
    e_hat = (-math.sin(ra_r), math.cos(ra_r), 0.0)
    proj_n = best_cx*n_hat[0] + best_cy*n_hat[1] + best_cz*n_hat[2]
    proj_e = best_cx*e_hat[0] + best_cy*e_hat[1] + best_cz*e_hat[2]
    if abs(proj_n) < 1e-10 and abs(proj_e) < 1e-10:
        return None, None
    pa = math.degrees(math.atan2(proj_e, proj_n)) % 180
    return pa, f"{best_n1}x{best_n2}"

def nearest_active_pole_pos(l_src, b_src, which):
    """
    Gibt (ra_pole, dec_pole, polname) des geometrisch nächsten
    aktiven HTM-Rings zurück (oder des insgesamt nächsten Pols).
    """
    if which:
        sv = lb_xyz(l_src, b_src)
        best_d, best_pl, best_pb, best_n = 999, None, None, ''
        for pname, n, theta in which:
            for pl_l, pl_b, p in POLES_DEF:
                if p != pname:
                    continue
                pv = lb_xyz(pl_l, pl_b)
                d  = abs(ang_sep_xyz(pv, sv) - n*THETA0)
                if d < best_d:
                    best_d = d; best_pl = pl_l; best_pb = pl_b; best_n = pname
        if best_pl is not None:
            ra_p, dec_p = lb_to_radec(best_pl, best_pb)
            return ra_p, dec_p, best_n
    # Fallback: geometrisch nächster Pol
    sv = lb_xyz(l_src, b_src)
    best_d, best = 999, None
    for pl_l, pl_b, pname in POLES_DEF:
        pv = lb_xyz(pl_l, pl_b)
        d  = ang_sep_xyz(pv, sv)
        if d < best_d:
            best_d = d; best = (pl_l, pl_b, pname)
    ra_p, dec_p = lb_to_radec(best[0], best[1])
    return ra_p, dec_p, best[2]

# ── Hardcoded Jet-PA Literaturwerte ───────────────────────────────────────────
# jet_pa = Richtung des Jets auf dem Himmel (PA, Grad N→E, 0-360°)
# Quellen: VLBI / VLA Literatur
JET_PA_LIT = {
    # Virgo-Bereich
    'M 87':      {'pa': 291.0, 'ref': 'Biretta+1999'},
    'NGC 4486':  {'pa': 291.0, 'ref': 'Biretta+1999'},
    'M 84':      {'pa': 308.0, 'ref': 'Laing+1986'},
    'NGC 4374':  {'pa': 308.0, 'ref': 'Laing+1986'},
    'NGC 4261':  {'pa':  88.0, 'ref': 'Biretta+1992'},
    'NGC 4278':  {'pa': 116.0, 'ref': 'Giroletti+2005'},
    'NGC 4552':  {'pa': 110.0, 'ref': 'Filho+2004'},
    'M 89':      {'pa': 110.0, 'ref': 'Filho+2004'},
    'NGC 4649':  {'pa':  83.0, 'ref': 'Shurkin+2008'},
    'NGC 4636':  {'pa': 142.0, 'ref': 'Jones+1997'},
    # Perseus
    'NGC 1275':  {'pa': 160.0, 'ref': 'Walker+1994'},
    'NGC 1052':  {'pa':  65.0, 'ref': 'Vermeulen+2003'},
    # Centaurus
    'NGC 5128':  {'pa':  51.0, 'ref': 'Burns+1983'},
    'Cen A':     {'pa':  51.0, 'ref': 'Burns+1983'},
    # Andere wohlbekannte Radio-Galaxien
    'NGC 6251':  {'pa': 251.0, 'ref': 'Sudou+2000'},
    'NGC 315':   {'pa': 175.0, 'ref': 'Cotton+1999'},
    'NGC 383':   {'pa':  15.0, 'ref': 'Worrall+2007'},
    'NGC 1265':  {'pa': 125.0, 'ref': "O'Dea+1987"},
    'NGC 1600':  {'pa':  68.0, 'ref': 'Dunn+2010'},
    'NGC 3115':  {'pa':  50.0, 'ref': 'Wrobel+2012'},
    'NGC 3379':  {'pa': 155.0, 'ref': 'Nyland+2016'},
    'NGC 3608':  {'pa': 100.0, 'ref': 'Nyland+2016'},
    'NGC 4168':  {'pa':  36.0, 'ref': 'Capetti+2000'},
    'NGC 4459':  {'pa':  90.0, 'ref': 'Nyland+2016'},
    'NGC 4473':  {'pa':  92.0, 'ref': 'Nyland+2016'},
    'NGC 4564':  {'pa':  47.0, 'ref': 'Krajnovic+2011'},
    'NGC 4594':  {'pa': 100.0, 'ref': 'Bajaja+1988'},
    'M 104':     {'pa': 100.0, 'ref': 'Bajaja+1988'},
    'NGC 4697':  {'pa':  65.0, 'ref': 'Sarazin+2000'},
    'NGC 4742':  {'pa':   5.0, 'ref': 'Krajnovic+2011'},
    'NGC 4889':  {'pa':  78.0, 'ref': 'Sohrab+2015'},
    'NGC 5845':  {'pa':  81.0, 'ref': 'Nyland+2016'},
    'NGC 7052':  {'pa': 100.0, 'ref': 'vandenBosch+1998'},
    'NGC 7768':  {'pa':  30.0, 'ref': 'Bettoni+1990'},
    'NGC 821':   {'pa':  50.0, 'ref': 'Nyland+2016'},
}

# ── SMBH-Katalog laden ───────────────────────────────────────────────────────
print("=" * 65)
print("SMBH Jet-Spin ↔ HTM-Pol Alignment-Test (check_jet_htm.py)")
print("=" * 65)

cat_file = os.path.join(RES, "catalogs", "smbh_extended.csv")
smbhs = []
with open(cat_file, newline='', encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try:
            ra  = float(r['RA_deg'])
            dec = float(r['Dec_deg'])
        except Exception:
            continue
        l, b = radec_to_lb(ra, dec)
        nc, which = crossing_density(l, b)
        smbhs.append({'name': r['Name'].strip(), 'ra': ra, 'dec': dec,
                      'l': l, 'b': b, 'nc': nc, 'which': which,
                      'logM': r.get('logMbh', '').strip()})

print(f"SMBH-Katalog: {len(smbhs)} Eintraege geladen")

# ── VizieR: MOJAVE J/AJ/147/143 ──────────────────────────────────────────────
def vizier_query(source, out_cols, constraints='', max_out=3000):
    params = {'-source': source, '-out': out_cols,
              '-out.max': str(max_out), '-oc.form': 'dec'}
    url = 'https://vizier.cds.unistra.fr/viz-bin/votable?' + urllib.parse.urlencode(params)
    if constraints:
        url += '&' + constraints
    for attempt in range(3):
        try:
            if attempt:
                time.sleep(5*attempt)
            req = urllib.request.Request(url, headers={'User-Agent':'Python/HTMJet1.0'})
            with urllib.request.urlopen(req, timeout=90, context=_SSL) as r:
                raw = r.read().decode('utf-8','replace')
            rows = []
            for tm in re.finditer(r'<TABLE\b[^>]*>(.*?)</TABLE>', raw, re.DOTALL):
                tb = tm.group(1)
                fields = re.findall(r'<FIELD[^>]+name="([^"]+)"', tb)
                fields = [f for f in fields if f != 'recno']
                if len(fields) < 2:
                    continue
                for tr in re.finditer(r'<TR>(.*?)</TR>', tb, re.DOTALL):
                    tds = re.findall(r'<TD>(.*?)</TD>', tr.group(1))
                    if len(tds) < 2:
                        continue
                    rows.append({fields[i]: (tds[i].strip() if i < len(tds) else '')
                                 for i in range(len(fields))})
            return rows
        except Exception as e:
            print(f"  Versuch {attempt+1} fehlgeschlagen: {e}")
    return []

print("\n[VizieR] MOJAVE J/AJ/147/143 table1 (Quellpositionen)...")
moj1_rows = vizier_query('J/AJ/147/143/table1', 'Name _RA _DE z O', max_out=300)
print(f"  table1: {len(moj1_rows)} Zeilen")

print("[VizieR] MOJAVE J/AJ/147/143 table2 (Jet-Komponenten)...")
moj2_rows = vizier_query('J/AJ/147/143/table2', 'Name m_Name r PA', max_out=2000)
print(f"  table2: {len(moj2_rows)} Zeilen")

# Positionen aus table1
def to_float(s):
    try:
        return float(s.strip())
    except Exception:
        return None

moj_pos = {}
for row in moj1_rows:
    nm  = row.get('Name','').strip()
    ra  = to_float(row.get('_RA',''))
    dec = to_float(row.get('_DE',''))
    if nm and ra is not None and dec is not None:
        moj_pos[nm] = (ra, dec)

# Inner-Jet-PA aus table2: nächste Nicht-Kern-Komponente(n)
inner_pas = {}
for row in moj2_rows:
    nm  = row.get('Name','').strip()
    try:
        comp = int(row.get('m_Name','0') or '0')
    except Exception:
        comp = 0
    r_v = to_float(row.get('r',''))
    pa  = to_float(row.get('PA',''))
    if comp == 0 or r_v is None or r_v <= 0 or pa is None:
        continue
    inner_pas.setdefault(nm, []).append((r_v, pa))

jet_pa_mojave = {}
for nm, comps in inner_pas.items():
    comps.sort(key=lambda x: x[0])
    plist = [c[1] for c in comps[:3]]  # bis zu 3 nächste Komponenten
    if not plist:
        continue
    # Kreismittelwert: Achse (mod 180) via Verdopplung
    sins = sum(math.sin(2*math.radians(p)) for p in plist)
    coss = sum(math.cos(2*math.radians(p)) for p in plist)
    pa_ax = math.degrees(math.atan2(sins, coss)/2) % 180
    jet_pa_mojave[nm] = pa_ax

print(f"  Inner-Jet-PA berechnet fuer {len(jet_pa_mojave)} MOJAVE-Quellen")

# ── Zusammenführen ────────────────────────────────────────────────────────────
all_jets = []  # Liste von Dicts

# A) SMBH-Katalog × Literatur-PA
for s in smbhs:
    nm = s['name']
    pa_info = None
    nm_uc = nm.upper().replace(' ','').replace('_','')
    for key, val in JET_PA_LIT.items():
        if key.upper().replace(' ','') == nm_uc:
            pa_info = val
            break
    if pa_info:
        all_jets.append({
            'name': nm, 'ra': s['ra'], 'dec': s['dec'],
            'l': s['l'], 'b': s['b'], 'nc': s['nc'], 'which': s['which'],
            'jet_pa': pa_info['pa'] % 180,
            'logM': s['logM'],
            'src': 'lit:' + pa_info['ref']
        })

# B) MOJAVE-Quellen (nicht bereits in A)
existing = {j['name'].upper().replace(' ','') for j in all_jets}
for nm, (ra, dec) in moj_pos.items():
    if nm.upper().replace(' ','') in existing:
        continue
    if nm not in jet_pa_mojave:
        continue
    l, b = radec_to_lb(ra, dec)
    nc, which = crossing_density(l, b)
    all_jets.append({
        'name': nm, 'ra': ra, 'dec': dec,
        'l': l, 'b': b, 'nc': nc, 'which': which,
        'jet_pa': jet_pa_mojave[nm],
        'logM': '',
        'src': 'MOJAVE'
    })

print(f"\nGesamt: {len(all_jets)} Quellen mit Jet-PA")
print(f"  Literatur: {sum(1 for j in all_jets if j['src'].startswith('lit'))}")
print(f"  MOJAVE:    {sum(1 for j in all_jets if j['src'] == 'MOJAVE')}")

# ── Alignment-Test ────────────────────────────────────────────────────────────
print("\n[TEST] Berechne Alignment-Winkel jet_PA vs Pol-Projektion...")

for j in all_jets:
    ra_p, dec_p, pname = nearest_active_pole_pos(j['l'], j['b'], j['which'])
    j['ra_pole']     = ra_p
    j['dec_pole']    = dec_p
    j['pole_name']   = pname
    j['expected_pa'] = pole_projected_pa(j['ra'], j['dec'], ra_p, dec_p)
    j['delta']       = alignment_angle(j['jet_pa'], j['expected_pa'])
    # Kreuzprodukt-Test (Fragmentierungs-Modell: jet || P1×P2)
    cross_pa, cross_poles = crossproduct_expected_pa(j['ra'], j['dec'], j['which'])
    j['cross_pa']    = cross_pa
    j['cross_poles'] = cross_poles
    j['cross_delta'] = alignment_angle(j['jet_pa'], cross_pa) if cross_pa is not None else None

ALIGN_THRESH = 30  # Grad: unter dieser Schwelle = ausgerichtet

def group_stats(jets, label):
    if not jets:
        print(f"  {label}: N=0")
        return
    n     = len(jets)
    n_al  = sum(1 for j in jets if j['delta'] < ALIGN_THRESH)
    delts = sorted(j['delta'] for j in jets)
    m_d   = sum(delts)/n
    med   = delts[n//2]
    print(f"  {label:10s}: N={n:4d}  "
          f"ausgerichtet(<{ALIGN_THRESH}deg)={n_al:3d}/{n} ({100*n_al/n:4.1f}%)  "
          f"Mittel={m_d:5.1f}deg  Median={med:5.1f}deg")

print(f"\n  Alignment-Schwelle: < {ALIGN_THRESH}°  (Zufallserwartung: {100*ALIGN_THRESH/90:.1f}%)")
group_stats([j for j in all_jets if j['nc']==0], "nc=0")
group_stats([j for j in all_jets if j['nc']==1], "nc=1")
group_stats([j for j in all_jets if j['nc']>=2], "nc>=2")
group_stats([j for j in all_jets if j['nc']>=3], "nc>=3")
group_stats(all_jets, "alle")

# ── Monte Carlo ────────────────────────────────────────────────────────────────
random.seed(42)
N_MC = 50000
nc2 = [j for j in all_jets if j['nc'] >= 2]

mc_result = {}
if nc2:
    n_obs = sum(1 for j in nc2 if j['delta'] < ALIGN_THRESH)
    frac_obs = n_obs / len(nc2)
    # MC: zufaellige Jet-PA und Pol-PA → Referenzverteilung
    mc_hits = sum(
        1 for _ in range(N_MC)
        if alignment_angle(random.uniform(0,180), random.uniform(0,180)) < ALIGN_THRESH
    )
    frac_mc = mc_hits / N_MC
    ratio   = frac_obs / frac_mc if frac_mc > 0 else float('nan')
    mc_result = dict(frac_obs=frac_obs, frac_mc=frac_mc, ratio=ratio,
                     n_nc2=len(nc2), n_obs=n_obs)
    flag = "*** SIGNAL ***" if ratio > 1.5 else ("(kein Signal)" if ratio < 0.8 else "(schwach)")
    print(f"\n  MC (nc>=2): obs={100*frac_obs:.1f}%  MC={100*frac_mc:.1f}%  "
          f"Ratio={ratio:.2f}  {flag}")

# ── Kreuzprodukt-Test (Fragmentierungs-Modell) ───────────────────────────────
print(f"\n[KREUZPRODUKT-TEST] Fragmentierungs-Modell: Jet || P1 x P2")
print(f"  Vorhersage: Jet senkrecht zu beiden Schalen-Polen (freie Richtung)")

nc2_cross = [j for j in all_jets if j['nc'] >= 2 and j['cross_delta'] is not None]
nc0_cross  = [j for j in all_jets if j['nc'] == 0 and j['cross_delta'] is not None]
nc1_cross  = [j for j in all_jets if j['nc'] == 1 and j['cross_delta'] is not None]

def group_stats_cross(jets, label):
    if not jets:
        print(f"  {label}: N=0")
        return
    n    = len(jets)
    nal  = sum(1 for j in jets if j['cross_delta'] < ALIGN_THRESH)
    ds   = sorted(j['cross_delta'] for j in jets)
    md   = sum(ds)/n
    med  = ds[n//2]
    print(f"  {label:10s}: N={n:4d}  "
          f"ausgerichtet(<{ALIGN_THRESH}deg)={nal:3d}/{n} ({100*nal/n:4.1f}%)  "
          f"Mittel={md:5.1f}deg  Median={med:5.1f}deg")

print(f"  Zufallserwartung: {100*ALIGN_THRESH/90:.1f}%")
group_stats_cross(nc0_cross,  "nc=0  (ref)")
group_stats_cross(nc1_cross,  "nc=1  (ref)")
group_stats_cross(nc2_cross,  "nc>=2 (test)")

mc_cross = {}
if nc2_cross:
    n_obs_c  = sum(1 for j in nc2_cross if j['cross_delta'] < ALIGN_THRESH)
    frac_obs_c = n_obs_c / len(nc2_cross)
    mc_hits_c  = sum(
        1 for _ in range(N_MC)
        if alignment_angle(random.uniform(0,180), random.uniform(0,180)) < ALIGN_THRESH
    )
    frac_mc_c = mc_hits_c / N_MC
    ratio_c   = frac_obs_c / frac_mc_c if frac_mc_c > 0 else float('nan')
    mc_cross  = dict(frac_obs=frac_obs_c, frac_mc=frac_mc_c, ratio=ratio_c,
                     n=len(nc2_cross), n_obs=n_obs_c)
    flag_c = "*** SIGNAL ***" if ratio_c > 1.5 else ("(kein Signal)" if ratio_c < 0.8 else "(schwach)")
    print(f"\n  MC (nc>=2, Kreuzprodukt): obs={100*frac_obs_c:.1f}%  "
          f"MC={100*frac_mc_c:.1f}%  Ratio={ratio_c:.2f}  {flag_c}")

# Vergleichs-Tabelle: welcher Test ist besser?
print(f"\n  Modell-Vergleich (nc>=2):")
print(f"  {'Modell':35s}  {'Ratio':>6}  Interpretation")
if mc_result:
    f1 = "Jet -> naechster Pol (Tidal Torque)"
    print(f"  {f1:35s}  {mc_result['ratio']:6.2f}  "
          + ("Signal" if mc_result['ratio']>1.5 else
             "Anti-Alignment" if mc_result['ratio']<0.8 else "kein Signal"))
if mc_cross:
    f2 = "Jet || P1xP2 (Fragmentierungs-Modell)"
    print(f"  {f2:35s}  {mc_cross['ratio']:6.2f}  "
          + ("Signal" if mc_cross['ratio']>1.5 else
             "Anti-Alignment" if mc_cross['ratio']<0.8 else "kein Signal"))

# Detail nc>=2 Kreuzprodukt-Test
if nc2_cross:
    print(f"\n  nc>=2 Kreuzprodukt-Detail (nach cross_delta sortiert):")
    print(f"  {'Name':22s} nc  {'Pole':8s}  JetPA  P1xP2-PA  cross_delta  Quelle")
    for j in sorted(nc2_cross, key=lambda x: x['cross_delta']):
        mk = 'ok' if j['cross_delta'] < ALIGN_THRESH else '  '
        print(f"  {j['name']:22s} {j['nc']:2d}  {j['cross_poles']:8s}  "
              f"{j['jet_pa']:5.1f}  {j['cross_pa']:8.1f}  {j['cross_delta']:10.1f}  "
              f"{mk} {j['src']}")

# ── Detail-Ausgabe ────────────────────────────────────────────────────────────
nc2_sorted = sorted((j for j in all_jets if j['nc']>=2), key=lambda x: x['delta'])
if nc2_sorted:
    print(f"\n  nc>=2 Detail (nach δ sortiert):")
    print(f"  {'Name':22s} nc  {'Ringe':20s} JetPA  ExpPA  δ       Quelle")
    for j in nc2_sorted:
        rings = ','.join(f"{p}n{n}" for p,n,_ in j['which'])
        mk = '✓' if j['delta'] < ALIGN_THRESH else ' '
        print(f"  {j['name']:22s} {j['nc']:2d}  {rings:20s} "
              f"{j['jet_pa']:5.1f}° {j['expected_pa']:5.1f}° {j['delta']:5.1f}° "
              f"{mk} {j['src']}")

# Bin-Verteilung
print(f"\n  δ-Verteilung (bins 15°):")
bins = [(0,15),(15,30),(30,45),(45,60),(60,75),(75,90)]
g = {'nc=0':  [j for j in all_jets if j['nc']==0],
     'nc=1':  [j for j in all_jets if j['nc']==1],
     'nc>=2': [j for j in all_jets if j['nc']>=2]}
header = f"  {'Bin':9s}  " + "  ".join(f"{lbl:15s}" for lbl in g)
print(header)
for lo, hi in bins:
    parts = []
    for lbl, lst in g.items():
        n = len(lst)
        c = sum(1 for j in lst if lo <= j['delta'] < hi)
        parts.append(f"{c:3d}/{n}={100*c/max(1,n):4.0f}%")
    print(f"  {lo:2d}-{hi:2d}deg    " + "    ".join(parts))

# ── Ergebnisdatei schreiben ───────────────────────────────────────────────────
buf = io.StringIO()

def w(s=''):
    buf.write(str(s) + '\n')

w("=" * 72)
w("OT-Jet: SMBH Jet-Spinachsen vs HTM-Pol Alignment-Test")
w("=" * 72)
w(f"Datum:       {time.strftime('%Y-%m-%d %H:%M')}")
w(f"HTM:         theta0={THETA0}deg, TOL_ANG={TOL_ANG}deg")
w(f"Pole:        D(305,25), A(125,-25), S1(35,25), S2(215,-25), S3(215,25), S4(35,-25)")
w(f"Alignment-Schwelle: delta < {ALIGN_THRESH}deg  (Zufallserwartung: {100*ALIGN_THRESH/90:.1f}%)")
w()
w(f"Quellen mit Jet-PA gesamt: {len(all_jets)}")
w(f"  Literatur: {sum(1 for j in all_jets if j['src'].startswith('lit'))}")
w(f"  MOJAVE:    {sum(1 for j in all_jets if j['src']=='MOJAVE')}")
w()
w("-"*72)
w("Statistik:")
for grp, cond in [('nc=0', lambda j: j['nc']==0),
                  ('nc=1', lambda j: j['nc']==1),
                  ('nc>=2',lambda j: j['nc']>=2),
                  ('nc>=3',lambda j: j['nc']>=3),
                  ('alle', lambda j: True)]:
    lst = [j for j in all_jets if cond(j)]
    if not lst:
        continue
    n   = len(lst)
    nal = sum(1 for j in lst if j['delta'] < ALIGN_THRESH)
    ds  = sorted(j['delta'] for j in lst)
    w(f"  {grp:8s}: N={n:4d}  "
      f"ausgerichtet={nal}/{n}({100*nal/max(1,n):4.1f}%)  "
      f"Mittel={sum(ds)/n:.1f}deg  Median={ds[n//2]:.1f}deg")
w()
if mc_result:
    w(f"Monte Carlo (nc>=2, N_MC={N_MC}):")
    w(f"  Beobachtet:  {100*mc_result['frac_obs']:.1f}% ausgerichtet")
    w(f"  Zufall:      {100*mc_result['frac_mc']:.1f}% erwartet")
    w(f"  Ratio:       {mc_result['ratio']:.2f}x")
    flag2 = "SIGNAL (>1.5x)" if mc_result['ratio']>1.5 else (
            "anti-alignment (<0.8x)" if mc_result['ratio']<0.8 else "kein Signal")
    w(f"  Bewertung:   {flag2}")
w()
w("-"*72)
w("KREUZPRODUKT-TEST (Fragmentierungs-Modell: Jet || P1 x P2)")
w("Physik: Schalen-Kollision an nc>=2 Knoten → Drehimpuls = P1 x P2 = 'freie Richtung'")
w()
for grp, lst in [('nc=0 (ref)', nc0_cross),('nc=1 (ref)', nc1_cross),('nc>=2 (test)', nc2_cross)]:
    if not lst:
        continue
    n   = len(lst)
    nal = sum(1 for j in lst if j['cross_delta'] < ALIGN_THRESH)
    ds  = sorted(j['cross_delta'] for j in lst)
    w(f"  {grp:14s}: N={n:4d}  "
      f"ausgerichtet={nal}/{n}({100*nal/max(1,n):4.1f}%)  "
      f"Mittel={sum(ds)/n:.1f}deg  Median={ds[n//2]:.1f}deg")
w()
if mc_cross:
    w(f"Monte Carlo (nc>=2, Kreuzprodukt, N_MC={N_MC}):")
    w(f"  Beobachtet:  {100*mc_cross['frac_obs']:.1f}% ausgerichtet")
    w(f"  Zufall:      {100*mc_cross['frac_mc']:.1f}% erwartet")
    w(f"  Ratio:       {mc_cross['ratio']:.2f}x")
    flag3 = "SIGNAL (>1.5x)" if mc_cross['ratio']>1.5 else (
            "anti-alignment (<0.8x)" if mc_cross['ratio']<0.8 else "kein Signal")
    w(f"  Bewertung:   {flag3}")
w()
w("Modell-Vergleich (nc>=2):")
if mc_result:
    w(f"  Jet->Pol (Tidal Torque):          Ratio={mc_result['ratio']:.2f}")
if mc_cross:
    w(f"  Jet||P1xP2 (Fragmentierung):      Ratio={mc_cross['ratio']:.2f}")
w()
if nc2_cross:
    w(f"{'Name':22s} {'nc':>3}  {'Pole':8s}  {'JetPA':>7}  {'P1xP2-PA':>9}  {'cross_d':>7}  Quelle")
    w("-"*80)
    for j in sorted(nc2_cross, key=lambda x: x['cross_delta']):
        mk = '*' if j['cross_delta'] < ALIGN_THRESH else ' '
        w(f"  {j['name']:20s} {j['nc']:3d}  {j['cross_poles']:8s}  "
          f"{j['jet_pa']:6.1f}  {j['cross_pa']:8.1f}  {j['cross_delta']:6.1f}  {j['src']} {mk}")
w()
w("-"*72)
w(f"{'Name':22s} {'nc':>3}  {'Ringe':22s}  {'JetPA':>7}  {'ExpPA':>7}  {'delta':>6}  Quelle")
w("-"*80)
for j in sorted(all_jets, key=lambda x: (-x['nc'], x['delta'])):
    rings = ','.join(f"{p}n{n}" for p,n,_ in j['which']) if j['which'] else '-'
    mk = '*' if j['delta'] < ALIGN_THRESH and j['nc']>=2 else ' '
    w(f"  {j['name']:20s} {j['nc']:3d}  {rings:22s}  "
      f"{j['jet_pa']:6.1f}  {j['expected_pa']:6.1f}  {j['delta']:6.1f}  {j['src']} {mk}")

txt = buf.getvalue()
with open(OUT, 'w', encoding='utf-8') as fh:
    fh.write(txt)

print(f"\nAusgabe: {OUT}")
print("FERTIG.")
