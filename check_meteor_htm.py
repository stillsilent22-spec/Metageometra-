"""
Meteor / Fireball HTM-Knoten-Test
===================================
Findet "Nächte mit vielen Himmelskörpern" (Meteorschwärme + Feuerbälle) und
prüft, ob ihre RADIANTE-Positionen und PFADE bevorzugt an HTM-Winkelknoten liegen.

Test 1 – Schauer-Radiant-Test
    Etablierte IAU/IMO-Meteorschwärme (hardcoded aus MDC/IMO-Kalender):
    nc(Radiant) > Zufallserwartung?

Test 2 – Feuerball-Radiant-Test
    NASA CNEOS-Feuerbälle mit ECI-Geschwindigkeitskomponenten (alle ~350 Ereignisse):
    Radiant = Anti-Velocity-Richtung in J2000.
    nc(Radiant) + nc(Anti-Radiant) vs MC-Vergleich.

Test 3 – Pfad-Test (Großkreisbögen ab Radiant)
    Für jeden Feuerball: sample 36 Richtungen × 4 Abstände (15/30/45/60°) ab
    Radiant → Prüfe ob die Punkte entlang typischer Meteor-Pfade HTM-Knoten
    kreuzen häufiger als erwartet.

HTM-Rahmen: θ₀=58.65°, R_S=4228.3 Mpc, D-Pol=(305°,25°), TOL_ANG=4°

Daten:
    CNEOS: https://ssd-api.jpl.nasa.gov/fireball.api?vel-comp=true&req-vel-comp=true
    IMO/MDC: hardcoded (J2000 RA,Dec in Grad)
"""
import math, os, sys, json, random, urllib.request, time, ssl
from collections import Counter

# SSL-Kontext ohne Verifikation (OK für öffentliche NASA/JPL-API)
_SSL_CTX = ssl.create_default_context()
_SSL_CTX.check_hostname = False
_SSL_CTX.verify_mode    = ssl.CERT_NONE

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(RES, "OT_meteor_htm.txt")
os.makedirs(RES, exist_ok=True)

# ── HTM-Konstanten ────────────────────────────────────────────────────────────
THETA0  = 58.65
TOL_ANG = 4.0   # °

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
    """J2000 RA/Dec → Galaktisch l/b"""
    ra  = math.radians(ra)
    dc  = math.radians(dec)
    RN  = math.radians(192.85948)
    DN  = math.radians(27.12825)
    b   = math.asin(max(-1., min(1.,
          math.sin(dc)*math.sin(DN) + math.cos(dc)*math.cos(DN)*math.cos(ra-RN))))
    y   = math.cos(dc)*math.sin(ra-RN)
    x   = math.cos(dc)*math.sin(DN)*math.cos(ra-RN) - math.sin(dc)*math.cos(DN)
    l   = (math.degrees(math.atan2(y, x)) + 122.93192) % 360.0
    return l, math.degrees(b)

def lb_xyz(l, b):
    l, b = math.radians(l), math.radians(b)
    return (math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b))

def ang_sep_xyz(v1, v2):
    d = sum(a*b for a,b in zip(v1, v2))
    return math.degrees(math.acos(max(-1., min(1., d))))

def crossing_density(l, b):
    sv = lb_xyz(l, b)
    nc = 0
    for pl, pb, _ in POLES_DEF:
        pv = lb_xyz(pl, pb)
        theta = ang_sep_xyz(pv, sv)
        for n in [1, 2, 3]:
            if abs(theta - n*THETA0) < TOL_ANG:
                nc += 1
    return nc

def eci_to_radiant(vx, vy, vz):
    """ECI-Geschwindigkeit (km/s) → Radiant-Richtung J2000 (RA,Dec)°
    Radiant = Richtung, AUS der der Feuerball kommt = Anti-Velocity."""
    mag = math.sqrt(vx**2 + vy**2 + vz**2)
    if mag < 1e-9:
        return None
    rx, ry, rz = -vx/mag, -vy/mag, -vz/mag
    dec = math.degrees(math.asin(max(-1., min(1., rz))))
    ra  = math.degrees(math.atan2(ry, rx)) % 360.0
    return ra, dec

def rotate_on_sphere(r_xyz, phi_deg, theta_deg):
    """Drehe Einheitsvektor r_xyz um Winkel phi_deg in Richtung theta_deg
    (azimuthal ab r_xyz).
    Gibt neuen Einheitsvektor zurück."""
    rx, ry, rz = r_xyz
    phi = math.radians(phi_deg)
    # Tangenten-Basis um r_xyz
    if abs(rz) < 0.99:
        ax = math.sqrt(rx**2 + ry**2)
        t1 = (-ry/ax, rx/ax, 0.0)        # parallel zur Äquator-Ebene
    else:
        t1 = (1.0, 0.0, 0.0)
    # t2 = r × t1
    t2 = (ry*t1[2] - rz*t1[1],
           rz*t1[0] - rx*t1[2],
           rx*t1[1] - ry*t1[0])
    ct = math.cos(math.radians(theta_deg))
    st = math.sin(math.radians(theta_deg))
    tx = ct*t1[0] + st*t2[0]
    ty = ct*t1[1] + st*t2[1]
    tz = ct*t1[2] + st*t2[2]
    cp = math.cos(phi); sp = math.sin(phi)
    px = cp*rx + sp*tx
    py = cp*ry + sp*ty
    pz = cp*rz + sp*tz
    n  = math.sqrt(px**2 + py**2 + pz**2)
    return (px/n, py/n, pz/n)

def xyz_to_lb(x, y, z):
    b = math.degrees(math.asin(max(-1., min(1., z))))
    l = math.degrees(math.atan2(y, x)) % 360.0
    return l, b

# ── Hardcoded IAU/IMO Meteorschauer-Liste (etabliert) ────────────────────────
# Format: (Code, Name, RA°, Dec°, ZHR, lambda_peak°)
# Quellen: IMO-Kalender 2026, IAU MDC Katalog, Jenniskens (2006)
IAU_SHOWERS = [
    # ── IMO Hauptschwärme (12) ──
    ("QUA", "Quadrantids",            230.1,  48.5, 120, 283.2),
    ("LYR", "April Lyrids",           271.0,  33.6,  18,  32.0),
    ("ETA", "eta Aquariids",          338.0,  -1.0,  50,  45.5),
    ("SDA", "S delta Aquariids",      340.0, -16.4,  25, 125.0),
    ("CAP", "alpha Capricornids",     307.0, -10.0,   5, 127.0),
    ("PER", "Perseids",                48.2,  57.4, 100, 140.0),
    ("ORI", "Orionids",                95.0,  15.8,  20, 208.4),
    ("STA", "S Taurids",               52.0,  13.0,   5, 200.0),
    ("NTA", "N Taurids",               58.0,  22.2,   5, 222.0),
    ("LEO", "Leonids",                152.0,  21.8,  15, 235.3),
    ("GEM", "Geminids",               112.0,  33.0, 150, 261.0),
    ("URS", "Ursids",                 217.0,  75.5,  10, 270.7),
    # ── Weitere etablierte Schwärme (IAU MDC) ──
    ("KCG", "kappa Cygnids",          286.0,  51.0,   3, 145.0),
    ("DRA", "October Draconids",      262.0,  54.0,  10, 195.4),
    ("SPE", "Sep Epsilon-Perseids",    49.0,  40.0,   5, 167.0),
    ("NDA", "N delta Aquariids",      335.0,  -5.0,   4, 127.0),
    ("ERI", "eta Eridanids",           44.0, -13.0,   3,  41.0),
    ("TAH", "tau Herculids",          228.0,  39.0,   5,  68.0),
    ("JBO", "June Bootids",           219.0,  49.0,   5,  95.0),
    ("PHO", "Phoenicids",              15.0, -48.0,   5, 254.0),
    ("PPU", "phi Puppids",            110.0, -45.0,   5,  45.0),
    ("MON", "Dec Monocerotids",       100.0,   8.0,   2, 261.0),
    ("HYD", "sigma Hydrids",          127.0,   2.0,   7, 259.0),
    ("JLE", "June eps Lyrids",        280.0,  42.0,   2,  78.0),
    ("ZCS", "zeta Cygnids",           289.0,  54.0,   4, 143.0),
    ("AHY", "alpha Hydrids",          127.0, -11.0,   4, 266.0),
    ("COM", "Comae Berenicids",       175.0,  18.0,   5, 265.0),
    ("GDR", "gamma Delphinids",       311.0,  16.0,   3,  93.0),
    ("DAU", "delta Aurigids",          84.0,  52.0,   6, 158.0),
    ("XCB", "chi Capricornids",       300.0, -20.0,   3, 127.0),
    ("DSX", "Daytime Sextantids",     152.0,   0.0,  15, 190.6), # radio/daytime
    ("ARI", "Daytime Arietids",        45.0,  24.0,  54,  77.0), # radio/daytime
]

# ── CNEOS Feuerbälle laden ────────────────────────────────────────────────────
def fetch_cneos_fireballs():
    """Lädt alle CNEOS-Feuerbälle mit ECI-Geschwindigkeitskomponenten."""
    url = ('https://ssd-api.jpl.nasa.gov/fireball.api'
           '?vel-comp=true&req-vel-comp=true')
    print(f"  CNEOS API: {url}", flush=True)
    for attempt in range(3):
        try:
            if attempt > 0:
                print(f"  Retry {attempt}...", flush=True)
                time.sleep(5*attempt)
            req = urllib.request.Request(url,
                headers={'User-Agent': 'Python/HTMFireball 1.0'})
            with urllib.request.urlopen(req, timeout=60, context=_SSL_CTX) as r:
                data = json.loads(r.read().decode('utf-8'))
            fields  = data['fields']
            records = data['data']
            print(f"  {data['count']} Einträge, Felder: {fields}")
            return fields, records
        except Exception as e:
            print(f"  FEHLER: {e}")
    return [], []

# ── Monte-Carlo-Vergleich ─────────────────────────────────────────────────────
def mc_nc_distribution(b_samples, n_mc=20000, seed=42):
    """Zufällige Punkte mit gleicher b-Verteilung → nc-Verteilung."""
    random.seed(seed)
    mc = []
    for _ in range(n_mc):
        b = random.choice(b_samples)
        l = random.uniform(0, 360)
        mc.append(crossing_density(l, b))
    return Counter(mc), n_mc

# ── Pfad-Test ──────────────────────────────────────────────────────────────────
N_THETA   = 36         # Azimut-Richtungen
PHI_STEPS = [15, 30, 45, 60]   # Winkelabstände ab Radiant (°)

def path_nc_stats(radiant_lbs):
    """
    Für jeden Radiant (l,b): Taste 36×4 Punkte entlang möglicher Pfade ab.
    Gibt (total_points, Counter(nc)) zurück.
    """
    all_nc = []
    for l0, b0 in radiant_lbs:
        r_xyz = lb_xyz(l0, b0)
        for theta in range(0, 360, 360//N_THETA):
            for phi in PHI_STEPS:
                p_xyz = rotate_on_sphere(r_xyz, phi, theta)
                pl, pb = xyz_to_lb(*p_xyz)
                all_nc.append(crossing_density(pl, pb))
    return Counter(all_nc), len(all_nc)

def path_mc_stats(b_samples, n_mc=5000, seed=99):
    """MC-Referenz für Pfad-Test: gleiche Winkel ab zufälligen Radianten."""
    random.seed(seed)
    all_nc = []
    sample_bs = [random.choice(b_samples) for _ in range(n_mc)]
    for b0 in sample_bs:
        l0 = random.uniform(0, 360)
        r_xyz = lb_xyz(l0, b0)
        for theta in range(0, 360, 360//N_THETA):
            for phi in PHI_STEPS:
                p_xyz = rotate_on_sphere(r_xyz, phi, theta)
                pl, pb = xyz_to_lb(*p_xyz)
                all_nc.append(crossing_density(pl, pb))
    return Counter(all_nc), len(all_nc)

# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("METEOR / FIREBALL  HTM-Knoten-Test")
print("=" * 65)

# ── Teil 1: IAU/IMO Schauer-Radiante ─────────────────────────────────────────
print("\n[1] IAU/IMO Etablierte Meteorschwärme (Radiant-Test)")
print(f"    N = {len(IAU_SHOWERS)} Schwärme")

shower_lbs = []
shower_info = []
for code, name, ra, dec, zhr, lam in IAU_SHOWERS:
    l, b = radec_to_lb(ra, dec)
    nc   = crossing_density(l, b)
    shower_lbs.append((l, b))
    shower_info.append({'code': code, 'name': name,
                        'ra': ra, 'dec': dec, 'zhr': zhr,
                        'l': l, 'b': b, 'nc': nc})

shower_nc = Counter(s['nc'] for s in shower_info)
shower_bs = [s['b'] for s in shower_info]

mc_shower_cnt, mc_shower_n = mc_nc_distribution(shower_bs, n_mc=20000)

print("  nc-Verteilung Schauer | MC-Erwartung:")
for nc in sorted(set(shower_nc.keys()) | set(mc_shower_cnt.keys())):
    f_obs = 100*shower_nc.get(nc,0)/len(shower_info)
    f_mc  = 100*mc_shower_cnt.get(nc,0)/mc_shower_n
    if f_obs > 0 or f_mc > 0:
        ratio = f_obs/f_mc if f_mc > 0 else float('nan')
        print(f"    nc={nc}: {shower_nc.get(nc,0):3d}/{len(shower_info)} "
              f"({f_obs:5.1f}%)  MC={f_mc:.1f}%  ratio={ratio:.2f}")

# Welche Schwärme sitzen bei nc>=2?
print("\n  Schwärme bei nc>=2:")
for s in sorted(shower_info, key=lambda x: -x['nc']):
    if s['nc'] >= 2:
        print(f"    {s['code']:4s} | {s['name']:30s} | "
              f"l={s['l']:5.1f}° b={s['b']:+5.1f}° | nc={s['nc']} | ZHR≈{s['zhr']}")

# ZHR-gewichteter Test
zhr_nc = {}
zhr_all = {}
for s in shower_info:
    # Verwende log(ZHR+1) als Gewicht
    w = math.log(max(1, s['zhr']) + 1)
    for nc in range(s['nc']+1):
        pass
    zhr_nc[s['nc']] = zhr_nc.get(s['nc'], 0) + w
    zhr_all[0] = zhr_all.get(0, 0) + w
total_zhr = sum(zhr_nc.values())
print("\n  ZHR-gewichtete nc-Verteilung (hohe ZHR → größeres Gewicht):")
for nc in sorted(zhr_nc.keys()):
    print(f"    nc={nc}: {100*zhr_nc[nc]/total_zhr:.1f}%")

# ── Teil 2: CNEOS Feuerball-Radiante ──────────────────────────────────────────
print("\n[2] CNEOS Feuerbälle (ECI-Velocity → Radiant J2000)")

fields, records = fetch_cneos_fireballs()

fireballs = []
if fields and records:
    idx = {k: fields.index(k) for k in fields if k in fields}
    for row in records:
        try:
            vx = float(row[fields.index('vx')])
            vy = float(row[fields.index('vy')])
            vz = float(row[fields.index('vz')])
        except (ValueError, IndexError):
            continue
        rd = eci_to_radiant(vx, vy, vz)
        if rd is None:
            continue
        ra_r, dec_r = rd
        l, b = radec_to_lb(ra_r, dec_r)
        nc = crossing_density(l, b)
        # Anti-Radiant (Richtung, in die Meteore fliegen)
        l_anti = (l + 180) % 360
        b_anti = -b
        nc_anti = crossing_density(l_anti, b_anti)
        date = row[fields.index('date')] if 'date' in fields else '?'
        vel  = float(row[fields.index('vel')]) if 'vel' in fields else 0.
        e    = float(row[fields.index('energy')]) if 'energy' in fields else 0.
        fireballs.append({'date': date, 'vx': vx, 'vy': vy, 'vz': vz,
                          'vel': vel, 'energy': e,
                          'ra': ra_r, 'dec': dec_r,
                          'l': l, 'b': b, 'nc': nc,
                          'nc_anti': nc_anti})

print(f"  {len(fireballs)} Feuerbälle mit Radiant-Daten")

if fireballs:
    fb_nc  = Counter(f['nc'] for f in fireballs)
    fb_bs  = [f['b'] for f in fireballs]
    anti_nc = Counter(f['nc_anti'] for f in fireballs)

    mc_fb_cnt, mc_fb_n = mc_nc_distribution(fb_bs, n_mc=20000, seed=43)

    print("  Radiant nc-Verteilung vs MC:")
    for nc in sorted(set(fb_nc.keys()) | set(mc_fb_cnt.keys())):
        f_obs = 100*fb_nc.get(nc,0)/len(fireballs)
        f_mc  = 100*mc_fb_cnt.get(nc,0)/mc_fb_n
        if f_obs > 0 or f_mc > 0:
            ratio = f_obs/f_mc if f_mc > 0 else float('nan')
            print(f"    nc={nc}: {fb_nc.get(nc,0):3d}/{len(fireballs)} "
                  f"({f_obs:5.1f}%)  MC={f_mc:.1f}%  ratio={ratio:.2f}  "
                  + ("★" if ratio > 1.25 and f_obs > 5 else ""))

    print("  Anti-Radiant nc-Verteilung vs MC (Richtung, in die Meteore fliegen):")
    mc_anti_cnt, mc_anti_n = mc_nc_distribution(fb_bs, n_mc=20000, seed=44)
    for nc in sorted(set(anti_nc.keys()) | set(mc_anti_cnt.keys())):
        f_obs = 100*anti_nc.get(nc,0)/len(fireballs)
        f_mc  = 100*mc_anti_cnt.get(nc,0)/mc_anti_n
        if f_obs > 0 or f_mc > 0:
            ratio = f_obs/f_mc if f_mc > 0 else float('nan')
            print(f"    nc={nc}: {anti_nc.get(nc,0):3d}/{len(fireballs)} "
                  f"({f_obs:5.1f}%)  MC={f_mc:.1f}%  ratio={ratio:.2f}  "
                  + ("★" if ratio > 1.25 and f_obs > 5 else ""))

    # Energie-gewichteter Test (log-Energie-Gewicht)
    e_nc_sum = {}
    e_total  = 0
    for f in fireballs:
        w = math.log10(max(0.1, f['energy']) + 1)
        e_nc_sum[f['nc']] = e_nc_sum.get(f['nc'], 0) + w
        e_total += w
    print("  Energie-gewichtete Radiant nc-Verteilung:")
    for nc in sorted(e_nc_sum.keys()):
        print(f"    nc={nc}: {100*e_nc_sum[nc]/e_total:.1f}%")

    # Herausragende Feuerbälle bei nc>=2
    print("\n  Feuerbälle bei Radiant nc>=2 (sortiert nach Energie):")
    high_nc_fb = [f for f in fireballs if f['nc'] >= 2]
    high_nc_fb.sort(key=lambda x: -x['energy'])
    for f in high_nc_fb[:20]:
        print(f"    {f['date']}  E={f['energy']:7.1f}×10¹⁰J  "
              f"l={f['l']:5.1f}° b={f['b']:+5.1f}°  "
              f"nc={f['nc']} nc_anti={f['nc_anti']}")

    # ── Teil 3: Pfad-Test ───────────────────────────────────────────────────
    print("\n[3] Pfad-Test (Punkte entlang Großkreisbögen ab Radiant)")
    print(f"    {N_THETA} Azimute × {len(PHI_STEPS)} Abstände ({PHI_STEPS}°)")

    # Nutze alle Feuerball-Radiante
    fb_radiant_lbs = [(f['l'], f['b']) for f in fireballs]

    path_cnt, path_n = path_nc_stats(fb_radiant_lbs)
    mc_path_cnt, mc_path_n = path_mc_stats(fb_bs, n_mc=len(fireballs), seed=55)

    print("  Pfad-Punkte nc-Verteilung vs MC:")
    for nc in sorted(set(path_cnt.keys()) | set(mc_path_cnt.keys())):
        f_obs = 100*path_cnt.get(nc,0)/path_n
        f_mc  = 100*mc_path_cnt.get(nc,0)/mc_path_n
        if f_obs > 0 or f_mc > 0:
            ratio = f_obs/f_mc if f_mc > 0 else float('nan')
            print(f"    nc={nc}: {path_cnt.get(nc,0):6d}/{path_n} "
                  f"({f_obs:.2f}%)  MC={f_mc:.2f}%  ratio={ratio:.3f}  "
                  + ("★" if ratio > 1.10 and f_obs > 2 else ""))

    # Separiere nach phi-Schritt (kurze vs. lange Pfade)
    print("\n  Pfad-Test separiert nach Angular-Abstand φ:")
    for phi in PHI_STEPS:
        phi_nc = []
        for l0, b0 in fb_radiant_lbs:
            r_xyz = lb_xyz(l0, b0)
            for theta in range(0, 360, 360//N_THETA):
                p_xyz = rotate_on_sphere(r_xyz, phi, theta)
                pl, pb = xyz_to_lb(*p_xyz)
                phi_nc.append(crossing_density(pl, pb))
        cnt_phi = Counter(phi_nc)
        frac_nc1 = 100*sum(v for k,v in cnt_phi.items() if k>=1)/len(phi_nc)
        frac_nc2 = 100*sum(v for k,v in cnt_phi.items() if k>=2)/len(phi_nc)

        # MC mit gleicher b-Verteilung
        random.seed(77 + phi)
        mc_phi_nc = []
        for _ in range(len(fb_radiant_lbs)*N_THETA):
            b_mc = random.choice(fb_bs)
            l_mc = random.uniform(0, 360)
            mc_phi_nc.append(crossing_density(l_mc, b_mc))
        mc_cnt_phi = Counter(mc_phi_nc)
        mc_frac1 = 100*sum(v for k,v in mc_cnt_phi.items() if k>=1)/len(mc_phi_nc)
        mc_frac2 = 100*sum(v for k,v in mc_cnt_phi.items() if k>=2)/len(mc_phi_nc)
        print(f"    φ={phi:2d}°:  nc≥1: {frac_nc1:.2f}% (MC={mc_frac1:.2f}%)"
              f"  nc≥2: {frac_nc2:.2f}% (MC={mc_frac2:.2f}%)"
              f"  Ratio-nc≥2: {frac_nc2/mc_frac2:.3f}" if mc_frac2 > 0
              else f"    φ={phi:2d}°: nc≥1: {frac_nc1:.2f}% (MC={mc_frac1:.2f}%)")

    # ── Schauer Pfad-Test ────────────────────────────────────────────────────
    print("\n  Pfad-Test für Meteor-Schauer-Radiante:")
    sh_path_cnt, sh_path_n = path_nc_stats(shower_lbs)
    mc_sh_path, mc_sh_n    = path_mc_stats(shower_bs, n_mc=len(shower_lbs), seed=66)
    for nc in sorted(set(sh_path_cnt.keys()) | set(mc_sh_path.keys())):
        f_obs = 100*sh_path_cnt.get(nc,0)/sh_path_n
        f_mc  = 100*mc_sh_path.get(nc,0)/mc_sh_n
        if f_obs > 0 or f_mc > 0:
            ratio = f_obs/f_mc if f_mc > 0 else float('nan')
            print(f"    nc={nc}: {f_obs:.2f}%  MC={f_mc:.2f}%  ratio={ratio:.3f}  "
                  + ("★" if ratio > 1.10 and f_obs > 2 else ""))

# ── Zusammenfassung ────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("ZUSAMMENFASSUNG")
print("=" * 65)

def fmt_sep(obs_pct, mc_pct, label):
    ratio = obs_pct/mc_pct if mc_pct > 0 else float('nan')
    flag  = " ★ SIGNAL" if ratio > 1.2 else " (erwartet)" if ratio < 1.05 else ""
    print(f"  {label}: {obs_pct:.1f}% vs MC {mc_pct:.1f}%  → {ratio:.2f}×{flag}")

# Schauer nc>=2
s_nc2 = 100*sum(v for k,v in shower_nc.items() if k>=2)/len(shower_info)
m_nc2 = 100*sum(v for k,v in mc_shower_cnt.items() if k>=2)/mc_shower_n
fmt_sep(s_nc2, m_nc2, "Schauer nc≥2 (von Radiante)")

# Schauer nc>=1
s_nc1 = 100*sum(v for k,v in shower_nc.items() if k>=1)/len(shower_info)
m_nc1 = 100*sum(v for k,v in mc_shower_cnt.items() if k>=1)/mc_shower_n
fmt_sep(s_nc1, m_nc1, "Schauer nc≥1 (von Radiante)")

if fireballs:
    # Feuerbälle nc>=2
    f_nc2 = 100*sum(v for k,v in fb_nc.items() if k>=2)/len(fireballs)
    fm_nc2= 100*sum(v for k,v in mc_fb_cnt.items() if k>=2)/mc_fb_n
    fmt_sep(f_nc2, fm_nc2, "Feuerball-Radiant nc≥2")

    f_nc1 = 100*sum(v for k,v in fb_nc.items() if k>=1)/len(fireballs)
    fm_nc1= 100*sum(v for k,v in mc_fb_cnt.items() if k>=1)/mc_fb_n
    fmt_sep(f_nc1, fm_nc1, "Feuerball-Radiant nc≥1")

    fa_nc2 = 100*sum(v for k,v in anti_nc.items() if k>=2)/len(fireballs)
    fmt_sep(fa_nc2, fm_nc2, "Feuerball Anti-Radiant nc≥2")

# ── Schreibe Ausgabe-Datei ─────────────────────────────────────────────────────
print(f"\nSchreibe Ergebnisse → {OUT}")

import io
buf = io.StringIO()
sys.stdout = buf
# ─ Wiederhole Hauptausgabe in Datei ─ (vereinfacht)
print("=" * 65)
print("METEOR / FIREBALL  HTM-Knoten-Test")
print("=" * 65)

print(f"\nDatumstempel: {time.strftime('%Y-%m-%d %H:%M')}")
print(f"HTM: θ₀={THETA0}°, TOL_ANG={TOL_ANG}°")
print(f"Poles: D(305,25), A(125,-25), S1(35,25), S2(215,-25), S3(215,25), S4(35,-25)")

print(f"\n─── TEST 1: IAU/IMO Meteor-Schauer-Radiante ─────────────────")
print(f"N = {len(shower_info)} etablierte Schwärme")
print(f"\nSchauer-Code | Name                         | l° | b° | nc | ZHR")
for s in sorted(shower_info, key=lambda x: -x['nc']):
    print(f"  {s['code']:4s} | {s['name']:30s} | {s['l']:5.1f} | {s['b']:+5.1f} | "
          f"{s['nc']} | {s['zhr']}")

print(f"\nnc-Verteilung:")
print(f"  {'nc':>4}  {'Schauer':>8}  {'Frac%':>7}  {'MC%':>6}  {'Ratio':>6}")
all_nc_vals = sorted(set(shower_nc.keys()) | set(mc_shower_cnt.keys()))
for nc in all_nc_vals:
    f_obs = 100*shower_nc.get(nc,0)/len(shower_info)
    f_mc  = 100*mc_shower_cnt.get(nc,0)/mc_shower_n
    ratio = f_obs/f_mc if f_mc > 0 else float('nan')
    print(f"  {nc:4d}  {shower_nc.get(nc,0):8d}  {f_obs:7.2f}  {f_mc:6.2f}  {ratio:6.3f}")

print(f"\n─── TEST 2: CNEOS Feuerball-Radiante ─────────────────────────")
print(f"N = {len(fireballs)} Feuerbälle mit ECI-Velocity (req-vel-comp=true)")
if fireballs:
    print(f"\nRadiant nc-Verteilung:")
    print(f"  {'nc':>4}  {'Feuerball':>9}  {'Frac%':>7}  {'MC%':>6}  {'Ratio':>6}")
    for nc in sorted(set(fb_nc.keys()) | set(mc_fb_cnt.keys())):
        f_obs = 100*fb_nc.get(nc,0)/len(fireballs)
        f_mc  = 100*mc_fb_cnt.get(nc,0)/mc_fb_n
        ratio = f_obs/f_mc if f_mc > 0 else float('nan')
        print(f"  {nc:4d}  {fb_nc.get(nc,0):9d}  {f_obs:7.2f}  {f_mc:6.2f}  {ratio:6.3f}"
              + ("  ★" if ratio > 1.2 and f_obs > 5 else ""))

    print(f"\nAnti-Radiant nc-Verteilung:")
    for nc in sorted(set(anti_nc.keys()) | set(mc_anti_cnt.keys())):
        f_obs = 100*anti_nc.get(nc,0)/len(fireballs)
        f_mc  = 100*mc_anti_cnt.get(nc,0)/mc_anti_n
        ratio = f_obs/f_mc if f_mc > 0 else float('nan')
        print(f"  {nc:4d}  {anti_nc.get(nc,0):9d}  {f_obs:7.2f}  {f_mc:6.2f}  {ratio:6.3f}"
              + ("  ★" if ratio > 1.2 and f_obs > 5 else ""))

    print(f"\nFeuerbälle bei Radiant nc≥2 (sortiert nach Energie):")
    for f in sorted([x for x in fireballs if x['nc']>=2], key=lambda x:-x['energy'])[:30]:
        print(f"  {f['date']}  E={f['energy']:8.1f}×10¹⁰J  "
              f"vel={f['vel']:.1f}km/s  "
              f"Radiant: l={f['l']:5.1f}° b={f['b']:+5.1f}°  "
              f"nc={f['nc']} anti-nc={f['nc_anti']}")

print(f"\n─── TEST 3: Pfad-Test ─────────────────────────────────────────")
print(f"Stichprobenpunkte entlang Großkreisbögen ab Radiant")
print(f"({N_THETA} Azimute × {len(PHI_STEPS)} φ-Abstände = {N_THETA*len(PHI_STEPS)} Punkte/Radiant)")
if fireballs:
    print(f"\nAlle φ kombiniert:")
    total_path_n = path_n
    for nc in sorted(set(path_cnt.keys()) | set(mc_path_cnt.keys())):
        f_obs = 100*path_cnt.get(nc,0)/path_n
        f_mc  = 100*mc_path_cnt.get(nc,0)/mc_path_n
        ratio = f_obs/f_mc if f_mc > 0 else float('nan')
        print(f"  nc={nc}: {f_obs:.3f}%  MC={f_mc:.3f}%  ratio={ratio:.4f}")

print(f"\n─── ZUSAMMENFASSUNG ──────────────────────────────────────────")
print(f"  Schauer nc≥1:  {s_nc1:.1f}%  MC={m_nc1:.1f}%  → {s_nc1/m_nc1:.2f}×")
print(f"  Schauer nc≥2:  {s_nc2:.1f}%  MC={m_nc2:.1f}%  → {s_nc2/m_nc2:.2f}×")
if fireballs:
    print(f"  Feuerball nc≥1:{f_nc1:.1f}%  MC={fm_nc1:.1f}%  → {f_nc1/fm_nc1:.2f}×")
    print(f"  Feuerball nc≥2:{f_nc2:.1f}%  MC={fm_nc2:.1f}%  → {f_nc2/fm_nc2:.2f}×")
    print(f"  Anti-Rad. nc≥2:{fa_nc2:.1f}%  MC={fm_nc2:.1f}%  → {fa_nc2/fm_nc2:.2f}×")

sys.stdout = sys.__stdout__
txt = buf.getvalue()
with open(OUT, 'w', encoding='utf-8') as fh:
    fh.write(txt)

print(txt[:1000])
print(f"\n── Vollständige Ausgabe in: {OUT}")
print("FERTIG.")
