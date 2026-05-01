"""
METAGEOMETRA OT-6 Extended — Dualitätssphäre & Skalierende Entfernung
=======================================================================
V19.1: Full SMBH shell analysis incorporating:
  1) Erweiterter Katalog (van den Bosch 2016 + McConnell 2013)
  2) Physikalische Abstände via SIMBAD Radialgeschwindigkeiten (skalierende Entfernung)
  3) Dualitätssphäre: D-Pol UND Anti-D-Pol gemeinsam getestet
  4) Physikalischer Schalenabstand in Mpc (nicht nur Grad)
  5) Anderson-Darling + Rayleigh Statistik

Metageometra-Konstanten:
  θ₀ = 58.65°  — fundamentales Winkelquant
  D-Pol: l=305°, b=+25°  — Dualitätspol
  A-Pol: l=125°, b=-25°  — Antipodaler Dual-Pol (Dualitätssphäre)
  r_s = 4227.9 Mpc        — SRM-Skalenradius (skalierende Entfernung)
  H₀ = 70 km/s/Mpc        — Hubble-Konstante

Author: Kevin Hannemann (via Metageometra V19 automated pipeline)
"""

import os, math, csv, time, urllib.request, urllib.parse, io

# ── Konstanten ────────────────────────────────────────────────────────────────
H0          = 70.0          # km/s/Mpc
C_KMS       = 2.998e5       # km/s
D_POLE_L    = 305.0         # galaktische Länge D-Pol
D_POLE_B    =  25.0         # galaktische Breite D-Pol
A_POLE_L    = 125.0         # Anti-D-Pol (Dualitätssphäre)
A_POLE_B    = -25.0
THETA_0     =  58.65        # Grad
R_S_MPC     = 4227.9        # Mpc — SRM-Skalenradius = c²/(2π·a₀)
SHELL_TOL   =   5.0         # Grad — Schalentoleranz
OUT = os.path.join(os.path.dirname(__file__), "results", "catalogs")
os.makedirs(OUT, exist_ok=True)

# ── Bekannte Literaturdistanzen (Mpc) für nahe Objekte ───────────────────────
# Für d < 30 Mpc ist cz/H₀ wegen Eigenbewegungen unzuverlässig
LIT_DIST_MPC = {
    "Sgr A*"    : 0.00822,   # 8.22 kpc (EHT 2022)
    "M 31"      : 0.785,     # Andromeda (McConnell 2013)
    "M 32"      : 0.785,     # M32 = NGC 221 (Andromeda-Gruppe)
    "NGC 404"   : 3.06,      # TRGB (Seth 2010)
    "NGC 5128"  : 3.78,      # Cen A (Harris 2010)
    "Circinus"  : 4.00,      # ESO 097-G013 (Freeman 1977)
    "M 81"      : 3.63,      # Bode's Galaxy (Freedman 1994)
    "M 77"      : 14.4,      # NGC 1068 (NED)
    "M 106"     : 7.60,      # NGC 4258 Maser (Herrnstein 1999)
    "NGC 3115"  : 9.68,      # SBF (Tonry 2001)
    "M 104"     : 9.55,      # Sombrero (Ford 1996)
    "NGC 3379"  : 10.57,     # Leo-Gruppe (Tonry 2001)
    "NGC 3384"  : 10.57,     # Leo-Gruppe
    "NGC 3377"  : 10.90,     # Leo-Gruppe
    "NGC 3489"  : 11.70,     # Leo-Gruppe
    "NGC 3608"  : 22.30,
    "NGC 4339"  : 16.0,
    "NGC 4342"  : 16.0,
    "NGC 4350"  : 16.8,
    "M 84"      : 17.93,     # NGC 4374 (Virgo)
    "M 85"      : 17.0,      # NGC 4382 (Virgo)
    "M 49"      : 14.90,     # NGC 4472 (Virgo)
    "M 87"      : 16.40,     # NGC 4486 (EHT 2019)
    "M 89"      : 15.90,     # NGC 4552 (Virgo)
    "M 59"      : 16.02,     # NGC 4621 (Virgo)
    "M 60"      : 17.56,     # NGC 4649 (Virgo)
    "NGC 4261"  : 31.60,
    "NGC 4473"  : 15.25,
    "NGC 4459"  : 16.00,
    "NGC 4486B" : 16.40,
    "NGC 4564"  : 15.95,
    "NGC 4596"  : 16.50,
    "NGC 4762"  : 22.60,
    "NGC 4697"  : 12.54,
    "NGC 5128"  : 3.78,
    "NGC 5576"  : 26.00,
    "NGC 5813"  : 29.90,
    "NGC 5845"  : 28.70,
    "NGC 5846"  : 26.30,
    "IC 1459"   : 29.20,
}

# ── Erweiterte Objektliste: van den Bosch (2016) Zusätze ─────────────────────
# Objekte die NICHT in McConnell & Ma (2013) sind
VDBOSCH_EXTRA = [
    # Name              logM  sigma  Quelle
    ("NGC 307",         8.31, 185,  "vdBosch2016"),
    ("NGC 524",         8.92, 237,  "vdBosch2016"),
    ("NGC 1271",        9.43, 257,  "vdBosch2016"),
    ("NGC 1281",        8.17, 190,  "vdBosch2016"),
    ("NGC 2748",        7.69,  90,  "vdBosch2016"),  # upd.
    ("NGC 3585",        8.48, 213,  "vdBosch2016"),
    ("NGC 4278",        8.62, 258,  "vdBosch2016"),
    ("NGC 4371",        7.49, 117,  "vdBosch2016"),
    ("NGC 4526",        8.65, 252,  "vdBosch2016"),
    ("NGC 4552",        8.67, 260,  "vdBosch2016"),  # = M 89 alternative
    ("NGC 6861",        9.31, 390,  "vdBosch2016"),
    ("NGC 7049",        9.10, 251,  "vdBosch2016"),
    ("NGC 7619",        9.40, 324,  "vdBosch2016"),
    ("NGC 7626",        8.93, 272,  "vdBosch2016"),
    ("NGC 4435",        7.80, 160,  "vdBosch2016"),
    ("NGC 4564",        7.90, 162,  "vdBosch2016"),  # upd.
    ("NGC 4649",        9.67, 385,  "vdBosch2016"),  # = M 60
    ("NGC 4621",        8.55, 230,  "vdBosch2016"),  # = M 59
    ("UGC 3789",        7.04,  68,  "vdBosch2016"),  # Maser
    ("NGC 2273",        6.75, 144,  "vdBosch2016"),  # Maser
    ("NGC 4388",        6.93, 107,  "vdBosch2016"),
    ("NGC 1275",        8.89, 270,  "vdBosch2016"),  # Perseus A
    ("NGC 1600",       10.23, 348,  "vdBosch2016"),  # massive SMBH
    ("NGC 4889",       10.32, 347,  "vdBosch2016"),  # Coma cluster
]

# ── Koordinaten-Funktionen ────────────────────────────────────────────────────
def eq_to_gal(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    ra_ngp = math.radians(192.859508)
    dec_ngp = math.radians(27.128336)
    l_ncp = math.radians(122.932)
    sin_b = math.sin(dec)*math.sin(dec_ngp) + math.cos(dec)*math.cos(dec_ngp)*math.cos(ra-ra_ngp)
    b = math.asin(max(-1., min(1., sin_b)))
    x = math.cos(dec)*math.sin(ra-ra_ngp)
    y = math.sin(dec)*math.cos(dec_ngp) - math.cos(dec)*math.sin(dec_ngp)*math.cos(ra-ra_ngp)
    l = (math.degrees(l_ncp - math.atan2(x, y))) % 360
    return l, math.degrees(b)

def gcd_deg(l1, b1, l2, b2):
    l1r, b1r, l2r, b2r = map(math.radians, [l1, b1, l2, b2])
    cos_c = math.sin(b1r)*math.sin(b2r) + math.cos(b1r)*math.cos(b2r)*math.cos(l1r-l2r)
    return math.degrees(math.acos(max(-1., min(1., cos_c))))

def theta_from_pole(ra, dec, pl, pb):
    l, b = eq_to_gal(ra, dec)
    return gcd_deg(pl, pb, l, b)

def nearest_shell(theta):
    best_n, best_d = 1, 999.
    for n in range(1, 7):
        d = abs(theta - THETA_0*n)
        if d < best_d: best_d, best_n = d, n
    return best_n, best_d

# Dual-Pol: teste beide Pole
def dual_shell_test(theta_d, theta_a):
    """
    Dualitätssphäre: Hit wenn ENTWEDER D-Pol ODER A-Pol innerhalb SHELL_TOL.
    Gibt zurück: (hit_bool, pole, n, delta)
    """
    n_d, dd = nearest_shell(theta_d)
    n_a, da = nearest_shell(theta_a)
    if dd <= da:
        hit = dd <= SHELL_TOL
        return hit, "D", n_d, dd
    else:
        hit = da <= SHELL_TOL
        return hit, "A", n_a, da

# ── SIMBAD Einzelabfrage: Koordinaten + Radialgeschwindigkeit ─────────────────
SIMBAD_ALIASES = {
    "Circinus"   : "ESO 097-G013",
    "NGC 4552"   : "M 89",           # Duplikat-Vermeidung
    "NGC 4649"   : "M 60",
    "NGC 4621"   : "M 59",
}
FALLBACK_COORDS = {
    "Circinus": (213.2917, -65.3392),
}

def simbad_one_rv(name):
    """SIMBAD-Abfrage: liefert (ra, dec, rv_kms) oder None.
    rv_kms ist None wenn nicht bekannt."""
    query_name = SIMBAD_ALIASES.get(name, name)
    script = (
        'output console=off script=off\n'
        'format object "%COO(d6;A)|%COO(d6;D)|%RV(V)"\n'
        f'query id {query_name}\n'
    )
    url = "http://simbad.cds.unistra.fr/simbad/sim-script"
    data = urllib.parse.urlencode({"script": script}).encode()
    req = urllib.request.Request(url, data=data, headers={
        "User-Agent": "Python/metageometra-ot6ext",
        "Content-Type": "application/x-www-form-urlencoded"
    })
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            raw = r.read().decode("utf-8", "replace")
    except Exception:
        return None
    for line in raw.split("\n"):
        line = line.strip()
        if "|" in line and not line.startswith(":"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2:
                try:
                    ra  = float(parts[0])
                    dec = float(parts[1])
                    rv  = float(parts[2]) if len(parts) >= 3 and parts[2] not in ("", "~", "---") else None
                    return ra, dec, rv
                except (ValueError, IndexError):
                    pass
    return None

def dist_mpc(name, rv_kms):
    """Physikalische Distanz in Mpc. Priorität: Literatur > SIMBAD-RV > None."""
    if name in LIT_DIST_MPC:
        return LIT_DIST_MPC[name], "lit"
    if rv_kms is not None and rv_kms > 200:   # < 200 km/s: Eigenbewegung dominiert
        return rv_kms / H0, "rv"
    return None, "unknown"

# ── Katalog laden / aufbauen ──────────────────────────────────────────────────
def load_base_catalog():
    """Liest smbh_catalog_combined.csv und gibt Liste von Dicts zurück."""
    path = os.path.join(OUT, "smbh_catalog_combined.csv")
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows

# ── Hauptprogramm ─────────────────────────────────────────────────────────────
def main():
    print("\n" + "="*70)
    print("  METAGEOMETRA OT-6 Extended")
    print("  Dualitätssphäre + Skalierende Entfernung")
    print("="*70)

    # Schlüsselerkenntnis: A-Pol = exakter Antipode des D-Pols
    # → für jeden Punkt P gilt: θ_A = 180° - θ_D
    # → keine extra SIMBAD-Abfragen nötig für Basis-Objekte!

    def norm(n): return n.strip().upper().replace(" ", "")

    # 1. Basis-Katalog laden + θ_A = 180° - θ_D berechnen
    base = load_base_catalog()
    print(f"\n  Basis-Katalog: {len(base)} Objekte aus smbh_catalog_combined.csv")
    for row in base:
        td = float(row["theta_dpole"])
        ta = 180.0 - td                    # Antipodensymmetrie
        row["theta_apole"] = f"{ta:.3f}"
        d_mpc, d_src = dist_mpc(row["Name"], None)
        row["dist_mpc"] = f"{d_mpc:.2f}" if d_mpc else ""
        row["dist_src"] = d_src
        row["rv_kms"] = ""
    base_norm = {norm(r["Name"]) for r in base}

    # 2. Neue van den Bosch Objekte (nur echte Neuzugänge)
    extra_new = [(n,m,s,src) for n,m,s,src in VDBOSCH_EXTRA
                 if norm(n) not in base_norm]
    print(f"  Neue vdBosch-2016-Objekte: {len(extra_new)} (SIMBAD-Abfrage...)")

    # 3. SIMBAD nur für die ~22 neuen Objekte
    new_rows = []
    not_found = []
    for i, (name, logM, sigma, src) in enumerate(extra_new):
        result = simbad_one_rv(name)
        if result is None and name in FALLBACK_COORDS:
            ra, dec = FALLBACK_COORDS[name]; rv = None; result = (ra, dec, rv)
        if result:
            ra, dec, rv = result
            th_d = theta_from_pole(ra, dec, D_POLE_L, D_POLE_B)
            th_a = 180.0 - th_d            # Antipodensymmetrie
            n_d, dd = nearest_shell(th_d)
            d_mpc, d_src = dist_mpc(name, rv)
            new_rows.append({
                "Name": name, "RA_deg": f"{ra:.5f}", "Dec_deg": f"{dec:.5f}",
                "logMbh": logM, "sigma_kms": sigma, "source": src,
                "theta_dpole": f"{th_d:.3f}", "shell_n": n_d,
                "shell_delta": f"{dd:.3f}",
                "hit_5deg": "YES" if dd <= SHELL_TOL else "no",
                "rv_kms": f"{rv:.1f}" if rv else "",
                "dist_mpc": f"{d_mpc:.2f}" if d_mpc else "",
                "dist_src": d_src,
                "theta_apole": f"{th_a:.3f}",
            })
            flag = "YES" if dd <= SHELL_TOL else "   "
            print(f"  [{i+1:2d}/{len(extra_new)}] {flag}  {name:<20} "
                  f"θD={th_d:6.2f}°  θA={th_a:6.2f}°  "
                  f"d={'??' if d_mpc is None else f'{d_mpc:.0f}':>6} Mpc")
        else:
            not_found.append(name)
            print(f"  [{i+1:2d}/{len(extra_new)}] ---  {name} (nicht gefunden)")
        time.sleep(0.3)

    # 4. Alle Zeilen zusammenführen + Duplikate bereinigen
    all_rows = list(base) + new_rows
    seen = {}
    for r in all_rows:
        seen[norm(r["Name"])] = r
    all_rows = list(seen.values())
    print(f"\n  Gesamt nach Merge: {len(all_rows)} Objekte")

    # 6. CSV schreiben (erweitertes Format)
    cols = ["Name","RA_deg","Dec_deg","logMbh","sigma_kms","source",
            "theta_dpole","theta_apole","shell_n","shell_delta","hit_5deg",
            "rv_kms","dist_mpc","dist_src"]
    out_path = os.path.join(OUT, "smbh_extended.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"  Gespeichert: {out_path}")

    # ── Statistiken ──────────────────────────────────────────────────────────
    import scipy.stats as st

    # A) D-Pol Standard (wie bisher)
    hits_d   = [r for r in all_rows if r.get("hit_5deg","no") == "YES"]
    n_total  = len(all_rows)
    frac_exp = 6 * 10 / 180  # 6 Schalen × 10° / 180° = 33.3%

    # B) Dualitätssphäre: D-Pol + A-Pol kombiniert
    # Effektive Schalen in θ_D-Koordinaten:
    # D-Pol: 58.65°, 117.30°, 175.95°
    # A-Pol: 180°-58.65°=121.35°, 180°-117.30°=62.70°, 180°-175.95°=4.05°
    # Zusammengeführte Coverage-Berechnung:
    shell_centers = sorted([
        THETA_0*1,          # D n=1 = 58.65
        THETA_0*2,          # D n=2 = 117.30
        THETA_0*3,          # D n=3 = 175.95
        180 - THETA_0*1,    # A n=1 = 121.35
        180 - THETA_0*2,    # A n=2 = 62.70
        180 - THETA_0*3,    # A n=3 = 4.05
    ])
    # Coverage = Länge der Vereinigung aller [c-5, c+5] ∩ [0,180]
    intervals = [(max(0, c-SHELL_TOL), min(180, c+SHELL_TOL)) for c in shell_centers]
    intervals.sort()
    merged = []
    for lo, hi in intervals:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append([lo, hi])
    dual_coverage = sum(hi-lo for lo,hi in merged)
    frac_dual = dual_coverage / 180.0

    dual_hits = []
    for r in all_rows:
        td = float(r["theta_dpole"])
        ta = float(r["theta_apole"])          # = 180° - td (Antipodensymmetrie)
        _, n_d, dd = nearest_shell(td)[1], *nearest_shell(td)  # unpack
        n_dd, delta_d = nearest_shell(td)
        n_da, delta_a = nearest_shell(ta)
        hit_d = delta_d <= SHELL_TOL
        hit_a = delta_a <= SHELL_TOL
        if hit_d or hit_a:
            dual_hits.append(r)

    # C) Skalierende Entfernung: physikalischer Schalenabstand
    rows_with_dist = [r for r in all_rows if r.get("dist_mpc","")]
    phys_offsets = []
    for r in rows_with_dist:
        d = float(r["dist_mpc"])
        delta_deg = float(r["shell_delta"])
        delta_rad = math.radians(delta_deg)
        phys_mpc = d * math.sin(delta_rad)  # Bogenabstand = d × sin(Δθ)
        phys_offsets.append((r["Name"], d, delta_deg, phys_mpc,
                              r.get("hit_5deg","no") == "YES"))

    # ── Ausgabe ──────────────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("  ERGEBNISSE — DUALITÄTSSPHÄRE & SKALIERENDE ENTFERNUNG")
    print("="*70)

    print(f"\n  A) D-Pol Standard (θ₀ = {THETA_0}°, Toleranz ±{SHELL_TOL}°):")
    print(f"     Treffer: {len(hits_d)}/{n_total} = {len(hits_d)/n_total:.1%}")
    print(f"     Zufallserwartung: {frac_exp:.1%}")
    b_ex = st.binomtest(len(hits_d), n_total, frac_exp, alternative='greater')
    b_av = st.binomtest(len(hits_d), n_total, frac_exp, alternative='less')
    b_2s = st.binomtest(len(hits_d), n_total, frac_exp, alternative='two-sided')
    print(f"     p(Überschuss):    {b_ex.pvalue:.4f}  {'*' if b_ex.pvalue < 0.05 else ''}")
    print(f"     p(Vermeidung):    {b_av.pvalue:.4f}  {'*' if b_av.pvalue < 0.05 else ''}")
    print(f"     p(zweiseitig):    {b_2s.pvalue:.4f}  {'*' if b_2s.pvalue < 0.05 else ''}")

    print(f"\n  B) Dualitätssphäre (D-Pol ∪ A-Pol, Antipodenpol l=125°,b=-25°):")
    print(f"     Effektive Schalenabdeckung: {dual_coverage:.1f}° / 180° = {frac_dual:.1%}")
    print(f"     Effektive Schalenzentren (θ_D): {[f'{c:.1f}' for c in shell_centers]}")
    print(f"     Verchmelzte Intervalle: {merged}")
    print(f"     Treffer (D oder A): {len(dual_hits)}/{n_total} = {len(dual_hits)/n_total:.1%}")
    print(f"     Zufallserwartung (dual): {frac_dual:.1%}")
    db_ex = st.binomtest(len(dual_hits), n_total, frac_dual, alternative='greater')
    db_av = st.binomtest(len(dual_hits), n_total, frac_dual, alternative='less')
    db_2s = st.binomtest(len(dual_hits), n_total, frac_dual, alternative='two-sided')
    print(f"     p(Überschuss):    {db_ex.pvalue:.4f}  {'*' if db_ex.pvalue < 0.05 else ''}")
    print(f"     p(Vermeidung):    {db_av.pvalue:.4f}  {'*' if db_av.pvalue < 0.05 else ''}")
    print(f"     p(zweiseitig):    {db_2s.pvalue:.4f}  {'*' if db_2s.pvalue < 0.05 else ''}")

    # Interpretiere Dualitätssphäre
    if len(dual_hits)/n_total > frac_dual:
        interp_dual = "SMBH bevorzugen Schalenpositionen (Dualitätssphäre)"
    else:
        interp_dual = "SMBH meiden Schalenpositionen (Dualitätssphäre)"
    print(f"\n     >> {interp_dual}")

    print(f"\n  C) Skalierende Entfernung (SRM r_s = {R_S_MPC} Mpc):")
    print(f"     Objekte mit Distanzmessung: {len(rows_with_dist)}/{n_total}")
    if phys_offsets:
        hits_phys  = [x for x in phys_offsets if x[4]]
        miss_phys  = [x for x in phys_offsets if not x[4]]
        mean_hit   = sum(x[3] for x in hits_phys)/len(hits_phys)  if hits_phys  else 0
        mean_miss  = sum(x[3] for x in miss_phys)/len(miss_phys)  if miss_phys  else 0
        mean_all   = sum(x[3] for x in phys_offsets)/len(phys_offsets)
        print(f"\n     Mittl. physik. Schalenabstand d×sin(Δθ):")
        print(f"       Alle Objekte:        {mean_all:8.1f} Mpc")
        print(f"       Treffer (Δθ<5°):     {mean_hit:8.1f} Mpc")
        print(f"       Nicht-Treffer:       {mean_miss:8.1f} Mpc")
        print(f"\n     d_skaliert = d / r_s (SRM-Einheiten):")
        print(f"     {'Name':<22} {'d(Mpc)':>8} {'d/r_s':>7} {'Δθ(°)':>7} {'Δr(Mpc)':>9}  hit")
        for nm, d, dt, dp, hit in sorted(phys_offsets, key=lambda x: x[1]):
            d_scaled = d / R_S_MPC
            print(f"     {nm:<22} {d:>8.1f} {d_scaled:>7.5f} {dt:>7.2f} {dp:>9.1f}  {'HIT' if hit else ''}")

    print(f"\n  Dualitätssphäre Treffer-Tabelle:")
    print(f"  {'Name':<22} {'θ_D':>7}  {'θ_A':>7}  {'hit_D':>5}  {'hit_A':>5}  d(Mpc)")
    for r in sorted(all_rows, key=lambda x: float(x["theta_dpole"])):
        td = float(r["theta_dpole"])
        ta = float(r["theta_apole"])
        _, dd = nearest_shell(td)
        _, da = nearest_shell(ta)
        hit_d = "HIT" if dd <= SHELL_TOL else f"{dd:.1f}°"
        hit_a = "HIT" if da <= SHELL_TOL else f"{da:.1f}°"
        d_str = r.get("dist_mpc", "") or "-"
        print(f"  {r['Name']:<22} {td:>7.2f}  {ta:>7.2f}  {hit_d:>5}  {hit_a:>5}  {d_str}")

    # ── D) MASSENGRADIENTEN-TEST (V19 Torsionsschock-Modell) ─────────────────
    # Vorhersage: SMBHs kondensieren an Kreuzpunkten der Fraktalschalen.
    # Knotenenergie nimmt mit Abstand von der Dualitaetssphaere (D-Pol) ab.
    # => M_BH sollte mit theta_D NEGATIV korrelieren (naeher am D-Pol -> massiver).
    #
    # Entstehungsmodell: Zwei invertiert-rotierende SMBHs kollidieren im Meta-Vakuum
    # -> Torsionsschock -> Rift entsteht -> Fraktalschalen als Stehwellen.
    # SMBHs = eingefrorene Knotenenergie. Ringstrukturen halten den Rift offen
    # und verhindern Big Rip / Kaeltetod durch Torsions-Gegendruck.
    print("\n" + "="*70)
    print("  D) MASSENGRADIENTEN-TEST: M_BH vs Abstand zur Dualitaetssphaere")
    print("="*70)
    print("  Vorhersage (V19 Torsionsschock): M_BH sinkt mit theta_D (D-Pol-Abstand)")
    print("  Entstehung: SMBHs = Knotenenergie der Schockwellen-Kreuzpunkte")

    grad_data = []
    for r in all_rows:
        try:
            log_m = float(r["logMbh"])
            th_d  = float(r["theta_dpole"])
            if log_m > 0 and 0 < th_d < 180:
                grad_data.append((th_d, log_m, r["Name"]))
        except (ValueError, TypeError, KeyError):
            continue

    if len(grad_data) >= 5:
        import scipy.stats as _st2
        thetas  = [x[0] for x in grad_data]
        log_ms  = [x[1] for x in grad_data]
        sp_r, sp_p = _st2.spearmanr(thetas, log_ms)
        pe_r, pe_p = _st2.pearsonr(thetas, log_ms)
        # Lineare Regression: log_M = a * theta + b
        slope, intercept, r_val, p_lin, se = _st2.linregress(thetas, log_ms)

        print(f"\n  Stichprobe: {len(grad_data)} SMBHs mit log(M_BH) und theta_D")
        print(f"  Spearman r = {sp_r:+.4f},  p = {sp_p:.4f}  {'<- signifikant *' if sp_p < 0.05 else ''}")
        print(f"  Pearson  r = {pe_r:+.4f},  p = {pe_p:.4f}  {'<- signifikant *' if pe_p < 0.05 else ''}")
        print(f"  Linreg:   log(M) = {slope:+.5f} * theta_D + {intercept:.3f}")
        print(f"            R^2 = {r_val**2:.4f},  p = {p_lin:.4f}")

        dml_per_10deg = slope * 10
        print(f"\n  Massengradienten: {dml_per_10deg:+.3f} dex pro 10 Grad theta_D")

        if sp_p < 0.05 and sp_r < 0:
            verdict = "BESTAETIGT: SMBHs naeher am D-Pol sind massiver (V19-Vorhersage)"
        elif sp_p < 0.05 and sp_r > 0:
            verdict = "UMGEKEHRT: SMBHs ferner vom D-Pol sind massiver (Modell-Widerspruch)"
        else:
            verdict = f"KEIN SIGNAL (p={sp_p:.3f}) — kein signifikanter Massengradient"
        print(f"\n  BEWERTUNG: {verdict}")

        # Top 5 massereichsten vs naechste zum D-Pol
        sorted_by_mass = sorted(grad_data, key=lambda x: x[1], reverse=True)
        sorted_by_dist = sorted(grad_data, key=lambda x: x[0])
        print(f"\n  Top-5 massereichste SMBHs:")
        for td, lm, nm in sorted_by_mass[:5]:
            print(f"    {nm:<22}  log(M)={lm:.2f}  theta_D={td:.1f}°")
        print(f"\n  Top-5 D-Pol-naechste SMBHs:")
        for td, lm, nm in sorted_by_dist[:5]:
            print(f"    {nm:<22}  log(M)={lm:.2f}  theta_D={td:.1f}°")
    else:
        print(f"  Zu wenige Datenpunkte ({len(grad_data)}) fuer Gradientenanalyse.")

    print("\n" + "="*70)
    print(f"  Ausgabe: {out_path}")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()
