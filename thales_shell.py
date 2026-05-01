"""
Satz des Thales — Vakuumschalen-Radius aus SMBH-Positionen
===========================================================
Geometrie:
  Wir (Erde) sitzen bei Ursprung (0,0).
  D-Pol zeigt in Richtung (l=305, b=25) — die Schock-Achse.
  Die Vakuumschale ist ein Kreis (Querschnitt) mit Durchmesser R entlang dieser Achse.

  Satz des Thales: Jeder Punkt P auf einem Kreis mit Durchmesser von O nach (R,0)
  erfuellt: x^2 + y^2 = R * x
  => R_Thales = d^2 / d_parallel = d / cos(theta_D)

  Falls SMBHs auf der Vakuumschale kondensiert sind, muss fuer alle gelten:
  R_Thales = const = Radius der Vakuumschale in Mpc.
"""
import math, csv, os, sys
import statistics

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
CSV_PATH  = os.path.join(RESULTS, "catalogs", "smbh_extended.csv")
DPOLE_L, DPOLE_B = 305.0, 25.0

def radec_to_lb(ra_deg, dec_deg):
    ra_r = math.radians(ra_deg); dc_r = math.radians(dec_deg)
    RA_NGP = math.radians(192.85948); DEC_NGP = math.radians(27.12825)
    b_r = math.asin(math.sin(dc_r)*math.sin(DEC_NGP) +
                    math.cos(dc_r)*math.cos(DEC_NGP)*math.cos(ra_r-RA_NGP))
    l_r = math.atan2(
        math.cos(dc_r)*math.sin(ra_r-RA_NGP),
        math.cos(dc_r)*math.sin(DEC_NGP)*math.cos(ra_r-RA_NGP) - math.sin(dc_r)*math.cos(DEC_NGP))
    return (math.degrees(l_r) + 122.93192) % 360.0, math.degrees(b_r)

def angle_dpole(l_deg, b_deg):
    l1,b1 = math.radians(l_deg), math.radians(b_deg)
    l2,b2 = math.radians(DPOLE_L), math.radians(DPOLE_B)
    c = math.sin(b1)*math.sin(b2) + math.cos(b1)*math.cos(b2)*math.cos(l1-l2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))

# ── Katalog laden ──────────────────────────────────────────────────────────────
rows = []
with open(CSV_PATH, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        dist_str = r.get('dist_mpc','').strip()
        if not dist_str:
            continue
        try:
            d_mpc = float(dist_str)
            if d_mpc <= 0:
                continue
            ra  = float(r['RA_deg'])
            dc  = float(r['Dec_deg'])
            lm  = float(r['logMbh'])
            l_deg, b_deg = radec_to_lb(ra, dc)
            theta = angle_dpole(l_deg, b_deg)
            theta_r = math.radians(theta)

            # Thales-Zerlegung in Koordinatensystem entlang D-Pol-Achse
            d_parallel = d_mpc * math.cos(theta_r)  # Projektion auf D-Pol-Achse
            d_perp     = d_mpc * math.sin(theta_r)  # senkrecht dazu

            # Satz des Thales: R = d^2 / d_parallel  (nur wenn d_parallel > 0)
            if d_parallel > 0:
                R_thales = (d_mpc ** 2) / d_parallel
            else:
                R_thales = float('nan')

            rows.append({
                'name':      r['Name'],
                'logM':      lm,
                'd_mpc':     d_mpc,
                'theta':     theta,
                'd_par':     d_parallel,
                'd_perp':    d_perp,
                'R_thales':  R_thales,
                'l': l_deg, 'b': b_deg
            })
        except Exception:
            continue

rows.sort(key=lambda r: r['theta'])
print(f"SMBHs mit Distanz: {len(rows)}")
print()

# ── Ausgabe: alle Objekte ──────────────────────────────────────────────────────
print(f"{'Name':<22} {'d[Mpc]':>8} {'theta':>7} {'d_par':>8} {'d_perp':>8} {'R_Thales':>10} {'log(M)':>7}")
print("-" * 80)
for r in rows:
    rt = f"{r['R_thales']:.1f}" if not math.isnan(r['R_thales']) else "  ---"
    print(f"  {r['name']:<20} {r['d_mpc']:>8.2f} {r['theta']:>6.1f}° "
          f"{r['d_par']:>8.2f} {r['d_perp']:>8.2f} {rt:>10} {r['logM']:>7.2f}")

# ── Statistik der R_Thales Werte ───────────────────────────────────────────────
valid = [r for r in rows if not math.isnan(r['R_thales']) and r['R_thales'] > 0]
R_vals = [r['R_thales'] for r in valid]

print()
print("=" * 60)
print("THALES-STATISTIK")
print("=" * 60)
print(f"  N (mit Distanz + theta > 0):  {len(valid)}")
print(f"  Median R_Thales:  {statistics.median(R_vals):.1f} Mpc")
print(f"  Mittel  R_Thales:  {statistics.mean(R_vals):.1f} Mpc")
print(f"  Std-Abw:           {statistics.stdev(R_vals):.1f} Mpc")
print(f"  Min:               {min(R_vals):.1f} Mpc")
print(f"  Max:               {max(R_vals):.1f} Mpc")
print()

# Cluster um Median (±30%)
R_med = statistics.median(R_vals)
cluster = [r for r in valid if abs(r['R_thales'] - R_med)/R_med < 0.30]
print(f"  Cluster um Median ±30%:  {len(cluster)}/{len(valid)} Objekte")
print(f"  Cluster-Median:    {statistics.median(r['R_thales'] for r in cluster):.1f} Mpc")
print()

# ── Die naechsten SMBHs zum D-Pol (kleinste theta) ────────────────────────────
naechste = sorted(valid, key=lambda r: r['theta'])[:8]
print("Naechste SMBHs zur D-Pol-Richtung:")
print(f"{'Name':<22} {'theta':>7} {'d[Mpc]':>8} {'R_Thales':>10} {'log(M)':>7}")
print("-" * 58)
for r in naechste:
    print(f"  {r['name']:<20} {r['theta']:>6.1f}° {r['d_mpc']:>8.2f} {r['R_thales']:>10.1f} {r['logM']:>7.2f}")

# ── Thales-Geometrie erklaert ─────────────────────────────────────────────────
print()
print("=" * 60)
print("GEOMETRISCHE INTERPRETATION")
print("=" * 60)
print(f"""
  Koordinatensystem (Ursprung = Erde):
    x-Achse → D-Pol (l=305°, b=25°)
    y-Achse ⊥ dazu

  SMBH bei Winkel theta und Distanz d:
    x = d*cos(theta)  [entlang D-Pol]
    y = d*sin(theta)  [senkrecht]

  Satz des Thales:
    Kreis mit Durchmesser R (von (0,0) bis (R,0)):
    x^2 + y^2 = R*x
    => R = d^2 / (d*cos(theta)) = d / cos(theta)

  Wenn alle SMBHs auf der Vakuumschale liegen:
    R_Thales muss fuer alle ≈ gleich sein.

  Streuung der R_Thales-Werte zeigt:
    - Kleine Streuung -> SMBHs liegen tatsaechlich auf Schale
    - Grosse Streuung -> SMBHs von verschiedenen Schalenebenen ODER
                         Sekundaerschalen (Hopf-Subharmoniken)
""")

# Interessant: Verhältnis R / a0
A0_HTM = 1.0964e-10  # m/s^2
C_LIGHT = 2.998e8    # m/s
MPC_IN_M = 3.0857e22
# R_s = c^2 / (2pi * a0) in Mpc
R_SRM = (C_LIGHT**2) / (2*math.pi*A0_HTM) / MPC_IN_M
print(f"  Vergleich mit HTM SRM-Skalenradius:")
print(f"    r_s = c^2/(2pi*a0) = {R_SRM:.1f} Mpc")
print(f"    Median R_Thales    = {R_med:.1f} Mpc")
print(f"    Verhältnis:          {R_med/R_SRM:.4f}")
print()

# Nahe D-pol SMBHs: Mittelwert
R_nahe = statistics.mean(r['R_thales'] for r in naechste)
print(f"  R_Thales (8 D-Pol-naechste):  {R_nahe:.1f} Mpc")
print(f"  Verhältnis zu r_s:             {R_nahe/R_SRM:.4f}")

# Ergebnis speichern
out = os.path.join(RESULTS, "thales_shell_radius.txt")
with open(out, 'w', encoding='utf-8') as f:
    f.write("SATZ DES THALES — Vakuumschalen-Radius\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"D-Pol: l={DPOLE_L}, b={DPOLE_B}\n")
    f.write(f"Methode: R_Thales = d / cos(theta_D)\n\n")
    f.write(f"N Objekte:        {len(valid)}\n")
    f.write(f"Median R_Thales:  {R_med:.1f} Mpc\n")
    f.write(f"Std-Abw:          {statistics.stdev(R_vals):.1f} Mpc\n")
    f.write(f"r_s (HTM):        {R_SRM:.1f} Mpc\n")
    f.write(f"R_Thales / r_s:   {R_med/R_SRM:.4f}\n\n")
    f.write("Naechste SMBHs zum D-Pol:\n")
    for r in naechste:
        f.write(f"  {r['name']:<22} theta={r['theta']:.1f}°  d={r['d_mpc']:.1f} Mpc  R={r['R_thales']:.1f} Mpc\n")
print(f"\nGespeichert: {out}")
