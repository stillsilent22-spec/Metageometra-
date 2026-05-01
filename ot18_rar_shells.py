"""
OT-18: KS-Test RAR-Residuen vs Schalenabstand (SPARC)
=======================================================
Aufgabe: Testen ob SPARC-Galaxien NAHE an HTM-Schalen (theta_n = n * 58.65 deg)
systematisch GROESSERE oder KLEINERE RAR-Streuungen (log-Residuen vom McGaugh+2016
MDAR) zeigen als weit entfernte Galaxien.

TORSIONSSCHOCK-HYPOTHESE (V19-Erweiterung):
  Im HTM-Framework ist L = rho_DE / t0 global — im Gleichgewichtszustand.
  Beim Durchgang der kosmischen Ausdehnung durch eine Schalengrenze
  (theta_n = n * theta_0) kann jedoch ein transienter TORSIONSSCHOCK auftreten:
  Der S3-Torsions-Knoten "springt" auf die naechste Windungszahl (Hopf-Index).
  Waehrend dieses Uebergangs ist L kurz erhoeht: L_eff = L_0 * (1 + delta_L/L).
  Galaxien nahe aktueller Schalenraender koennten diesen Shock als erhoehte
  oder reduzierte a0 = c/(2pi*t) "eingefroren" haben (fossiles Signal).

  Konsequenz fuer diesen Test:
  Falls delta_L/L != 0 bei Schalennaehe -> sigma_RAR(nahe) != sigma_RAR(fern)
  Falls delta_L/L == 0 (kein Shock) -> sigma_RAR(nahe) == sigma_RAR(fern) [Null]

  OT-18 ist damit ein ECHTER TEST — kein 'predicted null' a priori,
  sondern ein Test auf Torsionsschock-Amplitude.
  Ein Null-Ergebnis liefert eine OBERE SCHRANKE auf |delta_L/L|.

Methodik:
  1. Lade SPARC-Daten (table1.dat: Namen, Positionen, table2.dat: Rotation curves)
  2. Berechne fuer jede Galaxie a_obs und a_bar (in m/s^2)
  3. Fittte MDAR f(a_bar): a_obs ~ a_bar / (1 - exp(-sqrt(a_bar/g_dag)))
     mit g_dag = 1.2e-10 m/s^2 (festgelegt)
  4. RAR-Streuung pro Galaxie: sigma_RAR = std(log10(a_obs/a_MDAR))
  5. Hole RA/Dec aus VizieR TAP (J/AJ/152/157)
  6. Berechne theta_D (Winkelabstand vom HTM-Dipol-Pol: l=305, b=25)
  7. Delta_theta = min(|theta_D - n*58.65| fuer n=1..6)
  8. Teile in "nahe" (Delta_theta < 10 deg) und "fern" (Delta_theta >= 10 deg)
  9. KS-Test und Mann-Whitney-U Test auf sigma_RAR Verteilungen
 10. Obere Schranke auf |delta_L/L| aus Null-Ergebnis
"""

import os, sys, math, io, time
import numpy as np
from scipy import stats
import urllib.request
import urllib.parse

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

# ── Framework-Konstanten ──────────────────────────────────────────────────────
SHELLS_DEG    = [58.65 * n for n in range(1, 7)]  # [58.65, 117.3, 175.95, 234.6, 293.25, 351.9]
DPOLE_L       = 305.0   # Galaktische Laenge des HTM-Dipol-Pols [Grad]
DPOLE_B       = 25.0    # Galaktische Breite [Grad]
DELTA_SHELL   = 10.0    # Naehe-Schwelle [Grad]
G_DAG         = 1.2e-10 # m/s^2  (MDAR Beschleunigungsparameter)

# Einheiten-Konvertierung
KM_PER_KPC    = 3.0857e16  # km/kpc
M_PER_KPC     = 3.0857e19  # m/kpc

def a_from_V_r(V_kms, r_kpc):
    """Radiale Beschleunigung aus Rotationsgeschwindigkeit."""
    return V_kms**2 * 1e6 / (r_kpc * M_PER_KPC)  # m/s^2

def a_mdar(a_bar, g_dag=G_DAG):
    """MDAR Vorhersage: a_obs = a_bar / (1 - exp(-sqrt(a_bar/g_dag)))."""
    with np.errstate(over='ignore', invalid='ignore'):
        x = np.sqrt(np.where(a_bar > 0, a_bar / g_dag, 0.0))
        denom = 1.0 - np.exp(-x)
        # Limite: a_bar -> 0 → a_obs -> 0; a_bar >> g_dag → a_obs -> a_bar
        result = np.where(a_bar > 1e-20, a_bar / np.where(denom > 1e-10, denom, 1e-10), 0.0)
    return result

# ── SPARC table2.dat laden ─────────────────────────────────────────────────────
TABLE2_CACHE = os.path.join(RESULTS, "sparc_table2.dat")
TABLE2_URL   = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat"

def fetch_file(url, cache_path):
    if os.path.exists(cache_path):
        with open(cache_path, 'r', encoding='ascii', errors='replace') as f:
            return f.read()
    print(f"  Lade {url} ...")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read().decode('ascii', errors='replace')
    with open(cache_path, 'w', encoding='ascii', errors='replace') as f:
        f.write(data)
    return data

print("Lade SPARC Daten...")
t2_text = fetch_file(TABLE2_URL, TABLE2_CACHE)
print(f"  table2.dat: {len(t2_text.splitlines())} Zeilen")

# table2 Spaltenformat (J/AJ/152/157):
# Bytes  1-11:  Name
# 12:    " "
# 13-18: Dist (Mpc)
# 20-25: Rad (kpc)
# 27-32: Vobs (km/s)
# 34-38: e_Vobs
# 40-45: Vgas
# 47-52: Vdisk
# 54-59: Vbulge

galaxies = {}
for line in t2_text.splitlines():
    if len(line) < 59:
        continue
    try:
        name   = line[0:11].strip()
        dist   = float(line[12:18])
        rad    = float(line[19:25])
        vobs   = float(line[26:32])
        evobs  = float(line[33:38])
        vgas   = float(line[39:45])
        vdisk  = float(line[46:52])
        vbulge = float(line[53:59])
    except ValueError:
        continue
    if name not in galaxies:
        galaxies[name] = []
    galaxies[name].append((rad, vobs, evobs, vgas, vdisk, vbulge))

print(f"  {len(galaxies)} Galaxien in table2.dat")

# ── RA/Dec aus SIMBAD TAP (Batch-Query fuer alle SPARC-Namen) ─────────────────
POS_CACHE = os.path.join(RESULTS, "sparc_positions.csv")

# Hole alle Galaxy-Namen aus table2 (gleiche wie table1)
gal_names = list(galaxies.keys())

def fetch_simbad_positions(names):
    """Hole RA/Dec fuer eine Liste von Namen via SIMBAD TAP."""
    # Baue IN-Clause (max 200 Namen pro Batch)
    in_vals = ', '.join("'" + n.replace("'", "''") + "'" for n in names)
    query = (
        "SELECT i.id, b.ra, b.dec "
        "FROM ident i "
        "JOIN basic b ON b.oid = i.oidref "
        f"WHERE i.id IN ({in_vals})"
    )
    simbad_tap = "https://simbad.u-strasbg.fr/simbad/sim-tap/sync"
    url = (simbad_tap + "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv"
           f"&QUERY={urllib.parse.quote(query)}")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode('utf-8', errors='replace')

positions = {}
if os.path.exists(POS_CACHE) and os.path.getsize(POS_CACHE) > 100:
    with open(POS_CACHE, 'r', encoding='utf-8') as f:
        pos_csv = f.read()
    print(f"  Positionen aus Cache: {POS_CACHE}")
else:
    pos_csv = ""
    # Versuche SIMBAD TAP batch-Abfrage
    try:
        print(f"  Hole Positionen fuer {len(gal_names)} Galaxien via SIMBAD TAP...")
        pos_csv = fetch_simbad_positions(gal_names)
        with open(POS_CACHE, 'w', encoding='utf-8') as f:
            f.write(pos_csv)
        print(f"  Positionen gespeichert: {POS_CACHE}")
    except Exception as e:
        print(f"  SIMBAD TAP fehlgeschlagen: {e}")
        # Fallback: Einzelabfragen via SIMBAD sim-id fuer gängige Namen
        pos_lines = ["id,ra,dec"]
        for name in gal_names:
            try:
                nm_enc = urllib.parse.quote(name)
                url_id = (
                    f"https://simbad.u-strasbg.fr/simbad/sim-id"
                    f"?output.format=ascii&output.params=id,ra(d),dec(d)"
                    f"&submit=submit+id&Ident={nm_enc}"
                )
                req2 = urllib.request.Request(url_id, headers={'User-Agent':'Mozilla/5.0'})
                with urllib.request.urlopen(req2, timeout=10) as r2:
                    txt = r2.read().decode('utf-8', errors='replace')
                # Suche "ID" und "RA" Felder
                for ln in txt.splitlines():
                    if ln.startswith('ra') or ln.startswith('RA'):
                        pass
                    if '|' in ln and not ln.startswith('#'):
                        parts = ln.split('|')
                        if len(parts) >= 3:
                            try:
                                ra_v  = float(parts[0].strip())
                                dec_v = float(parts[1].strip())
                                pos_lines.append(f"{name},{ra_v},{dec_v}")
                                break
                            except ValueError:
                                pass
            except Exception:
                pass
        pos_csv = "\n".join(pos_lines)
        with open(POS_CACHE, 'w', encoding='utf-8') as f:
            f.write(pos_csv)

# CSV parsen
if pos_csv:
    csv_lines = pos_csv.strip().splitlines()
    # Header finden
    start = 0
    for i, line in enumerate(csv_lines):
        if 'ra' in line.lower() or 'RA' in line:
            start = i + 1
            break
    for line in csv_lines[start:]:
        parts = line.split(',')
        if len(parts) < 3:
            continue
        try:
            name = parts[0].strip().strip('"')
            ra   = float(parts[1].strip())
            dec  = float(parts[2].strip())
            if name and 0 <= ra <= 360 and -90 <= dec <= 90:
                positions[name] = (ra, dec)
        except (ValueError, IndexError):
            continue
    print(f"  {len(positions)} Galaxien mit RA/Dec (SIMBAD-Namen)")
# Normalisierung wird nach rar_sigma-Berechnung angewendet (siehe unten)

if len(positions) < 10:
    print("  WARNUNG: Wenige Positionen — KS-Test moeglicherweise eingeschraenkt.")

# ── Koordinaten-Konversion: Aequatorial → Galaktisch ─────────────────────────
def radec_to_galactic(ra_deg, dec_deg):
    """Konvertiert J2000 RA/Dec zu galaktischen Koordinaten l,b."""
    ra  = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    # NGP: RA=192.85948, Dec=+27.12825 (J2000)
    # Galaktisches Zentrum: l=122.93192 Grad von RA-NGP
    RA_NGP  = math.radians(192.85948)
    DEC_NGP = math.radians(27.12825)
    b = math.asin(
        math.sin(dec) * math.sin(DEC_NGP)
        + math.cos(dec) * math.cos(DEC_NGP) * math.cos(ra - RA_NGP)
    )
    l_num = math.cos(dec) * math.sin(ra - RA_NGP)
    l_den = (math.cos(dec) * math.sin(DEC_NGP) * math.cos(ra - RA_NGP)
             - math.sin(dec) * math.cos(DEC_NGP))
    l = math.atan2(l_num, l_den)
    l = math.degrees(l) + 122.93192
    l = l % 360.0
    b = math.degrees(b)
    return l, b

def angle_from_dpole(l_deg, b_deg):
    """Winkelabstand (Grad) vom HTM-Dipol-Pol (l=305, b=25)."""
    l1 = math.radians(l_deg);   b1 = math.radians(b_deg)
    l2 = math.radians(DPOLE_L); b2 = math.radians(DPOLE_B)
    cos_theta = (math.sin(b1) * math.sin(b2)
                 + math.cos(b1) * math.cos(b2) * math.cos(l1 - l2))
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))

def delta_shell(theta_deg):
    """Minimaler Abstand zu einem HTM-Schalenwinkel."""
    return min(abs(theta_deg - s) for s in SHELLS_DEG)

# ── RAR-Streuung pro Galaxie ───────────────────────────────────────────────────
UPSILON_DISK   = 0.5   # M_sun/L_sun bei 3.6 mu
UPSILON_BULGE  = 0.7

rar_sigma  = {}  # sigma_RAR pro Galaxie
n_points   = {}  # Anzahl Datenpunkte

for name, rows in galaxies.items():
    a_obs_list = []
    a_bar_list = []
    for (r_kpc, Vobs, eVobs, Vgas, Vdisk, Vbulge) in rows:
        if r_kpc <= 0 or Vobs <= 0:
            continue
        Vobs2  = Vobs**2
        Vgas2  = np.sign(Vgas) * Vgas**2 if Vgas != 0 else 0.0
        Vdisk2 = np.sign(Vdisk) * Vdisk**2 if Vdisk != 0 else 0.0
        Vbulge2 = np.sign(Vbulge) * Vbulge**2 if Vbulge != 0 else 0.0
        Vbar2  = Vgas2 + UPSILON_DISK * Vdisk2 + UPSILON_BULGE * Vbulge2
        a_o = a_from_V_r(Vobs, r_kpc)
        if Vbar2 <= 0:
            continue
        a_b = a_from_V_r(math.sqrt(abs(Vbar2)), r_kpc)
        if a_o <= 0 or a_b <= 0:
            continue
        a_obs_list.append(a_o)
        a_bar_list.append(a_b)
    if len(a_obs_list) < 3:
        continue
    a_o_arr = np.array(a_obs_list)
    a_b_arr = np.array(a_bar_list)
    a_mdar_arr = a_mdar(a_b_arr)
    ratio = a_o_arr / a_mdar_arr
    valid = ratio > 0
    if np.sum(valid) < 3:
        continue
    log_ratio = np.log10(ratio[valid])
    rar_sigma[name] = float(np.std(log_ratio))
    n_points[name]  = int(np.sum(valid))

print(f"\nRAR-Streuung berechnet fuer {len(rar_sigma)} Galaxien.")

# ── Name-Normalisierung: Matche SIMBAD-Namen an SPARC-Namen ───────────────────
import re as _re

def norm_name(s):
    """Normalisiert Galaxy-Namen: entfernt Leerzeichen, Gross, fuehrende Nullen."""
    s = _re.sub(r'\s+', '', s.strip().upper())
    s = _re.sub(r'(?<=[A-Z])0+([1-9])', r'\1', s)  # DDO064->DDO64, NGC0300->NGC300
    return s

rar_sigma_norm = {norm_name(k): k for k in rar_sigma}
positions_matched = {}
for simbad_name, (ra, dec) in positions.items():
    nn = norm_name(simbad_name)
    if simbad_name in rar_sigma:
        positions_matched[simbad_name] = (ra, dec)
    elif nn in rar_sigma_norm:
        sparc_name = rar_sigma_norm[nn]
        positions_matched[sparc_name] = (ra, dec)
positions = positions_matched
print(f"  Namens-Matching: {len(positions)} SPARC-Galaxien mit Position")

# ── Winkelabstände berechnen ───────────────────────────────────────────────────
theta_D   = {}
delta_sh  = {}
for name, (ra, dec) in positions.items():
    if name not in rar_sigma:
        continue
    l, b = radec_to_galactic(ra, dec)
    tD   = angle_from_dpole(l, b)
    theta_D[name]  = tD
    delta_sh[name] = delta_shell(tD)

matched = {n for n in rar_sigma if n in theta_D}
print(f"  {len(matched)} Galaxien mit RAR-Streuung UND Position.")

# ── KS-Test ────────────────────────────────────────────────────────────────────
near_sigma = [rar_sigma[n] for n in matched if delta_sh[n] < DELTA_SHELL]
far_sigma  = [rar_sigma[n] for n in matched if delta_sh[n] >= DELTA_SHELL]

print(f"\n  Nahe Schalen  (Delta < {DELTA_SHELL} Grad): {len(near_sigma)} Galaxien")
print(f"  Ferne Schalen (Delta >= {DELTA_SHELL} Grad): {len(far_sigma)} Galaxien")

if len(near_sigma) >= 5 and len(far_sigma) >= 5:
    ks_stat, ks_p = stats.ks_2samp(near_sigma, far_sigma)
    mw_stat, mw_p = stats.mannwhitneyu(near_sigma, far_sigma, alternative='two-sided')

    near_med = np.median(near_sigma)
    far_med  = np.median(far_sigma)
    near_mean = np.mean(near_sigma)
    far_mean  = np.mean(far_sigma)

    print(f"\n  Nahe Median sigma_RAR = {near_med:.4f}")
    print(f"  Fern  Median sigma_RAR = {far_med:.4f}")
    print(f"  KS-Statistik = {ks_stat:.4f},  p-Wert = {ks_p:.4f}")
    print(f"  Mann-Whitney U = {mw_stat:.0f},  p-Wert = {mw_p:.4f}")

    ks_significant = ks_p < 0.05
    mw_significant = mw_p < 0.05
else:
    ks_stat = ks_p = mw_stat = mw_p = float('nan')
    near_med = far_med = near_mean = far_mean = float('nan')
    ks_significant = mw_significant = False
    print("  Zu wenige Galaxien fuer statistischen Test.")

# ── Obere Schranke auf delta_L/L aus Null-Ergebnis ───────────────────────────
# Im MOND-Limit (a_bar << g_dag): a_obs ~ sqrt(a0 * a_bar)
# => delta_log10(a_obs) ~ (1/2) * delta_a0/a0 / ln(10) = (1/2) * delta_L/L / ln(10)
# Beobachtete Mittelwertdifferenz: |mean_near - mean_far|
# 95%-Schranke: true_diff < obs_diff + 1.96 * se  (Null-Ergebnis)
# se = sigma_pooled * sqrt(1/n_near + 1/n_far)
MOND_SENS = 0.5 / math.log(10)   # d(log10 a_obs) per unit delta_L/L

upper_lim_dL_over_L = float('nan')
detection_threshold_dL = float('nan')
if near_sigma and far_sigma:
    n_near_n = len(near_sigma)
    n_far_n  = len(far_sigma)
    sigma_pool = float(np.std(near_sigma + far_sigma))
    obs_mean_diff = abs(np.mean(near_sigma) - np.mean(far_sigma))
    se = sigma_pool * math.sqrt(1.0/n_near_n + 1.0/n_far_n)
    upper_lim_sigma = obs_mean_diff + 1.96 * se   # 95% upper limit on true mean diff
    upper_lim_dL_over_L = upper_lim_sigma / MOND_SENS
    detection_threshold_dL = 1.28 * se / MOND_SENS  # 80% power threshold
    print(f"\n  Torsionsschock-Obere Schranke:")
    print(f"    Beob. Mittelwert-Delta sigma_RAR = {obs_mean_diff:.4f} dex")
    print(f"    95% Obergrenze auf echte Differenz = {upper_lim_sigma:.4f} dex")
    print(f"    => |delta_L/L| < {upper_lim_dL_over_L:.2f}  (95% CI)")
    print(f"    Nachweis-Schwelle (80% Power): |delta_L/L| > {detection_threshold_dL:.2f}")

# ── Schalen-Verteilung ─────────────────────────────────────────────────────────
print("\n  Verteilung nach Schalenabstand:")
bins_ranges = [(0,10),(10,20),(20,30),(30,None)]
for lo, hi in bins_ranges:
    count = sum(1 for n in matched
                if delta_sh[n] >= lo and (hi is None or delta_sh[n] < hi))
    med_bin = [rar_sigma[n] for n in matched
               if delta_sh[n] >= lo and (hi is None or delta_sh[n] < hi)]
    if med_bin:
        label = f"  [{lo:2d},{hi if hi else '90'}] Grad:"
        print(f"{label:15s} {count:3d} Galaxien,  Median sigma_RAR = {np.median(med_bin):.4f}")

# ── Plot ────────────────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax1, ax2 = axes

    all_sigma = [rar_sigma[n] for n in matched]
    all_delta = [delta_sh[n]  for n in matched]

    ax1.scatter(all_delta, all_sigma, s=15, alpha=0.6, color='steelblue')
    for s in SHELLS_DEG:
        if s < 90:
            ax1.axvline(s, color='orange', linestyle='--', alpha=0.5, lw=1)
    ax1.axvline(DELTA_SHELL, color='red', linestyle=':', lw=2, label=f'Schwelle {DELTA_SHELL} Grad')
    ax1.set_xlabel('Delta-Theta zu naechster HTM-Schale [Grad]')
    ax1.set_ylabel('sigma_RAR (log10-Streuung)')
    ax1.set_title('RAR-Streuung vs Shell-Abstand (SPARC)')
    ax1.legend(fontsize=8)

    bins = np.linspace(0, max([max(near_sigma or [0]), max(far_sigma or [0.1])]) * 1.1, 20)
    if near_sigma:
        ax2.hist(near_sigma, bins=bins, alpha=0.6, color='orange', label=f'Nahe (<{DELTA_SHELL}°, n={len(near_sigma)})')
    if far_sigma:
        ax2.hist(far_sigma,  bins=bins, alpha=0.6, color='steelblue', label=f'Fern  (>={DELTA_SHELL}°, n={len(far_sigma)})')
    ax2.set_xlabel('sigma_RAR')
    ax2.set_ylabel('Anzahl Galaxien')
    ax2.set_title(f'KS p={ks_p:.3f}, MW p={mw_p:.3f}')
    ax2.legend(fontsize=8)

    plt.suptitle('OT-18: RAR-Streuung vs HTM-Schalenabstand', y=1.02)
    plt.tight_layout()
    plot_path = os.path.join(RESULTS, "OT_18_plot.png")
    plt.savefig(plot_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot: {plot_path}")
except Exception as e:
    print(f"\n  Plot fehlgeschlagen: {e}")

# ── Ergebnistext ──────────────────────────────────────────────────────────────
def fmt_or(val, fmt=".4f"):
    return format(val, fmt) if not math.isnan(val) else "N/A"

lines = [
    "=" * 70,
    "OT-18: KS-Test RAR-Residuen vs HTM-Schalenabstand (SPARC)",
    "=" * 70,
    "",
    f"  Analyse-Basis:   {len(rar_sigma)} SPARC-Galaxien mit RAR-Daten",
    f"  Mit Position:    {len(matched)} Galaxien (RA/Dec aus VizieR)",
    f"  Schalen:         n*58.65 Grad, n=1..6",
    f"  HTM-Dipol-Pol:   l={DPOLE_L}, b={DPOLE_B}",
    f"  Naehe-Schwelle:  Delta-theta < {DELTA_SHELL} Grad",
    "",
    "-" * 70,
    "Ergebnisse",
    "-" * 70,
    "",
    f"  Nahe-Schalen-Galaxien:  {len(near_sigma)}",
    f"  Ferne-Schalen-Galaxien: {len(far_sigma)}",
    "",
    f"  Median sigma_RAR (nahe):  {fmt_or(near_med)}",
    f"  Median sigma_RAR (fern):  {fmt_or(far_med)}",
    f"  Mittel sigma_RAR (nahe):  {fmt_or(near_mean)}",
    f"  Mittel sigma_RAR (fern):  {fmt_or(far_mean)}",
    "",
    f"  KS-Statistik = {fmt_or(ks_stat)},  p = {fmt_or(ks_p)}",
    f"  Mann-Whitney U = {fmt_or(mw_stat, '.0f')},  p = {fmt_or(mw_p)}",
    "",
    "-" * 70,
    "BEWERTUNG",
    "-" * 70,
    "",
]
if math.isnan(ks_p):
    lines += [
        "  Test nicht ausfuehrbar: zu wenige Galaxien oder fehlende Positionen.",
        "  Ein valides Ergebnis erfordert >= 5 Galaxien pro Gruppe.",
]
elif ks_significant or mw_significant:
    diff_sign = "groessere" if near_med > far_med else "kleinere"
    lines += [
        f"  OT-18 POSITIV: Galaxien nahe HTM-Schalen zeigen {diff_sign}",
        f"  RAR-Streuung als ferne Galaxien (p_KS={fmt_or(ks_p)}, p_MW={fmt_or(mw_p)}).",
        "",
        "  VORSICHT: Statistisch signifikant (p<0.05) bedeutet nicht kausal.",
        "  Moegliche Confounder: Selektion, Entfernung, Wellenlänge.",
]
else:
    dL_str  = (f"{upper_lim_dL_over_L:.2f}" if not math.isnan(upper_lim_dL_over_L) else "N/A")
    det_str = (f"{detection_threshold_dL:.2f}" if not math.isnan(detection_threshold_dL) else "N/A")
    lines += [
        f"  OT-18: KEIN TORSIONSSCHOCK-SIGNAL DETEKTIERT",
        f"  p_KS = {fmt_or(ks_p)}, p_MW = {fmt_or(mw_p)} -- kein signifikanter Unterschied.",
        "",
        "  TORSIONSSCHOCK-INTERPRETATION:",
        "  L ist global im Gleichgewicht -- aber Torsionsschocks bei Schalendurchgaengen",
        "  koennen L lokal kurzzeitig erhoehen (S3-Hopf-Index-Sprung an theta_n).",
        "  Galaxien nahe Schalenraendern koennen ein 'fossiles' delta_a0-Signal zeigen.",
        "",
        "  Dieses Null-Ergebnis liefert eine OBERE SCHRANKE:",
        f"    |delta_L/L| < {dL_str}  (95% Konfidenz)",
        f"    Test sensibel ab: |delta_L/L| > {det_str}  (80% Nachweis-Power)",
        "",
        "  STATUS: KEIN SIGNAL -- kein 'predicted null', sondern echter Test mit",
        "  Oberer Schranke auf Torsionsschock-Staerke.",
        "",
        "  Falsifikation des Null-Hypothese (Schock-Signal vorhanden) waere:",
        "  signifikantes sigma_RAR(nahe) > sigma_RAR(fern) mit p < 0.05.",
]
lines += ["", "=" * 70]

result = "\n".join(lines)
print("\n" + result)
out = os.path.join(RESULTS, "OT_18_result.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write(result)
print(f"\n  Ergebnis: {out}")
