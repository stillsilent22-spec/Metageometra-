"""
OT-23: Kosmische Voids vs HTM-Schalenwinkel
============================================
Aufgabe: Teste ob SDSS-Void-ZENTREN bevorzugt in der Naehe von HTM-Schalen-
winkeln (theta_n = n * 58.65 Grad vom HTM-Dipol-Pol) liegen.

Datenbasis:
  Prioritaet 1: cosmicvoids.net SDSS DR7/DR10 Void-Katalog
  Prioritaet 2: Sutter et al. 2012, VizieR J/MNRAS/431/2307
  Prioritaet 3: Nordhaus & Brandt 2020 oder vergleichbar
  Fallback: synthetischer Test mit Uniform-Null-Hypothese

Methodik:
  Fuer jeden Void: RA/Dec → Galaktisch l,b → theta_D vom HTM-Pol → Delta_shell
  Null-Hypothese: Voids sind zufaellig verteilt → Delta_shell ~ uniform
  Test: KS-Test p(Delta_shell) gegen Uniform(0, 29.325)
  Zusaetzlich: Chi-Squared Test in Schalennaehe-Bins
"""

import os, sys, math, io, csv
import numpy as np
from scipy import stats
import urllib.request
import urllib.parse

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

# ── Framework-Konstanten ──────────────────────────────────────────────────────
SHELLS_DEG = [58.65 * n for n in range(1, 7)]
DPOLE_L    = 305.0
DPOLE_B    = 25.0
HALF_SHELL = 58.65 / 2.0   # halb Schneckabstand = 29.325 Grad

# ── Koordinaten-Hilfen ────────────────────────────────────────────────────────
def radec_to_galactic(ra_deg, dec_deg):
    ra  = math.radians(ra_deg)
    dec = math.radians(dec_deg)
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
    return l % 360.0, math.degrees(b)

def angle_from_dpole(l_deg, b_deg):
    l1 = math.radians(l_deg);   b1 = math.radians(b_deg)
    l2 = math.radians(DPOLE_L); b2 = math.radians(DPOLE_B)
    cos_t = (math.sin(b1)*math.sin(b2) + math.cos(b1)*math.cos(b2)*math.cos(l1-l2))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_t))))

def delta_shell(theta_deg):
    return min(abs(theta_deg - s) for s in SHELLS_DEG)

def fetch_url(url, timeout=30):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', errors='replace')

# ── Void-Kataloge versuchen ──────────────────────────────────────────────────
voids = []   # list of (ra_deg, dec_deg, radius_Mpc, z) or (ra, dec)
source_label = "kein Katalog"

# ── Versuch 1: Sutter et al. 2012 VizieR CDS (J/MNRAS/431/2307) ──────────────
SUTTER_CACHE = os.path.join(RESULTS, "sutter_voids.dat")
SUTTER_CACHE_CSV = os.path.join(RESULTS, "sutter_voids.csv")
SUTTER_URL   = "https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/431/2307/"

print("Versuche Sutter+2012 Void-Katalog (VizieR)...")
if os.path.exists(SUTTER_CACHE_CSV) and os.path.getsize(SUTTER_CACHE_CSV) > 100:
    with open(SUTTER_CACHE_CSV,'r',encoding='ascii',errors='replace') as f:
        sutter_data = f.read()
    print(f"  Aus CSV-Cache: {SUTTER_CACHE_CSV}")
elif os.path.exists(SUTTER_CACHE) and os.path.getsize(SUTTER_CACHE) > 100:
    with open(SUTTER_CACHE,'r',encoding='ascii',errors='replace') as f:
        sutter_data = f.read()
    print(f"  Aus Cache: {SUTTER_CACHE}")
else:
    try:
        # Versuche TAP query
        q = ('SELECT VoidID, RAJ2000, DEJ2000, z, Rvoid '
             'FROM "J/MNRAS/431/2307/table3"')
        tap_url = (
            "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
            "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv"
            f"&QUERY={urllib.parse.quote(q)}"
        )
        sutter_data = fetch_url(tap_url, timeout=45)
        with open(SUTTER_CACHE,'w',encoding='ascii',errors='replace') as f:
            f.write(sutter_data)
        print(f"  Sutter+2012 geladen via TAP")
    except Exception as e:
        print(f"  Sutter+2012 TAP: {e}")
        sutter_data = ""

# Sutter-Daten parsen (CSV mit VoidID, RAJ2000, DEJ2000, z, Rvoid)
if sutter_data and 'RAJ2000' in sutter_data:
    reader = csv.DictReader(io.StringIO(sutter_data))
    for row in reader:
        try:
            ra  = float(row.get('RAJ2000') or row.get('_RAJ2000',''))
            dec = float(row.get('DEJ2000') or row.get('_DEJ2000',''))
            z   = float(row.get('z','0') or 0)
            r   = float(row.get('Rvoid','0') or 0)
            voids.append((ra, dec, r, z))
        except (ValueError, TypeError):
            continue
    if voids:
        source_label = f"Sutter+2012 (J/MNRAS/431/2307) [{len(voids)} Voids]"
        print(f"  {len(voids)} Voids aus Sutter+2012")

# ── Versuch 2: Alternativer Sutter TAP table1 oder table2 ──────────────────
if not voids:
    for tbl in ["table1", "table2", "voids"]:
        try:
            q2 = f'SELECT * FROM "J/MNRAS/431/2307/{tbl}" MAXREC=2000'
            tap2 = (
                "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
                "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv"
                f"&QUERY={urllib.parse.quote(q2)}"
            )
            data2 = fetch_url(tap2, timeout=45)
            lines2 = data2.strip().splitlines()
            if len(lines2) > 2 and 'error' not in data2.lower()[:200]:
                print(f"  Tabelle J/MNRAS/431/2307/{tbl}: {len(lines2)} Zeilen")
                print(f"  Header: {lines2[0][:100]}")
                # Suche RA/Dec Spalten
                header2 = lines2[0].split(',')
                ra_col  = next((i for i,h in enumerate(header2) if 'RA' in h.upper()), None)
                dec_col = next((i for i,h in enumerate(header2) if 'DE' in h.upper() or 'DEC' in h.upper()), None)
                if ra_col is not None and dec_col is not None:
                    for line in lines2[1:]:
                        parts = line.split(',')
                        try:
                            ra  = float(parts[ra_col])
                            dec = float(parts[dec_col])
                            voids.append((ra, dec, 0, 0))
                        except (ValueError, IndexError):
                            continue
                    if voids:
                        source_label = f"Sutter+2012 {tbl} [{len(voids)} Voids]"
                        print(f"  {len(voids)} Voids aus {tbl}")
                        break
        except Exception as e2:
            print(f"  {tbl}: {e2}")

# ── Versuch 3: Hovell & Vogeley 2004 (J/AJ/128/1) ──────────────────────────
if not voids:
    print("Versuche Hoyle & Vogeley 2004 (J/AJ/128/1)...")
    try:
        q3 = 'SELECT * FROM "J/AJ/128/1/table1" MAXREC=2000'
        tap3 = (
            "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
            "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv"
            f"&QUERY={urllib.parse.quote(q3)}"
        )
        data3 = fetch_url(tap3, timeout=45)
        lines3 = data3.strip().splitlines()
        if len(lines3) > 2:
            print(f"  J/AJ/128/1: {len(lines3)} Zeilen, Header: {lines3[0][:100]}")
            header3 = lines3[0].split(',')
            ra_col  = next((i for i,h in enumerate(header3) if 'RA' in h.upper()), None)
            dec_col = next((i for i,h in enumerate(header3) if 'DE' in h.upper()), None)
            if ra_col is not None and dec_col is not None:
                for line in lines3[1:]:
                    parts = line.split(',')
                    try:
                        voids.append((float(parts[ra_col]), float(parts[dec_col]), 0, 0))
                    except Exception:
                        continue
                source_label = f"Hoyle & Vogeley 2004 [{len(voids)} Voids]"
                print(f"  {len(voids)} Voids geladen")
    except Exception as e3:
        print(f"  Hoyle & Vogeley: {e3}")

# ── Versuch 4: cosmicvoids.net public files ──────────────────────────────────
if not voids:
    print("Versuche cosmicvoids.net...")
    for cv_url in [
        "https://www.cosmicvoids.net/voids/SDSS/DR7/dim3/central_RA_Dec.dat",
        "https://www.cosmicvoids.net/voids/SDSS/DR7/dim3/voids.dat",
        "https://www.cosmicvoids.net/data/voids/SDSS_DR7_voids.csv",
    ]:
        try:
            data_cv = fetch_url(cv_url, timeout=15)
            lines_cv = data_cv.strip().splitlines()
            print(f"  {cv_url}: {len(lines_cv)} Zeilen")
            for line in lines_cv[:5]:
                print(f"    {line[:80]}")
            for line in lines_cv:
                parts = line.strip().split()
                try:
                    ra, dec = float(parts[0]), float(parts[1])
                    if 0 <= ra <= 360 and -90 <= dec <= 90:
                        voids.append((ra, dec, 0, 0))
                except Exception:
                    continue
            if voids:
                source_label = f"cosmicvoids.net [{len(voids)} Voids]"
                break
        except Exception as e_cv:
            print(f"  {cv_url}: {e_cv}")

# ── Ergebnis ohne Daten ────────────────────────────────────────────────────────
if not voids:
    print("\n  KEIN Void-Katalog verfügbar. Erzeuge Null-Simulation zum Test.")
    # Simuliere: Was würde man bei uniformer Verteilung erwarten?
    np.random.seed(42)
    n_sim = 1000
    ra_sim  = np.random.uniform(0, 360, n_sim)
    dec_sim = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_sim)))
    for i in range(n_sim):
        voids.append((ra_sim[i], dec_sim[i], 0, 0))
    source_label = f"SIMULATION: {n_sim} uniform verteilte Punkte (kein echter Katalog)"

print(f"\n  Verwende: {source_label}")
print(f"  Anzahl Voids: {len(voids)}")

# ── Delta-Shell fuer alle Voids ───────────────────────────────────────────────
delta_shells = []
thetas_D     = []
for (ra, dec, r_v, z) in voids:
    l, b = radec_to_galactic(ra, dec)
    tD = angle_from_dpole(l, b)
    ds = delta_shell(tD)
    thetas_D.append(tD)
    delta_shells.append(ds)

delta_arr = np.array(delta_shells)
theta_arr = np.array(thetas_D)

# ── KS-Test gegen Null-Verteilung (numerisch berechnet) ──────────────────────
# Die erwartete Null-Verteilung von delta_shell fuer gleichmaessig verteilte
# Punkte auf der Kugel ist NICHT uniform, weil P(theta) = sin(theta)/2.
# Wir berechnen die erwartete CDF numerisch.
np.random.seed(0)
n_null = 200000
# Uniform auf Kugel: theta ~ arccos(U(-1,1)), phi ~ U(0,360)
null_theta_rad = np.arccos(np.random.uniform(-1, 1, n_null))
null_theta_deg = np.degrees(null_theta_rad)
null_delta = np.array([delta_shell(th) for th in null_theta_deg])

# KS-Test: real voids vs simulierte Null-Verteilung
ks_stat, ks_p = stats.ks_2samp(delta_arr, null_delta)

# Bins: 0-5, 5-10, 10-15, 15-20, 20-25, 25-29.3 Grad
bin_edges = [0, 5, 10, 15, 20, 25, HALF_SHELL]
observed_counts, _ = np.histogram(delta_arr, bins=bin_edges)
expected_counts, _ = np.histogram(null_delta, bins=bin_edges)
# Skaliere expected auf selbe Gesamtzahl (Anteil in-range)
n_obs_in = observed_counts.sum()
n_null_in = expected_counts.sum()
if n_null_in > 0:
    expected_counts_scaled = expected_counts * (n_obs_in / n_null_in)
else:
    expected_counts_scaled = expected_counts * 1.0

# Chi2-Test nur wenn genug Galaxien in jedem Bin
if np.all(expected_counts_scaled > 5) and n_obs_in >= 20:
    chi2_stat, chi2_p = stats.chisquare(observed_counts,
                                         f_exp=expected_counts_scaled)
else:
    chi2_stat, chi2_p = float('nan'), float('nan')

# Null-Verteilung: Erwartete in-range Fraktion
frac_in_null = n_null_in / n_null
frac_in_obs  = n_obs_in / len(delta_arr) if len(delta_arr) > 0 else 0.0

print(f"\n  KS-Test (2-Stichproben):   D = {ks_stat:.4f},  p = {ks_p:.4f}")
chi2_str = f"{chi2_stat:.2f}" if not math.isnan(chi2_stat) else "N/A"
chi2p_str = f"{chi2_p:.4f}" if not math.isnan(chi2_p) else "N/A"
print(f"  Chi2-Test (6 Bins):        chi2 = {chi2_str},  p = {chi2p_str}")
print(f"\n  Bin-Analyse:")
for i, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
    esci = expected_counts_scaled[i]
    print(f"    [{lo:.1f}, {hi:.1f}] Grad: {observed_counts[i]:5.0f} beob, {esci:5.1f} erw")

# ── Schalenwinkel-Bevorzugung ─────────────────────────────────────────────────
pct_near = 100 * np.mean(delta_arr < 10.0)
print(f"\n  Anteil Voids nahe Schale (<10 Grad): {pct_near:.1f}%")
print(f"  Erwartet bei uniform: {100*10.0/HALF_SHELL:.1f}%")

# ── Plot ────────────────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax1, ax2 = axes

    ax1.hist(delta_arr, bins=30, density=True, alpha=0.7, color='steelblue',
             label=f'Voids (n={len(voids)})')
    x_uni = np.linspace(0, HALF_SHELL, 100)
    ax1.plot(x_uni, np.ones_like(x_uni)/HALF_SHELL, 'r--', lw=2, label='Uniform H0')
    ax1.set_xlabel('Abstand zur naechsten HTM-Schale [Grad]')
    ax1.set_ylabel('Dichte')
    ax1.set_title(f'OT-23: Void-Verteilung vs HTM-Schalen\n{source_label[:50]}')
    ax1.legend(fontsize=8)

    ax2.hist(theta_arr, bins=36, color='steelblue', alpha=0.7)
    for s in SHELLS_DEG:
        ax2.axvline(s, color='orange', linestyle='--', alpha=0.7, lw=1.5,
                   label='HTM Schale' if s == SHELLS_DEG[0] else None)
    ax2.set_xlabel('theta_D vom HTM-Pol [Grad]')
    ax2.set_ylabel('Anzahl Voids')
    ax2.set_title(f'KS p={ks_p:.3f},  Chi2 p={chi2_p:.3f}')
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plot_path = os.path.join(RESULTS, "OT_23_plot.png")
    plt.savefig(plot_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot: {plot_path}")
except Exception as e:
    print(f"\n  Plot fehlgeschlagen: {e}")

# ── Ergebnistext ──────────────────────────────────────────────────────────────
is_sim = "SIMULATION" in source_label
sig_ks   = ks_p < 0.05
sig_chi2 = chi2_p < 0.05

lines = [
    "=" * 70,
    "OT-23: Kosmische Voids vs HTM-Schalenwinkel",
    "=" * 70,
    "",
    f"  Katalog:        {source_label}",
    f"  Voids gesamt:   {len(voids)}",
    f"  Schalen:        n*58.65 Grad, n=1..6",
    f"  HTM-Dipol-Pol:  l={DPOLE_L}, b={DPOLE_B}",
    "",
    "-" * 70,
    "Statistik",
    "-" * 70,
    "",
    f"  Anteil Voids nahe Schale (<10 Grad): {pct_near:.1f}%",
    f"  Erwartet (uniform):                  {100*10.0/HALF_SHELL:.1f}%",
    "",
    f"  KS-Test (2-Stichproben vs Null):      D = {ks_stat:.4f},  p = {ks_p:.4f}",
    f"  Chi2-Test (6 Bins):                   chi2 = {chi2_str},  p = {chi2p_str}",
    "",
    "  Bin breakdown (Delta-Schale):",
]
for i, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
    esc = expected_counts_scaled[i]
    pct_diff = 100*(observed_counts[i]-esc)/max(esc, 0.1)
    lines.append(
        f"    [{lo:.1f}-{hi:.1f}] Grad: {observed_counts[i]:5.0f} beob, "
        f"{esc:5.1f} erw  ({pct_diff:+.1f}%)"
    )
lines += [
    "",
    "-" * 70,
    "BEWERTUNG",
    "-" * 70,
    "",
]

if is_sim:
    lines += [
        "  OT-23 NICHT AUSFUEHRBAR: Kein echter Void-Katalog verfügbar.",
        "  cosmicvoids.net, VizieR J/MNRAS/431/2307 und J/AJ/128/1 waren",
        "  nicht erreichbar oder enthielten keine verwertbaren Positionen.",
        "  Die Simulation zeigt erwartungsgemäss: uniform verteilte Punkte",
        f"  ergeben p={ks_p:.3f} (kein Signal). Echter Test benoetigt echten Katalog.",
        "",
        "  EMPFEHLUNG: Sutter et al. 2012 SDSS void catalog manuell herunterladen",
        "  von cosmicvoids.net (nach Login) oder arXiv:1210.6446.",
]
elif sig_ks or sig_chi2:
    lines += [
        f"  OT-23 POSITIV: Voids haefen sich signifikant bei HTM-Schalenabstaenden",
        f"  (KS p={ks_p:.4f}, Chi2 p={chi2_p:.4f}).",
        "",
        "  VORSICHT: Selektionseffekte? SDSS-Footprint nicht uniform.",
        "  Der HTM-Pol liegt bei l=305, b=25 — Galaktische Ebene schneidet SDSS.",
]
else:
    lines += [
        f"  OT-23 NEGATIV: Keine signifikante Haeufung von Void-Zentren",
        f"  bei HTM-Schalenabstaenden (KS p={ks_p:.4f}, Chi2 p={chi2_p:.4f}).",
        "",
        "  Die Verteilung der kosmischen Voids ist vertraeglich mit",
        "  gleichmaessiger Winkelverteilung bezueglich des HTM-Pols.",
]
lines += ["", "=" * 70]

result = "\n".join(lines)
print("\n" + result)
out = os.path.join(RESULTS, "OT_23_result.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write(result)
print(f"\n  Ergebnis: {out}")
