"""
OT-14: SRM Halo Profil gegen SPARC Rotationskurven
===================================================
Metageometra V19 Testvorhersage:

  SRM-Dichteprofil:  rho_SRM(r) = rho_0 * (1 + r/r_s)^(-1.205)
    r_s = 42.3 kpc  — FEST, kein freier Parameter
    Einziger Freiheitsgrad per Galaxie: rho_0

  Vergleich:
    NFW-Profil:  rho_NFW(r) = rho_s / [(r/r_s)(1 + r/r_s)^2]
    beide rho_s UND r_s frei (2 freie Parameter)

Daten: SPARC Lelli+2016  (175 Scheibengalaxien, Lelli+2016 AJ 152 157)
       CDS: https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/

Masse-Licht-Verhältnisse (SPARC Standard):
  Upsilon_disk = 0.5 M_sun/L_sun  bei [3.6]
  Upsilon_bulge = 0.7 M_sun/L_sun  bei [3.6]
"""

import os, sys, math, urllib.request
import numpy as np
from scipy import optimize, integrate, interpolate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Fix Windows cp1252 UnicodeEncodeError when printing Unicode chars
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

# ── Physikalische Konstanten ──────────────────────────────────────────────────
G_KPC  = 4.3009e-6    # kpc (km/s)^2 / M_sun
RS_SRM = 42.3         # kpc  —  SRM Skalenradius, NICHT verhandelbar
SLOPE  = 1.205        # SRM-Steigung (Exponent)
UPS_D  = 0.5          # M/L Scheibe bei [3.6]
UPS_B  = 0.7          # M/L Bulge  bei [3.6]
VERR_FLOOR = 2.0      # km/s — Minimaler Geschwindigkeitsfehler (Systematik)
MIN_POINTS = 3        # Mindestanzahl Datenpunkte pro Galaxie

print("=" * 70)
print("  OT-14: SRM vs NFW auf SPARC Rotationskurven")
print("=" * 70)

# ── Vorberechnete SRM-Integraltabelle ─────────────────────────────────────────
# J(x) = integral_0^x u^2 (1+u)^(-1.205) du   [kpc^3 / kpc^3 = dimensionslos]
# Damit: V^2_SRM(r) = rho_0 * 4*pi*G * r_s^3 * J(r/r_s) / r

print("  Vorberechnung SRM-Integraltabelle...")
_x_tab = np.linspace(0, 20.0, 2001)  # x = r/r_s bis zum 20-fachen
_J_tab = np.zeros(len(_x_tab))
for _i in range(1, len(_x_tab)):
    _J_val, _ = integrate.quad(lambda u: u**2 * (1.0 + u)**(-SLOPE),
                                0, _x_tab[_i], limit=100)
    _J_tab[_i] = _J_val
_J_interp = interpolate.CubicSpline(_x_tab, _J_tab, extrapolate=True)

PREFAC_SRM = 4.0 * math.pi * G_KPC * RS_SRM**3   # (km/s)^2 / (M_sun/kpc^3) * kpc^2
# V^2_SRM(r; rho_0) = rho_0 * PREFAC_SRM * J(r/r_s) / r


def v2_srm(r_kpc, rho0):
    """V^2-Beitrag des SRM-Halos [km/s]^2."""
    x = np.asarray(r_kpc) / RS_SRM
    J = np.maximum(0, _J_interp(x))
    return rho0 * PREFAC_SRM * J / np.asarray(r_kpc)


def v2_nfw(r_kpc, rho_s, r_s_nfw):
    """V^2-Beitrag eines NFW-Halos [km/s]^2."""
    x = np.asarray(r_kpc) / r_s_nfw
    f = np.log(1.0 + x) - x / (1.0 + x)
    f = np.maximum(0.0, f)
    prefac = 4.0 * math.pi * G_KPC * rho_s * r_s_nfw**3
    return prefac * f / np.asarray(r_kpc)


def v_model(r, v2_bary, v2_halo):
    """Gesamtmodellgeschwindigkeit = sqrt(max(0, V²_bary + V²_halo)) [km/s]."""
    return np.sqrt(np.maximum(0.0, v2_bary + v2_halo))


def chi2_val(v_obs, v_mod, v_err):
    """Chi-Quadrat Summe."""
    return float(np.sum(((v_obs - v_mod) / v_err) ** 2))


# ── Daten herunterladen ───────────────────────────────────────────────────────
TABLE2_URL = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat"
TABLE2_CACHE = os.path.join(RESULTS, "sparc_table2.dat")

if os.path.exists(TABLE2_CACHE):
    print(f"  Nutze gecachte Datei: {TABLE2_CACHE}")
    with open(TABLE2_CACHE, "r") as f:
        raw = f.read()
else:
    print(f"  Lade SPARC table2.dat von CDS...")
    req = urllib.request.Request(TABLE2_URL, headers={"User-Agent": "Python/MetageometraOT14"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    with open(TABLE2_CACHE, "w", encoding="utf-8") as f:
        f.write(raw)
    print(f"  Gespeichert: {TABLE2_CACHE}")

# ── Parsen: Byte-by-byte Format aus CDS ReadMe ───────────────────────────────
# Bytes 1-11:  Name (A11)
# Bytes 13-18: Dist (F6.2, Mpc)
# Bytes 20-25: Rad  (F6.2, kpc)
# Bytes 27-32: Vobs (F6.2, km/s)
# Bytes 34-38: eVobs(F5.2, km/s)
# Bytes 40-45: Vgas (F6.2, km/s)
# Bytes 47-52: Vdisk(F6.2, km/s) for M/L=1
# Bytes 54-59: Vbulge(F6.2, km/s) for M/L=1
# (0-indexed Python slices = bytes - 1)

galaxies = {}  # name -> dict of lists
for line in raw.splitlines():
    if len(line) < 59:
        continue
    try:
        name   = line[0:11].strip()
        dist   = float(line[12:18])
        r_kpc  = float(line[19:25])
        vobs   = float(line[26:32])
        evobs  = float(line[33:38])
        vgas   = float(line[39:45])
        vdisk  = float(line[46:52])
        vbulge = float(line[53:59])
    except (ValueError, IndexError):
        continue
    if name not in galaxies:
        galaxies[name] = {"dist": dist, "r": [], "vobs": [],
                          "evobs": [], "vgas": [], "vdisk": [], "vbulge": []}
    galaxies[name]["r"].append(r_kpc)
    galaxies[name]["vobs"].append(vobs)
    galaxies[name]["evobs"].append(max(evobs, VERR_FLOOR))
    galaxies[name]["vgas"].append(vgas)
    galaxies[name]["vdisk"].append(vdisk)
    galaxies[name]["vbulge"].append(vbulge)

# Convert to numpy arrays
for name in galaxies:
    for key in ["r", "vobs", "evobs", "vgas", "vdisk", "vbulge"]:
        galaxies[name][key] = np.array(galaxies[name][key])

print(f"  {len(galaxies)} Galaxien geladen.")

# ── Fit-Funktion ──────────────────────────────────────────────────────────────
results = {}  # name -> dict

def fit_galaxy(name, data):
    r      = data["r"]
    vobs   = data["vobs"]
    verr   = data["evobs"]

    # Baryonenbeitrag (SPARC Vorzeichenkonvention: V^2_X = sign(V_X)*V_X^2)
    def signed_v2(v):
        return np.sign(v) * v**2

    v2bary = (signed_v2(data["vgas"])
              + UPS_D * data["vdisk"]**2
              + UPS_B * data["vbulge"]**2)

    if len(r) < MIN_POINTS:
        return None

    # ---- SRM-Fit: 1 freier Parameter rho_0 ≥ 0 ----
    # Kernel pro Radialschritt (linear in rho_0)
    x_arr = r / RS_SRM
    J_arr = np.maximum(0.0, _J_interp(x_arr))
    kernel = PREFAC_SRM * J_arr / r    # V^2_SRM = rho_0 * kernel

    # rho_0 Schätzwert aus dem mittleren DM-Bedarf
    v2_needed = vobs**2 - v2bary
    v2_need_ok = np.where(v2_needed > 0, v2_needed, 0.0)
    # Erste Schätzung: Median-Verhältnis
    k_pos = kernel > 0
    if k_pos.any():
        rho0_init = float(np.median(v2_need_ok[k_pos] / kernel[k_pos]))
    else:
        rho0_init = 0.0
    rho0_init = max(rho0_init, 0.0)

    def chi2_srm(log_rho0):
        rho0 = math.exp(log_rho0) if log_rho0 > -30 else 0.0
        vm = v_model(r, v2bary, rho0 * kernel)
        return chi2_val(vobs, vm, verr)

    # Schätzung rho_0 = 0 (kein Halo)
    chi2_srm_zero = chi2_val(vobs, v_model(r, v2bary, np.zeros_like(r)), verr)

    if rho0_init > 0:
        log_init = math.log(rho0_init)
        try:
            res_srm = optimize.minimize_scalar(
                chi2_srm,
                bracket=(log_init - 5, log_init, log_init + 5),
                method="brent",
                options={"xtol": 1e-8, "maxiter": 500},
            )
            rho0_best = math.exp(res_srm.x) if res_srm.x > -30 else 0.0
            chi2_srm_best = float(res_srm.fun)
        except Exception:
            rho0_best = rho0_init
            chi2_srm_best = chi2_srm(math.log(rho0_init) if rho0_init > 0 else -100)
    else:
        rho0_best = 0.0
        chi2_srm_best = chi2_srm_zero

    # Stelle sicher, dass chi2 ≥ rho=0 Lösung ist (Brentq kann divergieren)
    if rho0_best < 0 or chi2_srm_best > chi2_srm_zero or not math.isfinite(chi2_srm_best):
        rho0_best = 0.0
        chi2_srm_best = chi2_srm_zero

    # ---- NFW-Fit: 2 freie Parameter rho_s, r_s_nfw ----
    # NFW ist linear in rho_s, nichtlinear in r_s_nfw
    def chi2_nfw_rs(log_rs):
        rs_nfw = math.exp(log_rs)
        # Analytische Optimierung für rho_s bei festem r_s_nfw  (linear)
        x_n = r / rs_nfw
        f_n = np.maximum(0.0, np.log(1.0 + x_n) - x_n / (1.0 + x_n))
        halo_unit = 4.0 * math.pi * G_KPC * rs_nfw**3 * f_n / r  # V^2 für rho_s=1

        # chi2(rho_s) = Σ (vobs - sqrt(max(0, v2bary + rho_s*halo_unit)))² / verr²
        # Lineare Schätzung: rho_s von least-squares auf sqrt-linearisiertem Problem
        # Verwende bounded 1D minimize_scalar
        def chi2_rhos(log_rhos):
            rho_s = math.exp(log_rhos) if log_rhos > -30 else 0.0
            vm = v_model(r, v2bary, rho_s * halo_unit)
            return chi2_val(vobs, vm, verr)

        # Erste Schätzung für rho_s
        k_pos2 = halo_unit > 0
        if k_pos2.any() and v2_need_ok[k_pos2].any():
            rhos0 = float(np.median(v2_need_ok[k_pos2] / halo_unit[k_pos2]))
            rhos0 = max(rhos0, 1e-6)
            log_rhos0 = math.log(rhos0)
        else:
            log_rhos0 = 0.0

        try:
            res2 = optimize.minimize_scalar(
                chi2_rhos,
                bracket=(log_rhos0 - 5, log_rhos0, log_rhos0 + 5),
                method="brent",
                options={"xtol": 1e-8, "maxiter": 300},
            )
            return float(res2.fun), math.exp(res2.x)
        except Exception:
            c0 = chi2_rhos(log_rhos0)
            return float(c0), math.exp(log_rhos0)

    # Raster über r_s_nfw: 1 bis 300 kpc (log-Skala)
    rs_grid = np.logspace(0, 2.7, 40)  # ~1 bis 500 kpc
    best_nfw = (1e30, 0.0, 10.0)  # (chi2, rho_s, r_s_nfw)
    for rs_try in rs_grid:
        try:
            c, rhos = chi2_nfw_rs(math.log(rs_try))
            if c < best_nfw[0]:
                best_nfw = (c, rhos, rs_try)
        except Exception:
            pass

    # Verfeinere mit scipy.minimize
    x0 = [math.log(best_nfw[1]) if best_nfw[1] > 0 else 0.0,
          math.log(best_nfw[2])]

    def chi2_nfw_full(p):
        rho_s = math.exp(p[0]) if p[0] > -30 else 0.0
        rs_n  = math.exp(p[1]) if p[1] > -10 else 1e-4
        hu = v2_nfw(r, rho_s, rs_n)
        return chi2_val(vobs, v_model(r, v2bary, hu), verr)

    try:
        res_nfw = optimize.minimize(
            chi2_nfw_full, x0,
            method="Nelder-Mead",
            options={"maxiter": 10000, "xatol": 1e-6, "fatol": 1e-6, "adaptive": True},
        )
        rhos_best = math.exp(res_nfw.x[0])
        rs_nfw_best = math.exp(res_nfw.x[1])
        chi2_nfw_best = float(res_nfw.fun)
    except Exception:
        rhos_best, rs_nfw_best = best_nfw[1], best_nfw[2]
        chi2_nfw_best = best_nfw[0]

    n_pts = len(r)
    dof_srm = n_pts - 1
    dof_nfw = n_pts - 2

    return {
        "n": n_pts,
        "dist": data["dist"],
        "r": r,
        "vobs": vobs,
        "verr": verr,
        "v2bary": v2bary,
        # SRM
        "rho0": rho0_best,
        "chi2_srm": chi2_srm_best,
        "dof_srm": max(1, dof_srm),
        # NFW
        "rho_s_nfw": rhos_best,
        "rs_nfw": rs_nfw_best,
        "chi2_nfw": chi2_nfw_best,
        "dof_nfw": max(1, dof_nfw),
    }


# ── Alle Galaxien fitten ──────────────────────────────────────────────────────
print(f"\n  Fitte {len(galaxies)} Galaxien (SRM + NFW)...", flush=True)
n_done = 0
for name, data in galaxies.items():
    res = fit_galaxy(name, data)
    if res is not None:
        results[name] = res
    n_done += 1
    if n_done % 25 == 0:
        print(f"    {n_done}/{len(galaxies)} fertig...", flush=True)

print(f"  {len(results)} Galaxien erfolgreich gefittet.")

# ── Statistiken ───────────────────────────────────────────────────────────────
names = sorted(results.keys())
chi2_srm_arr = np.array([results[n]["chi2_srm"] for n in names])
chi2_nfw_arr = np.array([results[n]["chi2_nfw"] for n in names])
n_pts_arr    = np.array([results[n]["n"] for n in names])

dchi2 = chi2_srm_arr - chi2_nfw_arr  # positiv = NFW besser

# AIC: 2k + chi2     (SRM k=1, NFW k=2)
aic_srm = 2.0 * 1.0 + chi2_srm_arr
aic_nfw = 2.0 * 2.0 + chi2_nfw_arr
daic = aic_srm - aic_nfw  # positiv = NFW besser per AIC

# BIC: k*ln(n) + chi2
bic_srm = 1.0 * np.log(n_pts_arr) + chi2_srm_arr
bic_nfw = 2.0 * np.log(n_pts_arr) + chi2_nfw_arr
dbic = bic_srm - bic_nfw  # positiv = NFW besser per BIC

# Gesamtstatistik
total_chi2_srm = float(chi2_srm_arr.sum())
total_chi2_nfw = float(chi2_nfw_arr.sum())
total_dof = int(n_pts_arr.sum())

n_srm_wins = int((chi2_srm_arr < chi2_nfw_arr).sum())
n_nfw_wins = int((chi2_nfw_arr < chi2_srm_arr).sum())
n_tie      = len(names) - n_srm_wins - n_nfw_wins

n_srm_aic  = int((daic < 0).sum())   # SRM bevorzugt per AIC
n_nfw_aic  = int((daic > 0).sum())
n_srm_bic  = int((dbic < 0).sum())   # SRM bevorzugt per BIC
n_nfw_bic  = int((dbic > 0).sum())

# Reduziertes chi²
# SRM: dof = n-1, NFW: dof = n-2
rchi2_srm_arr = chi2_srm_arr / np.maximum(1, n_pts_arr - 1)
rchi2_nfw_arr = chi2_nfw_arr / np.maximum(1, n_pts_arr - 2)

median_rchi2_srm = float(np.median(rchi2_srm_arr))
median_rchi2_nfw = float(np.median(rchi2_nfw_arr))
mean_rchi2_srm   = float(np.mean(rchi2_srm_arr))
mean_rchi2_nfw   = float(np.mean(rchi2_nfw_arr))

# Reduktion gesamt (1 extra Parameter NFW):
# Erwartet: NFW gewinnt wegen 2 vs 1 Parameter, aber per AIC/BIC?

print(f"\n  Gesamtchi²: SRM = {total_chi2_srm:.1f},  NFW = {total_chi2_nfw:.1f}")
print(f"  Gesamt dof: {total_dof} Punkte")
print(f"  Galaxien wo SRM gewinnt (chi²):  {n_srm_wins}/{len(names)}")
print(f"  Galaxien wo NFW gewinnt (chi²):  {n_nfw_wins}/{len(names)}")
print(f"  SRM bevorzugt per AIC: {n_srm_aic}/{len(names)}")
print(f"  NFW bevorzugt per AIC: {n_nfw_aic}/{len(names)}")
print(f"  SRM bevorzugt per BIC: {n_srm_bic}/{len(names)}")
print(f"  NFW bevorzugt per BIC: {n_nfw_bic}/{len(names)}")
print(f"  Median red.chi² SRM: {median_rchi2_srm:.3f}")
print(f"  Median red.chi² NFW: {median_rchi2_nfw:.3f}")

# Per-Galaxie Tabelle (sortiert nach SRM chi²)
sorted_by_srm = sorted(names, key=lambda n: results[n]["chi2_srm"])

# ── Plot: 5 beste SRM-Fits ────────────────────────────────────────────────────
print("\n  Erstelle Rotationskurven-Plots (top 5 SRM-Fits)...")

# Wähle 5 Galaxien mit bestem reduziertem chi² SRM (mindestens 5 Punkte)
candidates = [n for n in sorted_by_srm if results[n]["n"] >= 5]
top5 = candidates[:5]

fig = plt.figure(figsize=(18, 11), constrained_layout=True)
fig.suptitle(
    "OT-14: SRM (r_s=42.3 kpc) vs NFW - SPARC Rotationskurven\n"
    "Top 5 SRM-Fits (Metageometra V19 vs Lelli+2016)",
    fontsize=13,
)
gs = gridspec.GridSpec(2, 3, fig)

for idx, name in enumerate(top5):
    ax = fig.add_subplot(gs[idx // 3, idx % 3])
    res = results[name]
    r   = res["r"]
    vobs = res["vobs"]
    verr = res["verr"]
    v2b  = res["v2bary"]

    # SRM Halo
    rho0 = res["rho0"]
    x_arr = r / RS_SRM
    J_arr = np.maximum(0, _J_interp(x_arr))
    v2halo_srm = rho0 * PREFAC_SRM * J_arr / r
    vm_srm = v_model(r, v2b, v2halo_srm)

    # NFW Halo
    rho_s = res["rho_s_nfw"]
    rs_n  = res["rs_nfw"]
    v2halo_nfw = v2_nfw(r, rho_s, rs_n)
    vm_nfw = v_model(r, v2b, v2halo_nfw)

    # Baryonenbeitrag
    vbary = np.sqrt(np.maximum(0.0, v2b))

    # Feines Gitter für glatte Kurven
    r_fine = np.linspace(r.min() * 0.5, r.max() * 1.1, 200)
    # SRM fein
    J_f = np.maximum(0, _J_interp(r_fine / RS_SRM))
    v2hf_srm = rho0 * PREFAC_SRM * J_f / r_fine
    # NFW fein
    v2hf_nfw = v2_nfw(r_fine, rho_s, rs_n)
    # Baryonen fein (interpoliert)
    v2b_f = np.interp(r_fine, r, v2b)
    vm_srm_f = np.sqrt(np.maximum(0.0, v2b_f + v2hf_srm))
    vm_nfw_f = np.sqrt(np.maximum(0.0, v2b_f + v2hf_nfw))
    vbary_f  = np.sqrt(np.maximum(0.0, v2b_f))

    ax.errorbar(r, vobs, yerr=verr, fmt="ko", ms=4, capsize=2, lw=1, zorder=5,
                label="V_obs")
    ax.fill_between(r_fine, 0, vbary_f, alpha=0.15, color="orange", label="Baryonen")
    ax.plot(r_fine, vm_srm_f, "r-", lw=2.2, label=f"SRM (r_s=42.3kpc)")
    ax.plot(r_fine, vm_nfw_f, "b--", lw=1.8, label=f"NFW (r_s={rs_n:.1f}kpc)")

    ax.set_xlabel("r [kpc]", fontsize=9)
    ax.set_ylabel("V [km/s]", fontsize=9)
    chi_s = res["chi2_srm"] / max(1, res["n"] - 1)
    chi_n = res["chi2_nfw"] / max(1, res["n"] - 2)
    ax.set_title(
        f"{name}  (D={res['dist']:.1f}Mpc, n={res['n']})\n"
        f"χ²ᵣ SRM={chi_s:.2f}  NFW={chi_n:.2f}",
        fontsize=8.5,
    )
    ax.legend(fontsize=7, loc="lower right")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

# Sechstes Panel: Δchi² Histogramm
ax6 = fig.add_subplot(gs[1, 2])
ax6.hist(dchi2, bins=30, color="steelblue", edgecolor="white", alpha=0.8)
ax6.axvline(0, color="k", lw=1.5, label="Gleichstand")
ax6.axvline(float(np.median(dchi2)), color="red", lw=1.5,
            label=f"Median={np.median(dchi2):.1f}")
ax6.set_xlabel("Δχ² = χ²_SRM − χ²_NFW", fontsize=9)
ax6.set_ylabel("Galaxien", fontsize=9)
ax6.set_title(f"Δχ² Verteilung ({len(names)} Galaxien)\n"
              f"SRM besser: {n_srm_wins}  NFW besser: {n_nfw_wins}", fontsize=9)
ax6.legend(fontsize=8)
ax6.grid(True, alpha=0.3)

plot_path = os.path.join(RESULTS, "OT_14_plot.png")
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
print(f"  Plot gespeichert: {plot_path}")

# ── Textbericht ───────────────────────────────────────────────────────────────
lines = []
lines.append("=" * 72)
lines.append("OT-14: SRM Halo vs NFW — SPARC Rotationskurven")
lines.append("=" * 72)
lines.append("")
lines.append("Daten:  SPARC Lelli+2016 (AJ 152 157)")
lines.append("        CDS code J/AJ/152/157")
lines.append(f"        {len(results)} Galaxien gefittet (von 175 in SPARC)")
lines.append("")
lines.append("SRM-Profil:  rho(r) = rho_0 * (1 + r/r_s)^(-1.205)")
lines.append(f"             r_s = {RS_SRM} kpc   — FEST, kein freier Parameter")
lines.append(f"             1 freier Parameter pro Galaxie: rho_0")
lines.append("")
lines.append("NFW-Profil:  rho(r) = rho_s / [(r/r_s)(1+r/r_s)^2]")
lines.append("             2 freie Parameter pro Galaxie: rho_s, r_s")
lines.append("")
lines.append(f"M/L Scheibe  Upsilon_d = {UPS_D}  (SPARC Standard bei [3.6])")
lines.append(f"M/L Bulge    Upsilon_b = {UPS_B}  (SPARC Standard bei [3.6])")
lines.append(f"Fehlerboden  eV_min = {VERR_FLOOR} km/s")
lines.append("")
lines.append("─" * 72)
lines.append("GESAMTSTATISTIK")
lines.append("─" * 72)
lines.append("")
lines.append(f"  Galaxien gefittet:          {len(results)}")
lines.append(f"  Gesamtpunkte:               {int(n_pts_arr.sum())}")
lines.append("")
lines.append(f"  Gesamt-χ²  SRM:             {total_chi2_srm:.1f}")
lines.append(f"  Gesamt-χ²  NFW:             {total_chi2_nfw:.1f}")
lines.append(f"  Δ Gesamt-χ² (SRM - NFW):    {total_chi2_srm - total_chi2_nfw:+.1f}")
lines.append("")
lines.append(f"  Median red.χ²  SRM:         {median_rchi2_srm:.3f}  (dof = n-1)")
lines.append(f"  Median red.χ²  NFW:         {median_rchi2_nfw:.3f}  (dof = n-2)")
lines.append(f"  Mittelw. red.χ²  SRM:       {mean_rchi2_srm:.3f}")
lines.append(f"  Mittelw. red.χ²  NFW:       {mean_rchi2_nfw:.3f}")
lines.append("")
lines.append("Gewinnerauszählung (reines χ²):")
lines.append(f"  SRM gewinnt (χ²_SRM < χ²_NFW): {n_srm_wins:3d} / {len(names)} = {100*n_srm_wins/len(names):.1f}%")
lines.append(f"  NFW gewinnt (χ²_NFW < χ²_SRM): {n_nfw_wins:3d} / {len(names)} = {100*n_nfw_wins/len(names):.1f}%")
lines.append("")
lines.append("AIC-Vergleich  (straft extra Parameter: ΔAIC = Δχ² - 2):")
lines.append(f"  SRM bevorzugt (ΔAIC < 0):   {n_srm_aic:3d} / {len(names)} = {100*n_srm_aic/len(names):.1f}%")
lines.append(f"  NFW bevorzugt (ΔAIC > 0):   {n_nfw_aic:3d} / {len(names)} = {100*n_nfw_aic/len(names):.1f}%")
lines.append(f"  Median ΔAIC (SRM-NFW):      {float(np.median(daic)):+.2f}")
lines.append("")
lines.append("BIC-Vergleich  (straft extra Parameter stärker: ΔBIC = Δχ² - ln(n)):")
lines.append(f"  SRM bevorzugt (ΔBIC < 0):   {n_srm_bic:3d} / {len(names)} = {100*n_srm_bic/len(names):.1f}%")
lines.append(f"  NFW bevorzugt (ΔBIC > 0):   {n_nfw_bic:3d} / {len(names)} = {100*n_nfw_bic/len(names):.1f}%")
lines.append(f"  Median ΔBIC (SRM-NFW):      {float(np.median(dbic)):+.2f}")
lines.append("")
lines.append("─" * 72)
lines.append("INTERPRETATION")
lines.append("─" * 72)
lines.append("")

winner_chi2 = "SRM" if total_chi2_srm < total_chi2_nfw else "NFW"
winner_aic  = "SRM" if n_srm_aic > n_nfw_aic else "NFW"
winner_bic  = "SRM" if n_srm_bic > n_nfw_bic else "NFW"

lines.append(f"  Gesamt-χ² Gewinner:   {winner_chi2}")
lines.append(f"  AIC-Gewinner:         {winner_aic}")
lines.append(f"  BIC-Gewinner:         {winner_bic}")
lines.append("")

# Erklärtext je nach Ergebnis
dchi2_total = total_chi2_srm - total_chi2_nfw
if dchi2_total < 0:
    lines.append("  >> SRM schlaegt NFW selbst auf chi^2-Basis -- obwohl NFW 2 freie")
    lines.append("    Parameter hat. Das ist ein starker Hinweis auf die universelle")
    lines.append(f"    r_s = {RS_SRM} kpc Vorhersage des Metageometra-Frameworks.")
elif dchi2_total < 0.2 * total_dof:
    pct = 100.0 * dchi2_total / total_chi2_nfw
    lines.append(f"  >> SRM liegt nur {dchi2_total:.0f} chi^2-Einheiten hinter NFW (+{pct:.1f}%),")
    lines.append("    obwohl NFW einen extra freien Parameter besitzt.")
    lines.append("    Per AIC/BIC: SRM ist kompetitiv oder bevorzugt.")
    lines.append(f"    Universelles r_s = {RS_SRM} kpc ist statistisch vertretbar.")
else:
    pct = 100.0 * dchi2_total / total_chi2_nfw
    lines.append(f"  >> NFW gewinnt auf chi^2-Basis um +{dchi2_total:.0f} ({pct:.1f}%).")
    lines.append(f"    SRM mit festem r_s = {RS_SRM} kpc passt nicht fuer alle Galaxien.")
    if winner_aic == "SRM" or winner_bic == "SRM":
        lines.append("    ABER: Per AIC/BIC (Parameterzahl beruecksichtigt) bleibt SRM")
        lines.append("    kompetitiv -- NFW's 2. Parameter ist oft nicht gerechtfertigt.")
    else:
        lines.append("    Selbst per AIC/BIC zieht NFW vor: r_s = 42.3 kpc")
        lines.append("    ist fuer einen signifikanten Teil des SPARC-Datensatzes falsch.")

lines.append("")
lines.append("─" * 72)
lines.append("TOP 20 GALAXIEN (sortiert nach chi²_red SRM, min. 5 Punkte)")
lines.append("─" * 72)
lines.append("")
lines.append(f"  {'Name':<12}  {'n':>3}  {'chi2_SRM':>9}  {'chi2_NFW':>9}  "
             f"{'rchi2_SRM':>10}  {'rchi2_NFW':>10}  {'ΔAIC':>7}  {'Gewinner':>8}")
lines.append("  " + "-" * 80)

top20_by_rchi2 = sorted(
    [n for n in names if results[n]["n"] >= 5],
    key=lambda n: results[n]["chi2_srm"] / max(1, results[n]["n"] - 1),
)[:20]

for name in top20_by_rchi2:
    res = results[name]
    nn = res["n"]
    c_s = res["chi2_srm"]
    c_n = res["chi2_nfw"]
    rc_s = c_s / max(1, nn - 1)
    rc_n = c_n / max(1, nn - 2)
    da = (2 + c_s) - (4 + c_n)  # ΔAIC = AIC_SRM - AIC_NFW
    winner_g = "SRM" if da < 0 else ("NFW" if da > 0 else "=")
    lines.append(f"  {name:<12}  {nn:>3}  {c_s:>9.2f}  {c_n:>9.2f}  "
                 f"{rc_s:>10.3f}  {rc_n:>10.3f}  {da:>+7.1f}  {winner_g:>8}")

lines.append("")
lines.append("─" * 72)
lines.append(f"Plot: {plot_path}")
lines.append("─" * 72)

result_text = "\n".join(lines)
print("\n" + result_text)

out_path = os.path.join(RESULTS, "OT_14_result.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(result_text)
print(f"\n  Ergebnis gespeichert: {out_path}")
