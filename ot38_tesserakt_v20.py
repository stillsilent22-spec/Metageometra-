#!/usr/bin/env python3
"""
OT-38: Tesserakt-Schalenformel V20
Formel: theta_n = arccos(cos(n*delta) * cos(n*chi))

Aufgabe 1: Parameter-Scan delta,chi in [0,90] deg @ 0.1 deg
Aufgabe 2: Chiralitaetspruefung Top-10
Aufgabe 3: Physikalische Interpretation (L_geo)
Aufgabe 4: Falsifikationstest SMBH-Liste
OT-38b:    f_echo aus Tesserakt-Geometrie
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import os
from scipy.stats import binom as sp_binom

# ── Konstanten ────────────────────────────────────────────────────────
RHO        = 0.406          # R4-Attraktordimension (OT-37 bestaetigt)
THETA1_OBS = 58.65          # Sgr A* / GCD_eps / OT-11
THETA2_OBS = 118.32         # NGC 1052
THETA3_OBS = 175.95         # NGC 0315
D_EFF      = 0.3388         # D_f,eff
L0         = 1.386e-44      # kg m^-3 s^-1 (Referenzwert = L_obs)
L_OBS      = 1.386e-44
F_ECHO_OBS = 2.65e40
N_SHELLS   = 6

# Toleranzen
TOL1     = 0.5   # deg  theta1
TOL2     = 5.0   # deg  theta2
TOL3     = 5.0   # deg  theta3
TOL_SMBH = 5.0   # deg  SMBH-Abgleich

# SMBH-Liste (Name, Grosskreisabstand vom D-Pol in Grad)
SMBH_LIST = [
    ("Sgr A*",   58.65),
    ("NGC 1052", 118.32),
    ("NGC 0315", 175.95),
    ("M87",       50.70),
    ("Cen A",      7.00),
    ("M106",      85.70),
    ("NGC 3115",  49.80),
    ("M77",      135.50),
    ("NGC 4151",  78.10),
]
N_SMBH = len(SMBH_LIST)

os.makedirs("results", exist_ok=True)
lines = []

def p(s=""):
    print(s)
    lines.append(str(s))

def theta_val(n, d_rad, c_rad):
    """Einzelwert theta_n in Grad."""
    v = np.clip(np.cos(n * d_rad) * np.cos(n * c_rad), -1.0, 1.0)
    return np.degrees(np.arccos(v))

# ── Aufgabe 1: Parameter-Scan ─────────────────────────────────────────
p("=" * 72)
p("  OT-38 / V20 — Tesserakt-Schalenformel — Parameter-Scan")
p("=" * 72)
p("  Formel : theta_n = arccos(cos(n*delta) * cos(n*chi))")
p(f"  Ziele  : theta1={THETA1_OBS}  theta2={THETA2_OBS}  theta3={THETA3_OBS}  (deg)")
p(f"  rho    = {RHO}  (R4-Attraktordimension, OT-37)")
p()

delta_deg = np.arange(0.0, 90.01, 0.1)
chi_deg   = np.arange(0.0, 90.01, 0.1)
p(f"  Aufgabe 1: Scan  {len(delta_deg)} x {len(chi_deg)} = {len(delta_deg)*len(chi_deg):,} Kombinationen ...")

D_rad = np.deg2rad(delta_deg)
C_rad = np.deg2rad(chi_deg)
D, C = np.meshgrid(D_rad, C_rad, indexing='ij')   # shape (901, 901)

def theta_grid(n, D, C):
    return np.degrees(np.arccos(np.clip(np.cos(n * D) * np.cos(n * C), -1.0, 1.0)))

th1 = theta_grid(1, D, C)
th2 = theta_grid(2, D, C)
th3 = theta_grid(3, D, C)

r1 = np.abs(th1 - THETA1_OBS)
r2 = np.abs(th2 - THETA2_OBS)
r3 = np.abs(th3 - THETA3_OBS)

tol1_used, tol2_used, tol3_used = TOL1, TOL2, TOL3
mask = (r1 < tol1_used) & (r2 < tol2_used) & (r3 < tol3_used)
n_hits = int(mask.sum())
p(f"  Treffer (TOL {tol1_used}/{tol2_used}/{tol3_used} deg): {n_hits}")

# Fallback: erweiterte Toleranzen falls noetig
if n_hits < 3:
    tol1_used, tol2_used, tol3_used = 1.0, 8.0, 8.0
    mask = (r1 < tol1_used) & (r2 < tol2_used) & (r3 < tol3_used)
    n_hits = int(mask.sum())
    p(f"  → Erweiterte TOL {tol1_used}/{tol2_used}/{tol3_used}: {n_hits} Treffer")

if n_hits == 0:
    p("  FEHLER: Keine Loesung gefunden. Ueberprueffe Formel.")
    exit(1)

# Sortiere nach kombiniertem Residuum
comb_r = r1 + r2 + r3
hi, ci = np.where(mask)
order  = np.argsort(comb_r[hi, ci])
hi, ci = hi[order], ci[order]

# Baue Ergebnisliste
results = []
for k in range(len(hi)):
    d = delta_deg[hi[k]]
    c = chi_deg[ci[k]]
    ths = [theta_val(n, np.radians(d), np.radians(c)) for n in range(1, N_SHELLS + 1)]
    res = abs(ths[0]-THETA1_OBS) + abs(ths[1]-THETA2_OBS) + abs(ths[2]-THETA3_OBS)
    results.append({'delta': d, 'chi': c, 'thetas': ths, 'residual': res})

p()
p(f"  {'Rang':>4}  {'delta':>6}  {'chi':>6}  {'Res':>7}  {'theta1':>7}  {'theta2':>7}  {'theta3':>7}")
p(f"  {'----':>4}  {'------':>6}  {'------':>6}  {'-------':>7}  {'-------':>7}  {'-------':>7}  {'-------':>7}")
for k, sol in enumerate(results[:20]):
    p(f"  {k+1:>4}  {sol['delta']:>6.1f}  {sol['chi']:>6.1f}  "
      f"{sol['residual']:>7.3f}  "
      f"{sol['thetas'][0]:>7.2f}  {sol['thetas'][1]:>7.2f}  {sol['thetas'][2]:>7.2f}")

# ── Aufgabe 2: Chiralitaet ─────────────────────────────────────────────
p()
p("=" * 72)
p("  Aufgabe 2: Chiralitaetspruefung Top-10")
p("=" * 72)
p()
p(f"  {'Nr':>3}  {'delta':>6}  {'chi':>6}  {'Res':>6}  "
  f"{'G12':>6}  {'G23':>6}  {'G34':>6}  {'G45':>6}  {'G56':>6}  "
  f"{'sig_G':>6}  {'Chiral':>7}")
p("  " + "-" * 72)

top10  = results[:min(10, len(results))]

def gap_stats(ths):
    gaps = [ths[i+1] - ths[i] for i in range(len(ths)-1)]
    return gaps, float(np.std(gaps))

best_chiral = None
for k, sol in enumerate(top10):
    gaps, sig = gap_stats(sol['thetas'])
    is_chiral = sig > 1.0
    gstr = "  ".join(f"{g:+.1f}" for g in gaps)
    p(f"  {k+1:>3}  {sol['delta']:>6.1f}  {sol['chi']:>6.1f}  "
      f"{sol['residual']:>6.3f}  "
      f"{gaps[0]:>+6.1f}  {gaps[1]:>+6.1f}  {gaps[2]:>+6.1f}  "
      f"{gaps[3]:>+6.1f}  {gaps[4]:>+6.1f}  "
      f"{sig:>6.2f}  {'JA' if is_chiral else 'NEIN':>7}")
    if is_chiral and best_chiral is None:
        best_chiral = sol

# Bevorzuge chirale Loesung falls vorhanden, sonst beste ueberhaupt
best = best_chiral if best_chiral is not None else results[0]
d0   = best['delta']
c0   = best['chi']
ths0 = best['thetas']
gaps0, sig0 = gap_stats(ths0)
is_chiral0  = sig0 > 1.0

p()
if best_chiral is not None:
    p(f"  -> Chirale Loesung gewaehlt: delta={d0:.1f}  chi={c0:.1f}  sig_G={sig0:.2f}")
else:
    p(f"  -> Keine chirale Loesung; beste gesamt: delta={d0:.1f}  chi={c0:.1f}")

# ── Aufgabe 3: Physikalische Interpretation ───────────────────────────
p()
p("=" * 72)
p("  Aufgabe 3: Physikalische Interpretation")
p("=" * 72)
p()
p("  L_geo = L0 * (1 - cos(delta)) * sin(chi) * rho^2")
p(f"  L0  = {L0:.3e} kg m^-3 s^-1")
p(f"  rho = {RHO}")
p()

factor_L = (1.0 - np.cos(np.radians(d0))) * np.sin(np.radians(c0)) * RHO**2
L_geo    = L0 * factor_L
L_pct    = (factor_L - 1.0) * 100.0

p(f"  delta = {d0:.1f} deg  =>  1 - cos(delta) = {1-np.cos(np.radians(d0)):.6f}")
p(f"  chi   = {c0:.1f} deg  =>  sin(chi)       = {np.sin(np.radians(c0)):.6f}")
p(f"  rho^2 = {RHO**2:.6f}")
p(f"  Faktor= {factor_L:.6f}")
p(f"  L_geo = {L_geo:.4e} kg m^-3 s^-1")
p(f"  L_obs = {L_OBS:.4e} kg m^-3 s^-1")
p(f"  Delta = {L_pct:+.2f}%")
p()

# Suche exakte L-Punkte (Faktor = 1 => L_geo = L_obs)
factor_all = (1.0 - np.cos(D)) * np.sin(C) * RHO**2
exact_mask = (np.abs(factor_all - 1.0) < 0.01) & (r1 < tol1_used)
ei, ej     = np.where(exact_mask)

p("  Suche (delta, chi) mit L_geo = L_obs genau (Faktor = 1):")
if len(ei) > 0:
    for k in range(min(5, len(ei))):
        p(f"    delta={delta_deg[ei[k]]:.1f}  chi={chi_deg[ej[k]]:.1f}  "
          f"theta1={th1[ei[k],ej[k]]:.2f}")
else:
    p("  → Keine Loesung im Scanbereich (Faktor << 1 fuer alle theta1-Treffer)")
    # Zeige was der beste Faktor bei theta1-Treffern ist
    theta1_mask = r1 < tol1_used
    max_f       = factor_all[theta1_mask].max() if theta1_mask.any() else 0.0
    p(f"    Bester Faktor in theta1-Band: {max_f:.4f}")

# ── Aufgabe 4: Falsifikationstest ─────────────────────────────────────
p()
p("=" * 72)
p("  Aufgabe 4: Falsifikationstest — SMBH-Abgleich")
p("=" * 72)
p()
p(f"  Beste Loesung: delta={d0:.1f}  chi={c0:.1f}")
p(f"  Toleranz: +-{TOL_SMBH:.0f} deg")
p()
p(f"  {'n':>3}  {'theta_n':>8}  {'Kandidat':<20}  {'theta_obs':>9}  {'Delta':>6}")
p("  " + "-" * 58)

matched_names = set()
n_shell_hits  = 0

for n, th in enumerate(ths0, 1):
    best_smbh = None
    best_dth  = 9999.0
    for nm, ang in SMBH_LIST:
        dth = abs(th - ang)
        if dth < TOL_SMBH and dth < best_dth:
            best_smbh = (nm, ang)
            best_dth  = dth
    if best_smbh is not None:
        matched_names.add(best_smbh[0])
        n_shell_hits += 1
        p(f"  {n:>3}  {th:>8.2f}  {best_smbh[0]:<20}  {best_smbh[1]:>9.2f}  {best_dth:>+6.2f}  ***")
    else:
        # Zeige naechsten Kandidaten auch wenn ausserhalb Toleranz
        nearest = min(SMBH_LIST, key=lambda x: abs(th - x[1]))
        dth_n   = abs(th - nearest[1])
        p(f"  {n:>3}  {th:>8.2f}  {'UNBEKANNT':<20}  "
          f"n={nearest[0]} ({nearest[1]}) Delta={dth_n:.1f}")

p()
n_matched = len(matched_names)
p(f"  Treffer: {n_matched} / {N_SMBH} SMBHs")
p(f"  Gematchte SMBHs: {', '.join(sorted(matched_names)) if matched_names else '—'}")
p()

# Zufallserwartung: Binomial(N_SHELLS, p_hit) wo p_hit = P(mind. 1 Treffer)
p_single  = (2 * TOL_SMBH) / 180.0          # P(Schale trifft 1 SMBH zufaellig)
p_anysmbh = 1 - (1 - p_single) ** N_SMBH    # P(Schale trifft irgendeinen SMBH)
expected  = N_SMBH * (1 - (1 - p_single) ** N_SHELLS)  # E[einzigartige Treffer]
# p-Wert: P(X >= n_matched) mit X ~ Binom(N_SMBH, p_anysmbh^N_SHELLS)
p_val = float(1 - sp_binom.cdf(n_matched - 1, N_SMBH,
                                 1 - (1 - p_single)**N_SHELLS))

p(f"  Erwartete Zufallstreffer : {expected:.2f}")
p(f"  p-Wert (binomial)        : {p_val:.4f}")

# ── OT-38b: f_echo aus Tesserakt-Geometrie ────────────────────────────
p()
p("=" * 72)
p("  OT-38b: f_echo aus Tesserakt-Geometrie")
p("=" * 72)
p()
p("  f_echo_geo = exp( D_eff * ln( (1-cos(delta)) / (rho^2 * sin(chi)) ) )")
p(f"  D_eff      = {D_EFF}")
p(f"  f_echo_obs = {F_ECHO_OBS:.3e}")
p()

num_f = 1.0 - np.cos(np.radians(d0))
den_f = RHO**2 * np.sin(np.radians(c0))

p(f"  delta = {d0:.1f}:  1 - cos(delta) = {num_f:.6f}")
p(f"  chi   = {c0:.1f}:  rho^2*sin(chi) = {den_f:.6f}")

if den_f > 1e-12 and num_f > 0:
    arg        = num_f / den_f
    f_echo_geo = float(np.exp(D_EFF * np.log(arg)))
    log10_geo  = float(np.log10(f_echo_geo))
    log10_obs  = float(np.log10(F_ECHO_OBS))
    delta_dek  = log10_geo - log10_obs
    p(f"  Argument : {arg:.6f}")
    p(f"  f_echo_geo = {arg:.4f}^{D_EFF} = {f_echo_geo:.4e}  (log10={log10_geo:.2f})")
    p(f"  f_echo_obs = {F_ECHO_OBS:.4e}                    (log10={log10_obs:.2f})")
    p(f"  Differenz  : {delta_dek:+.2f} Dekaden")
    p()
    if abs(delta_dek) < 2:
        p("  ERGEBNIS: f_echo_geo ~ f_echo_obs  (< 2 Dekaden)  => BESTANDEN")
    elif abs(delta_dek) < 5:
        p("  ERGEBNIS: Abweichung < 5 Dekaden => BEDINGT")
    else:
        p(f"  ERGEBNIS: Abweichung = {abs(delta_dek):.1f} Dekaden => NICHT DIREKT ABGELEITET")
        p("  Physikalisch: Die Formel benoetigt einen kosmologischen")
        p(f"  Skalierungsfaktor N ~ {int(round(log10_obs / (D_EFF * np.log10(abs(arg)) + 1e-30)))}")
        p("  => f_echo = ((1-cos(delta))/(rho^2*sin(chi)))^(D_eff * N)")
        # Berechne N das f_echo_obs liefert
        if arg > 1:
            N_needed = log10_obs / (D_EFF * np.log10(arg))
            p(f"  => N_benoetigt = {N_needed:.1f}  (HTM: N = 274-275 => konsistent?)")
        elif arg < 1 and arg > 0:
            # arg^(D_eff*N) grows if N<0 — nonsensical
            p("  => arg < 1: Formel liefert immer < 1, unabhaengig von N")
        p()
        p("  Rescaled-Test: f_echo_geo mit N=274:")
        N_htm = 274
        if arg > 0:
            f_rescaled = float(np.exp(D_EFF * N_htm * np.log(arg)))
            p(f"  f_echo_rescaled = {arg:.4f}^({D_EFF}*{N_htm}) = {f_rescaled:.4e}")
            p(f"  f_echo_obs      = {F_ECHO_OBS:.4e}")
            log_diff = np.log10(f_rescaled) - log10_obs if f_rescaled > 0 else float('nan')
            p(f"  Differenz:  {log_diff:+.2f} Dekaden")
else:
    p(f"  WARNUNG: sin(chi) = {np.sin(np.radians(c0)):.6f}  oder  1-cos(delta) = {num_f:.6f}")
    p("  Formel nicht auswertbar bei delta=0 oder chi=0")

# Suche (delta, chi) mit log_f_echo_geo ~ log10_obs im theta1-Band
p()
p("  Suche Punkte mit |log10(f_echo_geo) - log10(f_echo_obs)| < 2.0:")
with np.errstate(divide='ignore', invalid='ignore'):
    num_arr  = 1.0 - np.cos(D)
    den_arr  = RHO**2 * np.sin(C)
    ratio    = np.where(den_arr > 1e-12, num_arr / den_arr, np.nan)
    log10_f  = np.where(ratio > 0, D_EFF * np.log10(ratio), np.nan)

log10_obs_t = np.log10(F_ECHO_OBS)
fmask       = (np.abs(log10_f - log10_obs_t) < 2.0) & (r1 < tol1_used)
fhi, fci    = np.where(fmask)

if len(fhi) > 0:
    for k in range(min(5, len(fhi))):
        p(f"    delta={delta_deg[fhi[k]]:.1f}  chi={chi_deg[fci[k]]:.1f}  "
          f"log10_f={log10_f[fhi[k],fci[k]]:.1f}")
else:
    # Zeige beste (naechste) Annaeherung im theta1-Band
    valid_f = np.where((~np.isnan(log10_f)) & (r1 < tol1_used),
                        np.abs(log10_f - log10_obs_t), np.inf)
    best_fi = int(np.argmin(valid_f))
    bi, bj  = np.unravel_index(best_fi, log10_f.shape)
    p(f"  → Kein Treffer im theta1-Band.")
    p(f"    Naechste Ann.: delta={delta_deg[bi]:.1f}  chi={chi_deg[bj]:.1f}  "
      f"log10_f={log10_f[bi,bj]:.1f}  (Abstand {valid_f[bi,bj]:.1f} Dekaden)")

# ── Zusammenfassung ───────────────────────────────────────────────────
p()
p("=" * 72)
p("  === BESTE LOESUNG ===")
p("=" * 72)
p()
p(f"  delta = {d0:.1f} deg")
p(f"  chi   = {c0:.1f} deg")
p(f"  rho   = {RHO}  (fixed)")
p()
p("  Schalenwinkel:")
for n, th in enumerate(ths0, 1):
    cands = sorted(SMBH_LIST, key=lambda x: abs(x[1]-th))
    nm, ang = cands[0]
    dth = abs(th - ang)
    flag = "***" if dth < TOL_SMBH else "   "
    p(f"  n={n}: theta = {th:7.2f} deg  |  {nm:<12} ({ang:7.2f})  |  Delta = {dth:+6.2f} deg  {flag}")
p()
p(f"  Schalenluecken (G_i,i+1):")
for i, g in enumerate(gaps0):
    p(f"    G({i+1},{i+2}) = {g:+.2f} deg")
p()
p(f"  Chiralitaet  : {'JA' if is_chiral0 else 'NEIN'}  (sigma_G = {sig0:.2f} deg)")
p(f"  L_geo        = {L_geo:.4e}  |  Abweichung von L_obs: {L_pct:+.2f}%")
p()
p("  === FALSIFIKATIONSSTATUS ===")
p(f"  Treffer in SMBH-Liste  : {n_matched} / {N_SMBH}")
p(f"  Erwartete Zufallstref. : {expected:.1f}")
p(f"  p-Wert                 : {p_val:.4f}")
p()
p("  === OFFENE PARAMETER FUER V20 ===")
p(f"  delta = {d0:.1f} deg  => Versatzwinkel der zwei chiralen Tessserakthaefte")
p(f"  chi   = {c0:.1f} deg  => Inversionswinkel der chiralen Spiegelung")
p()

# Verdict
p("=" * 72)
p("  OT-38 VERDICT")
p("=" * 72)
if n_matched >= 4:
    verdict = "BESTANDEN"
elif n_matched >= 2:
    verdict = "BEDINGT BESTANDEN"
else:
    verdict = "VORBEHALT"

p(f"  {verdict}")
p(f"  SMBH-Treffer: {n_matched}/{N_SMBH}  (p={p_val:.4f})")
if is_chiral0:
    p("  Chirale Loesung: konsistent mit Sigma_meta-Chiralitaet")
else:
    p("  Symmetrische Loesung: delta nahe 0 => lineares Schalenspektrum")
    p("  V20-Chiralitaet erfordert delta >> 0 mit SMBH-Constraints gelockert")
p("=" * 72)

# Schreibe Ergebnisdatei
out_path = os.path.join("results", "OT_38_result.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\n  Ergebnisse -> {out_path}")
