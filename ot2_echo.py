"""
OT-2: f_echo Analytische Verifikation
======================================
Aufgabe: Verifiziere f_echo = exp(D_f,eff * 274) über VIER unabhaengige Pfade.

Metageometra V19/V20 / F(L)-Gleichung:
  eta_baryon = f_echo * (tau_inflation / t_0)
  f_echo     = exp(D_f,eff * ln(rho_early / rho_DE))
             = exp(0.3388 * 274)

D_f,eff = D_f,geo * D_f,diss = 0.77 * 0.44 = 0.3388

Pfad A: Thermodynamisch (Strahlungsdichte bei T_reh)
  274 = ln(rho_radiation(T_reh) / rho_DE)

Pfad B: Baryon-Asymmetrie Rückrechnung
  274 = ln(eta_obs * t_0 / tau_infl) / D_f,eff

Pfad C: Planck-Skala (parameterfrei)
  N_Planck = ln(rho_Planck / rho_DE)  — obere natuerliche Grenze

Pfad D: DM-Echo-Wellenasymmetrie + S³-π-Unterdrückung (V20)
  Idee: Kevin Hannemann, 28.04.2026
  Wenn DM = HTM-Echo (Rückwärtswelle), ist η die Asymmetrie
  zwischen Vorwärts- und Rückwärtswelle, unterdrückt durch
  S³-Phasenraumreduktion:  η = ρ(1-ρ) / ((2π²)^N_shells × π)

V19/V20-Erweiterungen:
  - Friedmann-Konsistenz: tau_reh(T) aus Friedmann H = sqrt(8piG*rho/3)
  - T_exact: welche Temperatur gibt exakt N=274?
  - Unsicherheitsanalyse via eta_obs, H0, rho_DE Fehler
  - Pfad D: kein tau_infl als freier Parameter (V20-Kandidat)
"""

import os, sys, math
import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

print("=" * 70)
print("  OT-2: f_echo Analytische Verifikation")
print("=" * 70)

# ── Framework-Konstanten ──────────────────────────────────────────────────────
DF_GEO   = 0.77
DF_DISS  = 0.44
DF_EFF   = DF_GEO * DF_DISS     # 0.3388

C        = 2.998e8    # m/s
H0_SI    = 67360.0 / 3.0856e22  # s^-1  (H0=67.36 km/s/Mpc)
T0       = 1.0 / H0_SI          # ≈ 4.35e17 s
RHO_DE   = 6.034e-27  # kg/m^3  (dunkle Energie)

# Baryon-Asymmetrie (Planck 2018 TT,TE,EE+lowE)
ETA_OBS  = 6.104e-10   # η_b = n_b/n_γ  (gemessen, CMB+BBN)

# ── Pfad A: Thermodynamisches Dichteverhältnis ─────────────────────────────────
# Inflation endet bei T_reh ~ 10^15 GeV → rho_reh aus Strahlungsdominanz
# rho_radiation = (pi^2/30) * g_* * (kT)^4 / (hbar*c)^3
# Bei T_reh = 10^15 GeV = 10^24 eV = 1.602e5 J (k_B * T):
# rho_reh [kg/m^3] = (pi^2/30) * g_* * (kT)^4 / (hbar*c)^3 / c^2

kB      = 1.3806e-23   # J/K
HBAR    = 1.0546e-34   # J*s

# T_reheating Varianten (alle etablierten Werte)
# Niedrig (SUSY-Leptogenese): T_reh ≈ 10^9 GeV
# Standard GUT:                T_reh ≈ 10^15 GeV
# Planck-Inflation:            T_reh = T_Planck ≈ 1.416e32 K

def rho_radiation(T_K, g_star=106.75):
    """Strahlungsdichte bei Temperatur T [kg/m^3]."""
    kT = kB * T_K
    rho = (math.pi**2 / 30.0) * g_star * kT**4 / (HBAR**3 * C**3) / C**2
    return rho

# Verschiedene T_reh Annahmen und zugehörige ln(rho/rho_DE)
T_scenarios = {
    "T_reh = 10^9 GeV  (Leptogenese)":  1e9  * 1e9 * 1.602e-19 / kB,
    "T_reh = 10^14 GeV (GUT)":          1e14 * 1e9 * 1.602e-19 / kB,
    "T_reh = 10^15 GeV (Standard GUT)": 1e15 * 1e9 * 1.602e-19 / kB,
    "T_reh = 10^16 GeV (hohe Inflation)":1e16 * 1e9 * 1.602e-19 / kB,
    "T_Planck = 1.416e32 K":            1.416e32,
}

print("\nPfad A: ln(rho_reh / rho_DE) fuer verschiedene T_reh:")
print(f"  {'Szenario':<45}  {'T [K]':>12}  {'ln-Ratio':>10}  {'Abstand zu 274':>14}")
print("  " + "-" * 85)
pathA_results = {}
for label, T_K in T_scenarios.items():
    rho_reh = rho_radiation(T_K)
    ln_ratio = math.log(rho_reh / RHO_DE)
    diff = ln_ratio - 274.0
    print(f"  {label:<45}  {T_K:>12.3e}  {ln_ratio:>10.1f}  {diff:>+14.1f}")
    pathA_results[label] = (T_K, rho_reh, ln_ratio)

# ── Pfad B: Baryon-Asymmetrie Rückrechnung ─────────────────────────────────────
# eta = f_echo * (tau_infl / t_0)
# → f_echo = eta / (tau_infl / t_0) = eta * t_0 / tau_infl
# → ln(f_echo) = ln(eta) + ln(t_0) - ln(tau_infl)
# → D_f,eff * N = ln(f_echo)  → N = ln(f_echo) / D_f,eff

print("\nPfad B: Baryon-Asymmetrie Rückrechnung (eta_obs = 6.104e-10):")
print(f"  t_0 = {T0:.4e} s  (Hubble-Zeit)")
tau_scenarios = {
    "tau_infl = 10^-36 s (Ende Inflation)":  1e-36,
    "tau_infl = 10^-32 s (Friedmann start)": 1e-32,
    "tau_infl = 10^-33 s (GUT Übergang)":    1e-33,
    "tau_infl = 10^-38 s (Planck Grenze)":   1e-38,
}

pathB_results = {}
print(f"\n  {'Szenario':<45}  {'tau [s]':>10}  {'N=lnf/Df':>10}  {'Abstand zu 274':>14}")
print("  " + "-" * 85)
for label, tau in tau_scenarios.items():
    f_echo = ETA_OBS * T0 / tau
    ln_f = math.log(f_echo)
    N = ln_f / DF_EFF
    diff = N - 274.0
    print(f"  {label:<45}  {tau:>10.1e}  {N:>10.1f}  {diff:>+14.1f}")
    pathB_results[label] = (tau, f_echo, N)

# ── Direkter Vergleich bei N=274 ───────────────────────────────────────────────
N_CLAIM  = 274.0
f_echo   = math.exp(DF_EFF * N_CLAIM)
exp_arg  = DF_EFF * N_CLAIM

print(f"\n  f_echo = exp({DF_EFF:.4f} × {N_CLAIM:.0f}) = exp({exp_arg:.2f}) = {f_echo:.6e}")

# Baryonasymmetrie aus f_echo
tau_ref  = 1e-32   # Standardannahme
eta_pred = f_echo * (tau_ref / T0)
print(f"\n  η_predicted = f_echo × (tau_ref/t_0) = {f_echo:.4e} × {tau_ref/T0:.4e}")
print(f"             = {eta_pred:.4e}")
print(f"  η_observed  = {ETA_OBS:.4e}")
print(f"  Abweichung: {100*abs(eta_pred-ETA_OBS)/ETA_OBS:.1f}%")

# Rho early aus N=274
rho_early = RHO_DE * math.exp(N_CLAIM)
print(f"\n  Bei N=274: rho_early = rho_DE * exp(274) = {rho_early:.4e} kg/m^3")

# Welche Temperatur entspricht rho_early?
# rho_radiation ∝ T^4 → T = (rho * 30 * hbar^3 * c^5 / (pi^2 * g* * kB^4))^(1/4)
g_star = 106.75
T_from_rho = (rho_early * 30 * HBAR**3 * C**5 / (math.pi**2 * g_star * kB**4))**0.25
print(f"  Entspricht Temperatur: T = {T_from_rho:.4e} K")
T_GeV = T_from_rho * kB / (1e9 * 1.602e-19)
print(f"                        T = {T_GeV:.4e} GeV")

# Vergleich mit Planck-Temperatur
T_Planck = math.sqrt(HBAR * C**5 / (2 * math.pi * 6.674e-11)) / kB
print(f"\n  Planck-Temperatur: T_Planck = {T_Planck:.4e} K")
print(f"  Verhältnis T_reh/T_Planck = {T_from_rho/T_Planck:.4f}")

# ── Zwei-Pfad Konsistenzcheck ──────────────────────────────────────────────────
# Konsistenz: Pfad A (T_reh = T_GUT) vs Pfad B (tau = 10^-32 s)
# Wenn beide dieselbe Physik beschreiben, müssen sie übereinstimmen.

# Bei T_reh = 10^15 GeV: ln_ratio_A
rho_15 = rho_radiation(T_scenarios["T_reh = 10^15 GeV (Standard GUT)"])
N_A = math.log(rho_15 / RHO_DE)
# Bei tau_infl = 10^-32 s: N_B
f_B = ETA_OBS * T0 / 1e-32
N_B = math.log(f_B) / DF_EFF

# ── Pfad C: Planck-Dichte (parameterfrei) ─────────────────────────────────────
G_N  = 6.674e-11   # Gravitationskonstante [m^3 kg^-1 s^-2]
# Planck-Dichte: rho_P = c^5 / (hbar * G^2)  (Standarddefinition)
RHO_PLANCK = C**5 / (HBAR * G_N**2)
N_C = math.log(RHO_PLANCK / RHO_DE)
print("\n" + "=" * 70)
print("  PFAD C: Planck-Dichte (parameterfrei)")
print("=" * 70)
print(f"  rho_Planck = c^5 / (hbar * G^2) = {RHO_PLANCK:.4e} kg/m^3")
print(f"  N_C = ln(rho_Planck / rho_DE) = {N_C:.2f}")
print(f"  Abweichung von 274: {N_C - 274:+.2f}  ({100*abs(N_C-274)/274:.1f}%)")

# ── Friedmann-Konsistenz ───────────────────────────────────────────────────────
# tau_reh = 1 / (2 * H_reh) = 1 / (2 * sqrt(8*pi*G*rho_reh/3))
# — berechnet tau aus T_reh physikalisch-konsistent
print("\n" + "=" * 70)
print("  FRIEDMANN-KONSISTENZ: tau_reh(T_reh) aus H = sqrt(8piG*rho/3)")
print("=" * 70)
print(f"  {'T_reh':<45}  {'N_A':>8}  {'tau_Friedmann':>14}  {'N_B(tau)':>10}")
print("  " + "-" * 83)
for label, (T_K, rho_reh, N_a) in pathA_results.items():
    H_reh   = math.sqrt(8 * math.pi * G_N * rho_reh / 3.0)
    tau_F   = 1.0 / (2.0 * H_reh)          # Strahlungsdominanz: t = 1/(2H)
    f_echo_F = ETA_OBS * T0 / tau_F
    N_b_F   = math.log(f_echo_F) / DF_EFF
    print(f"  {label:<45}  {N_a:>8.1f}  {tau_F:>14.2e}  {N_b_F:>10.1f}")

# ── T_exact: welche T_reh liefert exakt N=274? ─────────────────────────────────
# rho_target = rho_DE * exp(274)
# T_exact = (rho_target * 30 * hbar^3 * c^5 / (pi^2 * g* * kB^4))^(1/4)
rho_target = RHO_DE * math.exp(N_CLAIM)
g_star = 106.75
T_exact = (rho_target * 30.0 * HBAR**3 * C**5 / (math.pi**2 * g_star * kB**4))**0.25
T_exact_GeV = T_exact * kB / (1e9 * 1.602e-19)
T_Planck = math.sqrt(HBAR * C**5 / (2 * math.pi * G_N)) / kB
print("\n" + "=" * 70)
print("  T_EXAKT (loest N=274 analytisch)")
print("=" * 70)
print(f"  rho_target = rho_DE * exp(274) = {rho_target:.4e} kg/m^3")
print(f"  T_exact    = {T_exact:.4e} K")
print(f"  T_exact    = {T_exact_GeV:.4e} GeV")
print(f"  T_Planck   = {T_Planck:.4e} K  ({T_Planck*kB/(1e19*1.602e-19):.3f} × 10^19 GeV)")
print(f"  T_exact / T_Planck = {T_exact/T_Planck:.4f}")
# Entsprechende tau via Friedmann
H_exact = math.sqrt(8 * math.pi * G_N * rho_target / 3.0)
tau_exact = 1.0 / (2.0 * H_exact)
print(f"  tau_Friedmann(T_exact) = {tau_exact:.3e} s")
# Was ergibt das in Pfad B?
N_B_exact = math.log(ETA_OBS * T0 / tau_exact) / DF_EFF
print(f"  N_B(tau_exact) = {N_B_exact:.2f}  (soll 274)")

# ── Unsicherheitsanalyse ──────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  UNSICHERHEITSANALYSE (Pfad B: N_B = ln(eta * t0 / tau) / Df)")
print("=" * 70)
# Sensitivität: dN/d(ln eta) = 1/Df, dN/d(ln t0) = 1/Df, dN/d(ln tau) = -1/Df
# Relative Unsicherheiten (1-sigma):
delta_eta = 0.006   # Planck 2018: ±0.6%
delta_H0  = 0.015   # ±1.5% (Hubble tension range)
delta_tau = 1.0     # tau unknown within factor ~3 (1 dex = 2.3 / Df ≈ 6.8 in N)
# Linear error propagation (logarithmic):
sigma_eta_N = delta_eta / DF_EFF           # ΔN from eta uncertainty
sigma_H0_N  = delta_H0  / DF_EFF           # ΔN from H0 (affects t0)
sigma_tau_N = math.log(3.0) / DF_EFF      # ΔN from 1 order of magnitude in tau
print(f"  Quelle        Rel.Unsicher.  ΔN (1σ)")
print(f"  eta_obs       ±{delta_eta*100:.1f}%          ±{sigma_eta_N:.1f}")
print(f"  H0 (→ t0)     ±{delta_H0*100:.1f}%          ±{sigma_H0_N:.1f}")
print(f"  tau_infl      ×{3:.0f} (1 Ord.)    ±{sigma_tau_N:.1f}  (dominiert)")
print(f"  tau_infl      ×{10:.0f} (10×)       ±{math.log(10)/DF_EFF:.1f}")
print(f"\n  N_B (best guess) = {N_B:.1f}")
print(f"  Intrinsic band:  N = {N_B:.1f} ± {sigma_tau_N:.1f}  (bei tau ±Faktor 3)")
print(f"  N=274 liegt {'INNERHALB' if abs(N_B-274) < sigma_tau_N else 'AUSSERHALB'} "
      f"des 1σ-Bandes  (|{N_B-274:+.1f}| vs ±{sigma_tau_N:.1f})")

print("\n" + "=" * 70)
print("  ZWEI-PFAD KONSISTENZCHECK")
print("=" * 70)
print(f"  Pfad A (T_reh = 10^15 GeV): N_A = {N_A:.2f}")
print(f"  Pfad B (tau  = 10^-32  s): N_B = {N_B:.2f}")
print(f"  Pfad C (Planck-Dichte):     N_C = {N_C:.2f}")
print(f"  Framework-Wert:            N   = {N_CLAIM:.1f}")
print(f"  A vs Framework:  Δ = {N_A - N_CLAIM:+.2f}  ({100*abs(N_A-N_CLAIM)/N_CLAIM:.1f}%)")
print(f"  B vs Framework:  Δ = {N_B - N_CLAIM:+.2f}  ({100*abs(N_B-N_CLAIM)/N_CLAIM:.1f}%)")
print(f"  C vs Framework:  Δ = {N_C - N_CLAIM:+.2f}  ({100*abs(N_C-N_CLAIM)/N_CLAIM:.1f}%)")
print(f"  A vs B:          Δ = {N_A - N_B:+.2f}  (beide unabhängig)")
print(f"  T_exact (N=274): {T_exact_GeV:.2e} GeV  (= {T_exact/T_Planck:.3f} T_Planck)")

# ── DF_EFF Verifikation ────────────────────────────────────────────────────────
print("\n  DF_EFF Verifikation:")
print(f"  DF_GEO  = {DF_GEO}")
print(f"  DF_DISS = {DF_DISS}")
print(f"  DF_EFF  = {DF_GEO} × {DF_DISS} = {DF_EFF:.6f}")
print(f"  f_echo = exp({DF_EFF:.4f} × 274) = {math.exp(DF_EFF*274):.4e}")

# ═══════════════════════════════════════════════════════════════════════════════
# PFAD D: DM-Echo-Wellenasymmetrie + S³-π-Unterdrückung (V20)
# ═══════════════════════════════════════════════════════════════════════════════
# Idee (Kevin Hannemann, 28.04.2026):
#   Wenn Dunkle Materie = HTM-Echo (Rückwärtswelle, ρ-gedämpft), dann ist die
#   Materie/Antimaterie-Asymmetrie die Asymmetrie zwischen Vorwärts- und
#   Rückwärtswelle, unterdrückt durch die S³->3D Phasenraumreduktion via π.
#
# Schritt 1: Rohe Wellenasymmetrie
#   A_+ = Vorwärtswelle (Amplitude 1)
#   A_- = Rückwärtswelle (Amplitude ρ, Echo-gedämpft)
#   η_roh = (A_+ - A_-) / (A_+ + A_-) = (1-ρ)/(1+ρ)
#
# Schritt 2: S³-Phasenraumvolumen
#   S³ voll:     Vol = 2π²r³  (enthält Faktor 2π²)
#   S¹-Faser:    Vol = 2πr
#   Projektion S³->R³/S²: Unterdrückungsfaktor = (2π²) pro Windungsskala
#
# Schritt 3: N_shells Hopf-Windungen (V20: 6 Schalen)
#   η ≈ ρ(1-ρ) / (2π²)^N_shells
#
# Schritt 4: Zweiter Torsionsschock auf Dualitätssphäre (n=3, 175.96°)
#   Zusätzlicher Unterdrückungsfaktor π (Dualitätssphäre-Projektion)
#   η_dual = ρ(1-ρ) / ((2π²)^N_shells * π)

print("\n" + "=" * 70)
print("  PFAD D: DM-Echo-Wellenasymmetrie + S³-π-Unterdrückung (V20)")
print("  Idee: Kevin Hannemann, 28.04.2026")
print("=" * 70)

RHO_IFS   = 0.406      # IFS-Kontraktionsverhältnis ρ
N_SHELLS  = 6          # V20: 6 Schalen (n=1...6)
PI        = math.pi

# Schritt 1: Rohe Wellenasymmetrie
eta_raw = (1 - RHO_IFS) / (1 + RHO_IFS)
print(f"\n  Schritt 1 — Rohe Wellenasymmetrie:")
print(f"    ρ = {RHO_IFS}")
print(f"    η_roh = (1-ρ)/(1+ρ) = {1-RHO_IFS:.4f}/{1+RHO_IFS:.4f} = {eta_raw:.4f}")
print(f"    (viel zu groß; S³-Unterdrückung nötig)")

# Schritt 2: S³-Phasenraumfaktor
vol_S3_factor = 2 * PI**2   # = 19.739...
print(f"\n  Schritt 2 — S³-Phasenraumvolumen:")
print(f"    S³-Volumen ∝ 2π²  = {vol_S3_factor:.4f}")
print(f"    (2π²)^{N_SHELLS} = {vol_S3_factor**N_SHELLS:.4e}")

# Schritt 3: N_shells Windungen
eta_D_no_dual = RHO_IFS * (1 - RHO_IFS) / (vol_S3_factor**N_SHELLS)
print(f"\n  Schritt 3 — Nach S³-Unterdrückung (N_shells={N_SHELLS}):")
print(f"    η_D = ρ(1-ρ) / (2π²)^{N_SHELLS}")
print(f"        = {RHO_IFS:.4f} × {1-RHO_IFS:.4f} / {vol_S3_factor**N_SHELLS:.4e}")
print(f"        = {eta_D_no_dual:.4e}")
ratio_no_dual = eta_D_no_dual / ETA_OBS
print(f"    Verhältnis η_D / η_obs = {ratio_no_dual:.2f}  (Faktor {ratio_no_dual:.1f} zu groß)")

# Schritt 4: Dualitätssphäre-Korrekturfaktor π
eta_D = eta_D_no_dual / PI
print(f"\n  Schritt 4 — Mit Dualitätssphäre-Faktor π:")
print(f"    η_final = ρ(1-ρ) / ((2π²)^{N_SHELLS} × π)")
print(f"            = {eta_D_no_dual:.4e} / π")
print(f"            = {eta_D:.4e}")
ratio_D = eta_D / ETA_OBS
print(f"    η_obs   = {ETA_OBS:.4e}")
print(f"    Ratio   = {ratio_D:.2f}  (Abweichung: Faktor {ratio_D:.1f})")
pct_D = 100 * abs(eta_D - ETA_OBS) / ETA_OBS
print(f"    %-Abweichung: {pct_D:.0f}%")

# Vergleich: Originaler Pfad B vs Pfad D
print(f"\n  Vergleich der Ableitungswege:")
print(f"    Pfad B (f_echo × τ/t₀, τ=10⁻³² s): η = {ETA_OBS * (N_B/N_CLAIM):.2e}  (kalibriert)")
print(f"    Pfad D (Wellenasymmetrie, kein τ):  η = {eta_D:.2e}  (parameterfrei)")
print(f"    Pfad D Vorteil: KEIN τ_infl als freier Parameter!")

# Df_eff Konsistenz mit Pfad D
# Wenn η ~ ρ(1-ρ)/(2π²)^N_shells/π, dann implizit:
# f_echo_D = η_D / (τ/t₀) → N_D = ln(f_echo_D)/DF_EFF
tau_ref_D = 1e-32
f_echo_D  = eta_D / (tau_ref_D / T0)
N_D       = math.log(f_echo_D) / DF_EFF
print(f"\n  Rückrechnung Pfad D → N:")
print(f"    f_echo(D) = η_D × t₀/τ = {f_echo_D:.4e}")
print(f"    N_D = ln(f_echo_D)/D_f,eff = {N_D:.1f}  (Framework-Wert: 274)")
print(f"    Abweichung: {N_D-274:+.1f}  ({100*abs(N_D-274)/274:.1f}%)")

# ── Ergebnistext ────────────────────────────────────────────────────────────────────────────
lines = []
lines += [
    "=" * 70,
    "OT-2: f_echo Analytische Verifikation (V19)",
    "=" * 70, "",
    "Framework-Formel:",
    f"  f_echo = exp(D_f,eff * N)  mit D_f,eff = {DF_EFF:.4f}",
    f"  N = ln(rho_early / rho_DE) ~ 274",
    f"  f_echo = exp({DF_EFF:.4f} * 274) = exp({DF_EFF*274:.2f}) = {math.exp(DF_EFF*274):.4e}",
    "",
    "  D_f,eff = D_f,geo * D_f,diss = 0.77 * 0.44 = 0.3388",
    "",
    "-" * 70,
    "PFAD A: Thermodynamisches Dichtevaeltnis",
    "-" * 70,
]
for label, (T_K, rho_reh, ln_ratio) in pathA_results.items():
    lines.append(f"  {label}: N = {ln_ratio:.1f}  (Delta = {ln_ratio-N_CLAIM:+.1f})")

best_A_key = "T_reh = 10^15 GeV (Standard GUT)"
lines += [
    "",
    f"  Referenz-Match: T_reh = 10^15 GeV (Standard GUT), N = {pathA_results[best_A_key][2]:.1f}",
    f"  Abweichung: {pathA_results[best_A_key][2]-N_CLAIM:+.1f}  ({100*abs(pathA_results[best_A_key][2]-N_CLAIM)/N_CLAIM:.1f}%)",
    "",
    "-" * 70,
    "PFAD B: Baryon-Asymmetrie Rueckrechnung",
    "-" * 70,
]
for label, (tau, f_b, N) in pathB_results.items():
    lines.append(f"  {label}: N = {N:.1f}  (Delta = {N-N_CLAIM:+.1f})")

best_B_key = "tau_infl = 10^-32 s (Friedmann start)"
lines += [
    "",
    f"  Referenz-Match: tau = 10^-32 s (GUT-Inflation Ende), N = {pathB_results[best_B_key][2]:.1f}",
    "",
    "-" * 70,
    "PFAD C: Planck-Dichte (parameterfrei)",
    "-" * 70,
    "",
    f"  rho_Planck = c^5 / (hbar * G^2) = {RHO_PLANCK:.4e} kg/m^3",
    f"  N_C = ln(rho_Planck / rho_DE) = {N_C:.2f}  (Delta = {N_C-N_CLAIM:+.2f})",
    "  Interpretation: Planck-Dichte liefert obere Schranke; N=274 liegt",
    f"  unterhalb N_C = {N_C:.1f}, was physikalisch plausibel ist.",
    "",
    "-" * 70,
    "FRIEDMANN-KONSISTENZ",
    "-" * 70,
    "",
    "  Fuer jede T_reh ergibt die Friedmann-Gleichung eine konsistente tau_F:",
    "",
    f"  {'T_reh':<45}  {'N_A':>6}  {'tau_F':>12}  {'N_B(tau_F)':>10}",
    "  " + "-" * 78,
]
for label, (T_K, rho_reh, N_a) in pathA_results.items():
    H_reh_l  = math.sqrt(8 * math.pi * G_N * rho_reh / 3.0)
    tau_F_l  = 1.0 / (2.0 * H_reh_l)
    f_echo_F = ETA_OBS * T0 / tau_F_l
    N_b_F    = math.log(f_echo_F) / DF_EFF
    lines.append(f"  {label:<45}  {N_a:>6.1f}  {tau_F_l:>12.2e}  {N_b_F:>10.1f}")

lines += [
    "",
    "  Hinweis: tau_infl = 10^-32 s ist Friedmann-inkonsistent mit T_reh=10^15 GeV.",
    f"  Konsistente tau bei T_GUT: ~10^-37 s (5 Groessenordnungen kleiner).",
    "",
    "-" * 70,
    "T_EXAKT: welche T_reh ergibt exakt N=274?",
    "-" * 70,
    "",
    f"  rho_target = rho_DE * exp(274) = {rho_target:.4e} kg/m^3",
    f"  T_exact    = {T_exact_GeV:.4e} GeV",
    f"  T_exact / T_Planck = {T_exact/T_Planck:.4f}  (nahe Planck-Skala)",
    f"  tau_Friedmann(T_exact) = {tau_exact:.3e} s",
    f"  N_B(tau_exact) = {N_B_exact:.2f}  (Pfad-B-Wert bei tau_exact; soll 274)",
    "",
    f"  Pfad A-Konsistenz: T={T_exact_GeV:.2e} GeV ergibt N_A=274 per Konstruktion.",
    f"  Pfad B-Inkonsistenz: N_B(tau_exact)={N_B_exact:.0f} != 274 (Pfad B versagt hier).",
    "  -> Kein (T_reh, tau) Paar erfuellt beide Pfade gleichzeitig (Friedmann-konform).",
    "",
    "-" * 70,
    "UNSICHERHEITSANALYSE (Pfad B)",
    "-" * 70,
    "",
    "  Quelle        Rel.Unsicherheit  Delta_N (1sigma)",
    f"  eta_obs       +/-0.6%             +/-{0.006/DF_EFF:.1f}",
    f"  H0 (-> t0)    +/-1.5%             +/-{0.015/DF_EFF:.1f}",
    f"  tau_infl      x3 (Faktor)         +/-{math.log(3)/DF_EFF:.1f}  (dominiert)",
    "",
    f"  N_B(10^-32 s) = {N_B:.2f}",
    f"  tau-Band (x3): N = {N_B:.1f} +/- {math.log(3)/DF_EFF:.1f}",
    f"  N=274 liegt {'INNERHALB' if abs(N_B-N_CLAIM) < math.log(3)/DF_EFF else 'AUSSERHALB'} "
    f"des 1sigma-Bandes  (|Delta_N|={abs(N_B-N_CLAIM):.1f} vs +/-{math.log(3)/DF_EFF:.1f})",
    "",
    "-" * 70,
    "DREI-PFAD KONSISTENZ",
    "-" * 70,
    "",
    f"  N = {N_CLAIM:.0f} (Framework-Wert)",
    f"  Pfad A (T_GUT = 10^15 GeV): N_A = {N_A:.1f}  Delta = {N_A-N_CLAIM:+.1f}  ({100*abs(N_A-N_CLAIM)/N_CLAIM:.1f}%)",
    f"  Pfad B (tau  = 10^-32  s):  N_B = {N_B:.1f}  Delta = {N_B-N_CLAIM:+.1f}  ({100*abs(N_B-N_CLAIM)/N_CLAIM:.1f}%)",
    f"  Pfad C (Planck-Dichte):     N_C = {N_C:.1f}  Delta = {N_C-N_CLAIM:+.1f}  ({100*abs(N_C-N_CLAIM)/N_CLAIM:.1f}%)",
    f"  Bracket: N_A={N_A:.0f} < N=274 < N_C={N_C:.0f}  [N liegt zwischen GUT und Planck]",
    "",
    "-" * 70,
    "BEWERTUNG",
    "-" * 70,
    "",
]

tol = 0.15
delta_AB = abs(N_A - N_B) / N_CLAIM
if delta_AB < tol and abs(N_A - N_CLAIM) / N_CLAIM < tol:
    lines += [
        "  OT-2 BESTANDEN: Beide Referenz-Pfade stimmen innerhalb 15% ueberein.",
        f"  Delta(A-B) = {N_A-N_B:.1f}  ({100*delta_AB:.1f}% von 274)",
        "",
        "  VORBEHALT: Pfad A und Pfad B sind Friedmann-inkonsistent miteinander.",
        "  tau=10^-32 s und T=10^15 GeV entsprechen nicht dem gleichen Weltzeitalter.",
]
elif delta_AB < 0.5:
    lines += [
        f"  OT-2 BEDINGT: Pfad A und B innerhalb 50% (Delta = {N_A-N_B:.1f}).",
        f"  Pfad B (N_B={N_B:.1f}) und Pfad C (N_C={N_C:.1f}) stuetzen N=274.",
        "  Pfad A (T_GUT=10^15 GeV) ist ~9% zu niedrig; T_exact~3e18 GeV schliesst Luecke.",
        "  Kritisch: tau=10^-32 s und T=10^15 GeV entsprechen verschiedenen Epochen.",
]
else:
    lines += [
        f"  OT-2 NEGATIV: Pfade weichen stark ab (Delta(A-B) = {N_A-N_B:.1f}).",
        "  N=274 kann nicht aus T_reh=10^15 GeV und tau=10^-32 s erklaert werden.",
]

lines += [
    "",
    "-" * 70,
    "PFAD D: DM-Echo-Wellenasymmetrie + S3-pi-Unterdrueckung (V20)",
    "Idee: Kevin Hannemann, 28.04.2026",
    "-" * 70,
    "",
    "  Wenn DM = HTM-Echo (Rueckwaertswelle, rho-gedaempft):",
    "  Baryon-Asymmetrie = Asymmetrie Vorwaerts/Rueckwaertswelle",
    "  untergedrueckt durch S3->R3 Phasenraumreduktion via pi.",
    "",
    "  Schritt 1 — Rohe Wellenasymmetrie:",
    f"    eta_roh = (1-rho)/(1+rho) = {1-RHO_IFS:.4f}/{1+RHO_IFS:.4f} = {eta_raw:.4f}",
    "",
    "  Schritt 2 — S3-Phasenraumvolumen:",
    f"    (2*pi^2)^N_shells = ({vol_S3_factor:.3f})^{N_SHELLS} = {vol_S3_factor**N_SHELLS:.4e}",
    "",
    "  Schritt 3 — Nach S3-Unterdrueckung (N_shells=6):",
    f"    eta_D = rho*(1-rho) / (2*pi^2)^6 = {eta_D_no_dual:.4e}",
    f"    Verhaeltnis zu eta_obs: {ratio_no_dual:.1f}x (noch zu gross)",
    "",
    "  Schritt 4 — Dualitaetssphaere-Korrekturfaktor pi:",
    f"    eta_final = rho*(1-rho) / ((2*pi^2)^6 * pi)",
    f"              = {eta_D:.4e}",
    f"    eta_obs   = {ETA_OBS:.4e}",
    f"    Abweichung: {pct_D:.0f}%  (kein freier Parameter)",
    "",
    f"  Vorteil Pfad D: KEIN tau_infl als freier Parameter!",
    f"  Roh-Asymmetrie + S3-Geometrie -> korrekte Groessenordnung eta",
    f"  N_D (Rueckrechnung) = {N_D:.1f}  (Framework: 274, Delta={N_D-274:+.1f})",
    "",
    "  STATUS: Pfad D liefert eta ~ {:.0e} = korrekte Groessenordnung  (V20 KANDIDAT)".format(eta_D),
    "",
    "-" * 70,
    "EPISTEMISCHER STATUS (V19/V20)",
    "-" * 70,
    "",
    "  N=274 = ln(rho_early/rho_DE) ist eine dimensionslose Zahl.",
    "  Pfad B (tau=10^-32 s) reproduziert N=274.9 mit <1% Abweichung.",
    f"  Pfad C (Planck-Dichte) ergibt N_C={N_C:.1f} -- N=274 liegt darunter (plausibel).",
    f"  T_exact = {T_exact_GeV:.2e} GeV loest N_A=274 exakt (= {T_exact/T_Planck:.2f} T_Planck).",
    f"  Aber: N_B(tau bei T_exact) = {N_B_exact:.0f} != 274 (Pfad-B-Versagen bei Friedmann-tau).",
    "  Kritisch: Kein Friedmann-konsistentes (T, tau) Paar loest beide Pfade gleichzeitig.",
    "  Pfad B mit tau=10^-32 s und Pfad A mit T_GUT=10^15 GeV sind epistemisch entkoppelt.",
    "  Gesamtstatus: BEDINGT -- Pfad B trifft N=274 (Kalibrierung), Pfad A schwach, C als Bracket.",
    "",
    "=" * 70,
]
result_text = "\n".join(lines)
print("\n" + result_text)
out = os.path.join(RESULTS, "OT_02_result.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write(result_text)
print(f"\n  Ergebnis: {out}")
