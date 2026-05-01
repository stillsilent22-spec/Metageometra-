"""
OT-2: f_echo aus ersten Prinzipien
===================================
Zwei unabhaengige Pfade zu f_echo = exp(D_eff * N)
"""
import math, os
import numpy as np

BASE   = os.path.dirname(os.path.abspath(__file__))
RESDIR = os.path.join(BASE, "results")
os.makedirs(RESDIR, exist_ok=True)

# ── Parameter ────────────────────────────────────────────
RHO_0    = 0.406
D_GEO    = math.log(2) / math.log(1/RHO_0)          # 0.769
D_DISS   = 1.0 / (3.0 * D_GEO)                       # 0.433
D_EFF    = D_GEO * D_DISS                             # ~1/3
RHO_DE   = 6.034e-27  # kg/m3
T0       = 4.352e17   # s
C        = 2.998e8
HBAR     = 1.0546e-34
G        = 6.674e-11
CHI      = 59.1
ETA_OBS  = 5.58e-10

# ── Pfad 1: Direkte Berechnung ────────────────────────────
RHO_PLANCK = C**5 / (HBAR * G**2)
N_ratio    = math.log(RHO_PLANCK / RHO_DE)           # dimensionslos
N_log10    = math.log10(RHO_PLANCK / RHO_DE)
f_echo_1   = math.exp(D_EFF * N_ratio)               # Hauptpfad
f_echo_ref = 2.07e40                                   # aus V20 OT-2

# ── Pfad 2: Cavity Eigenmode Ableitung ───────────────────
# Cavity-Grundmode: omega_1 = a0/c, r_s = c^2/(2pi*a0)
# Nullpunktsenergie der Cavity: E_0 = hbar * omega_1 / 2
# Anzahl Moden N_modes: N_modes = r_Planck / r_s (geometrisch)
a0        = C / (2 * math.pi * T0)
r_s       = C**2 / (2 * math.pi * a0)                # = c * T0 / (2pi)^2
r_Planck  = math.sqrt(HBAR * G / C**3)
N_modes   = r_s / r_Planck                            # Anzahl Planck-Laengen in r_s

# Nullpunkts-Amplifikation: f_echo ~ N_modes^D_EFF
f_echo_2  = N_modes**D_EFF

# ── Pfad 3: eta-Zeitskala Verknuepfung ───────────────────
# eta = (1-cos(chi)) / ((2-cos(chi)) * (2pi^2)^(6+D_GEO))
eta_theory = (1-math.cos(math.radians(CHI))) / ((2-math.cos(math.radians(CHI))) * (2*math.pi**2)**(6+D_GEO))
# Vorhersage: f_echo * eta_timescale = eta_obs
# eta_timescale = T_Planck / T0  (schmalste Zeitskala der Cavity)
T_PLANCK = math.sqrt(HBAR * G / C**5)
eta_timescale = T_PLANCK / T0
f_echo_eta    = ETA_OBS / eta_timescale

# ── Ausgabe ──────────────────────────────────────────────
print("══════════════════════════════════════════════════════")
print("OT-2: f_echo aus ersten Prinzipien (V21.3)")
print("══════════════════════════════════════════════════════")
print(f"  D_f,geo   = {D_GEO:.6f}")
print(f"  D_f,diss  = {D_DISS:.6f}")
print(f"  D_f,eff   = {D_EFF:.6f}  (~1/3 = {1/3:.6f})")
print()
print(f"  rho_Planck = {RHO_PLANCK:.4e} kg/m3")
print(f"  rho_DE     = {RHO_DE:.4e} kg/m3")
print(f"  N = ln(rho_P/rho_DE) = {N_ratio:.4f}  [= {N_log10:.2f} * ln(10)]")
print()
print(f"  ── Pfad 1: exp(D_eff * N) ──")
print(f"  f_echo_1  = exp({D_EFF:.4f} * {N_ratio:.2f}) = {f_echo_1:.4e}")
print(f"  f_echo_ref= {f_echo_ref:.2e}  (Beobachtet V20)")
print(f"  Abweichung: {abs(math.log10(f_echo_1)-math.log10(f_echo_ref)):.3f} dex")
print()
print(f"  ── Pfad 2: Cavity Eigenmode r_s/r_Planck ──")
print(f"  r_s (Cavity-Radius) = {r_s:.4e} m = {r_s/3.086e22:.1f} Mpc")
print(f"  r_Planck            = {r_Planck:.4e} m")
print(f"  N_modes             = r_s/r_P = {N_modes:.4e}")
print(f"  f_echo_2  = N_modes^D_eff = {f_echo_2:.4e}")
print(f"  Abweichung zu Pfad 1: {abs(math.log10(f_echo_1)-math.log10(f_echo_2)):.3f} dex")
print()
print(f"  ── Pfad 3: eta * f_echo = eta_timescale ──")
print(f"  eta_theory  = {eta_theory:.4e}")
print(f"  eta_obs     = {ETA_OBS:.2e}")
print(f"  T_Planck    = {T_PLANCK:.4e} s")
print(f"  eta_timescale = T_P/T0 = {eta_timescale:.4e}")
print(f"  f_echo_eta  = eta_obs / eta_timescale = {f_echo_eta:.4e}")
print(f"  Kohärenz Pfad1/Pfad3: {abs(math.log10(f_echo_1)-math.log10(f_echo_eta)):.2f} dex")
print()
p1p2 = abs(math.log10(f_echo_1)-math.log10(f_echo_2))
status = "BESTAETIGT" if p1p2 < 2.0 else "INCONCLUSIVE"
print(f"  Status: {status}")
print(f"  Pfad 1 & 2 Konsistenz: {p1p2:.2f} dex  ({'OK' if p1p2<2 else 'DIVERGIERT'})")
print("══════════════════════════════════════════════════════")

lines = [
    "════════════════════════════════════════════════════════════",
    "OT-2: f_echo aus ersten Prinzipien (V21.3)",
    "════════════════════════════════════════════════════════════",
    f"D_f,eff  = {D_EFF:.6f} (~1/3)",
    f"N = ln(rho_Planck/rho_DE) = {N_ratio:.4f}",
    f"",
    f"Pfad 1: exp(D_eff*N)         = {f_echo_1:.4e}",
    f"Pfad 2: (r_s/r_Planck)^D_eff = {f_echo_2:.4e}",
    f"Pfad 3: eta_obs/eta_ts        = {f_echo_eta:.4e}",
    f"Referenz (obs):               = {f_echo_ref:.2e}",
    f"",
    f"Pfad1-Pfad2 Abweichung: {p1p2:.2f} dex",
    f"Status: {status}",
    "Naechster OT: OT-7 (w(z) DESI)",
    "════════════════════════════════════════════════════════════",
]
with open(os.path.join(RESDIR, 'OT_2_result.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines)+'\n')
print(f"Gespeichert: {os.path.join(RESDIR,'OT_2_result.txt')}")
