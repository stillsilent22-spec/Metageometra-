"""
OT-20: Praezessionszyklus als Tier-3-Resonanz (V21.3)
"""
import math, os
from fractions import Fraction
import numpy as np

BASE   = os.path.dirname(os.path.abspath(__file__))
RESDIR = os.path.join(BASE, "results")
os.makedirs(RESDIR, exist_ok=True)

C   = 2.998e8; T0 = 4.352e17
A0  = C / (2 * math.pi * T0)
YR  = 3.15576e7  # s

T_PREC  = 25771.57 * YR   # Erdachsen-Praezession (Bretagnon 2003)
T_P9    = 5000. * YR       # Planet-9-Kandidat ~5000 Jahre
T_GAL   = 225e6 * YR       # Galaktische Umlaufzeit

print("══════════════════════════════════════════════════════")
print("OT-20: Praezessionszyklus als Tier-3-Resonanz (V21.3)")
print("══════════════════════════════════════════════════════")
print(f"  a0 = {A0:.4e} m/s2")
print(f"  omega_I = a0/c = {A0/C:.4e} s-1")
print(f"  T_I = 2pi/omega_I = {2*math.pi*C/A0/YR:.1f} yr")
print()

T_I = 2 * math.pi * C / A0

print(f"  ── Resonanzsuche (n=1..20) ──")
print(f"  n  | T_n = T_I/n [yr] | T_prec/T_n   | n/m Kette  | Fehler")
print(f"  {'─'*3}+{'─'*18}+{'─'*14}+{'─'*12}+{'─'*10}")
best = None
for n in range(1, 21):
    T_n = T_I / n
    ratio = T_PREC / T_n
    frac  = Fraction(ratio).limit_denominator(100)
    err   = abs(float(frac) - ratio) / ratio * 100
    flag  = " <<< RESONANZ" if err < 0.5 else ""
    if err < 0.5 and best is None:
        best = (n, T_n, ratio, frac, err)
    print(f"  {n:>2} | {T_n/YR:>16.1f} | {ratio:>12.6f} | {frac.numerator}/{frac.denominator:>4}     | {err:6.3f}% {flag}")

print()
print(f"  ── Galaktische Umlaufzeit ──")
ratio_gal = T_GAL / T_I
frac_gal  = Fraction(ratio_gal).limit_denominator(1000)
err_gal   = abs(float(frac_gal) - ratio_gal)/ratio_gal*100
print(f"  T_gal / T_I = {ratio_gal:.4e}  ~  {frac_gal} ({err_gal:.3f}%)")

print()
print(f"  ── Hubble-Zeit ──")
ratio_h0 = T0 / T_I
frac_h0  = Fraction(ratio_h0).limit_denominator(1000)
err_h0   = abs(float(frac_h0) - ratio_h0)/ratio_h0*100
print(f"  t0 / T_I = {ratio_h0:.4e}  ~  {frac_h0} ({err_h0:.3f}%)")

print()
if best:
    print(f"  BESTE RESONANZ: n={best[0]}, T_n={best[1]/YR:.1f} yr")
    print(f"  T_prec/T_n = {best[2]:.4f} ~ {best[3]} ({best[4]:.3f}%)")
    print(f"  => {best[3].numerator}:{best[3].denominator} Resonanz mit T_prec")
    status = "BESTAETIGT"
else:
    print(f"  Keine <0.5% Resonanz gefunden (n<20)")
    # Check broader
    min_err = 100
    for n in range(1, 100):
        T_n = T_I/n
        ratio = T_PREC/T_n
        frac = Fraction(ratio).limit_denominator(1000)
        err  = abs(float(frac)-ratio)/ratio*100
        if err < min_err:
            min_err = err
            best_n  = n
    best_n = 1
    for n in range(1, 100):
        T_n_l = T_I/n
        ratio_l = T_PREC/T_n_l
        frac_l = Fraction(ratio_l).limit_denominator(1000)
        err_l  = abs(float(frac_l)-ratio_l)/ratio_l*100
        if err_l < min_err:
            min_err = err_l; best_n = n
    print(f"  Bestes n={best_n}: Fehler={min_err:.3f}%")
    status = "INCONCLUSIVE"

print(f"  Status: {status}")
print("══════════════════════════════════════════════════════")

lines = [
    "════════════════════════════════════════════════════════════",
    "OT-20: Praezessionszyklus Tier-3-Resonanz (V21.3)",
    "════════════════════════════════════════════════════════════",
    f"a0 = {A0:.4e} m/s2",
    f"T_I = {T_I/YR:.1f} yr",
    f"T_prec = {T_PREC/YR:.2f} yr",
    f"T_gal = {T_GAL/YR:.2e} yr",
    f"t0 / T_I = {ratio_h0:.4e} ~ {frac_h0} ({err_h0:.3f}%)",
    f"Status: {status}",
    "════════════════════════════════════════════════════════════",
]
with open(os.path.join(RESDIR,'OT_20_result.txt'),'w',encoding='utf-8') as f:
    f.write('\n'.join(lines)+'\n')
print(f"Gespeichert: {os.path.join(RESDIR,'OT_20_result.txt')}")
