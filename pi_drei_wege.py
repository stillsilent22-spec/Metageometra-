"""Pi berechnen auf vier unabhaengigen Wegen (inkl. Hausdorff-Box-Counting)."""
import math
import numpy as np

# ── Pfad A: Leibniz-Reihe (vektorisiert) ──────────────────────────────
# pi/4 = 1 - 1/3 + 1/5 - 1/7 + ...
def leibniz(n_terms=10_000_000):
    k = np.arange(n_terms, dtype=np.float64)
    return 4.0 * np.sum((-1.0)**k / (2*k + 1))

# ── Pfad B: Monte-Carlo (vektorisiert) ────────────────────────────────
# Punkte im Einheitsquadrat; Anteil im Einheitskreis -> pi/4
def monte_carlo(n_points=5_000_000, seed=42):
    rng = np.random.default_rng(seed)
    xy  = rng.random((n_points, 2))
    inside = np.sum(xy[:,0]**2 + xy[:,1]**2 <= 1.0)
    return 4 * inside / n_points

# ── Pfad C: Nilakantha-Reihe (vektorisiert) ───────────────────────────
# pi = 3 + sum_{k=1}^{N} (-1)^(k+1) * 4 / ((2k)(2k+1)(2k+2))
def nilakantha(n_terms=500_000):
    k     = np.arange(1, n_terms + 1, dtype=np.float64)
    denom = (2*k) * (2*k + 1) * (2*k + 2)
    return 3.0 + np.sum((-1.0)**(k+1) * 4.0 / denom)

# ── Pfad D: Archimedes / Hausdorff-H¹-Approximation ───────────────────
# H¹(Einheitskreis) = 2*pi.
# Einbeschriebenes n-Eck: Umfang = 2n*sin(pi/n) -> 2*pi.
# Halbwinkelrekursion (KEINE Trig-Funktion), numerisch STABIL:
#   Zustandsvariable c = 1 - cos(theta)  (Differenz statt Kosinuswert)
#   sin(theta/2)     = sqrt(c / 2)
#   1-cos(theta/2)   = (c/2) / (1 + sqrt(1 - c/2))   <- kein Ausloeschung
# Start: n=6, c=1-cos(pi/6)=1-sqrt(3)/2
# pi = n * sin(pi/n) = n * sin_half  (sin des aktuellen Halbwinkels)
def archimedes_hausdorff(n_doublings=46):
    n        = 6
    c        = 1.0 - math.sqrt(3) / 2          # 1 - cos(pi/6)
    sin_half = 0.5                              # sin(pi/6), Startwert
    for _ in range(n_doublings):
        sin_half = math.sqrt(c / 2)             # sin(pi/(2n))  [aktueller Halbwinkel]
        c        = (c / 2) / (1.0 + math.sqrt(1.0 - c / 2))  # 1-cos(pi/(2n))
        n       *= 2
    # n*sin_half = n * sin(pi/n) -> pi
    return n * sin_half

# ── Auswertung ─────────────────────────────────────────────────────────
PI_REF = math.pi

print("=" * 60)
print("Pi-Berechnung auf vier unabhaengigen Wegen")
print("=" * 60)

for name, func in [
    ("Pfad A: Leibniz-Reihe  (10^7 Terme)",      leibniz),
    ("Pfad B: Monte-Carlo    (5*10^6 Punkte)",    monte_carlo),
    ("Pfad C: Nilakantha     (5*10^5 Terme)",     nilakantha),
    ("Pfad D: Archimedes H¹ (46 Verdopplungen, stabil)",  archimedes_hausdorff),
]:
    pi_val = func()
    err    = abs(pi_val - PI_REF)
    print(f"\n  {name}")
    print(f"  pi ~ {pi_val:.10f}")
    print(f"  Fehler: {err:.2e}  ({100*err/PI_REF:.6f}%)")

print(f"\n  Referenz (math.pi): {PI_REF:.10f}")
print("=" * 60)
