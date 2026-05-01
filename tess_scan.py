import math

chi = math.radians(59.1)
rho = 0.406
Df  = math.log(2) / math.log(1/rho)
Va  = 1 - math.cos(chi)   # asymmetrischer Anteil (versetzt)
Vs  = math.cos(chi)        # symmetrischer Anteil (Überlapp)
denom = (2 * math.pi**2) ** (6 + Df)
eta_obs = 6.104e-10

print(f"chi={math.degrees(chi):.1f}°  rho={rho}  Df_geo={Df:.4f}")
print(f"Va = 1-cos(chi) = {Va:.4f}  (asymm. / versetzter Teil)")
print(f"Vs = cos(chi)   = {Vs:.4f}  (symm. Überlapp)")
print(f"Denom (2pi^2)^(6+Df) = {denom:.4e}")
print()

cases = [
    ("Va*(1-rho)",         Va*(1-rho)),
    ("Vs*(1-rho)",         Vs*(1-rho)),
    ("Va*rho*(1-rho)",     Va*rho*(1-rho)),
    ("2*Va*rho*(1-rho)",   2*Va*rho*(1-rho)),
    ("Va*Vs",              Va*Vs),
    ("Va*Vs*(1-rho)",      Va*Vs*(1-rho)),
    ("Va^2",               Va**2),
    ("Vs*(1-Vs)",          Vs*(1-Vs)),
    ("Va/(1+Va)",          Va/(1+Va)),
    ("(1-rho)/(1+rho)*Va", (1-rho)/(1+rho)*Va),
    ("sin^2(chi/2)",       math.sin(chi/2)**2),
    ("sin^2*rho*(1-rho)",  math.sin(chi/2)**2 * rho*(1-rho)),
    ("(1-cos^2(chi))",     1 - math.cos(chi)**2),
    ("sin(chi)*(1-rho)",   math.sin(chi)*(1-rho)),
    ("sin(chi)*rho*(1-rho)",math.sin(chi)*rho*(1-rho)),
]

print(f"{'Zähler':<30}  {'eta':>12}  {'ratio':>8}")
print("-"*60)
for label, num in cases:
    eta = num / denom
    ratio = eta / eta_obs
    marker = " <-- TREFFER" if 0.8 < ratio < 1.2 else ""
    print(f"{label:<30}  {eta:>12.4e}  {ratio:>8.3f}{marker}")
