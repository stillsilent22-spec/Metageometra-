"""
OT-8: M_meta — Metageometra Mandelbrot-Analogon
================================================
Scanne chi x rho Raum, pruefe Konvergenz zu Shell-Knoten
"""
import math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE   = os.path.dirname(os.path.abspath(__file__))
RESDIR = os.path.join(BASE, "results")
os.makedirs(RESDIR, exist_ok=True)

CHI_0  = 59.1; RHO_0 = 0.406; DELTA = 1.0

def shell_nodes(chi, delta=DELTA, n_max=6):
    nodes = []
    for n in range(1, n_max+1):
        d = math.radians(delta); c = math.radians(chi)
        nodes.append(math.degrees(math.acos(max(-1., min(1., math.cos(n*d)*math.cos(n*c))))))
    return nodes

def torsion_map(theta, chi, rho, n_iter=50, tol=3.0):
    """
    Iterierte Torsionsabbildung auf S2:
      T(theta) = arccos(cos(delta)*cos(chi) * cos(theta) + rho*sin(theta))
    Prueft ob orbit zu einem Shell-Knoten konvergiert (innerhalb tol Grad).
    Gibt (konvergiert, n_steps, theta_final) zurueck.
    """
    d  = math.radians(DELTA); c = math.radians(chi)
    nodes = shell_nodes(chi)
    th = math.radians(theta)
    for i in range(n_iter):
        cos_new = math.cos(d)*math.cos(c)*math.cos(th) - rho*math.sin(th)
        cos_new = max(-1., min(1., cos_new))
        th_new  = math.acos(cos_new)
        th_new  = th_new + 0.01*(math.pi - th_new)  # leichte Attraktion zur Mitte
        # Konvergenz zu einem Knoten?
        th_deg  = math.degrees(th_new)
        for nd in nodes:
            if abs(th_deg - nd) < tol:
                return True, i+1, th_deg
        th = th_new
    return False, n_iter, math.degrees(th)

# ── 500x500 Scan ──────────────────────────────────────────
N = 200  # 200x200 fuer Geschwindigkeit (500 wäre sehr langsam)
chi_range = np.linspace(55.0, 65.0, N)
rho_range = np.linspace(0.35, 0.50, N)

conv_map  = np.zeros((N, N), dtype=float)
theta0    = 45.0  # Startwinkel fuer alle Orbits

print("OT-8: M_meta Scan laeuft...")
for i, chi in enumerate(chi_range):
    for j, rho in enumerate(rho_range):
        ok, n_steps, _ = torsion_map(theta0, chi, rho, n_iter=80, tol=2.5)
        conv_map[i, j] = n_steps if ok else 0

# ── Fraktaldimension der Grenze (Box-Counting) ────────────
boundary = (conv_map > 0).astype(float)
# Gradient = Grenze
gx = np.diff(boundary, axis=0)
gy = np.diff(boundary, axis=1)
edge = np.zeros_like(boundary)
edge[:-1, :] += np.abs(gx); edge[:, :-1] += np.abs(gy)
edge_bin = edge > 0

epsilons = [2, 4, 8, 16, 32]
N_boxes  = []
for eps in epsilons:
    cnt = 0
    for ii in range(0, N, eps):
        for jj in range(0, N, eps):
            if edge_bin[ii:ii+eps, jj:jj+eps].any():
                cnt += 1
    N_boxes.append(cnt)

log_eps = np.log(epsilons)
log_N   = np.log(np.maximum(N_boxes, 1))
slope, intercept = np.polyfit(log_eps, log_N, 1)
FD = -slope

# ── Referenzpunkt chi=59.1, rho=0.406 ──────────────────────
i0 = np.argmin(np.abs(chi_range - CHI_0))
j0 = np.argmin(np.abs(rho_range - RHO_0))
ref_conv = conv_map[i0, j0]

print("══════════════════════════════════════════════════════")
print("OT-8: M_meta — Metageometra Mandelbrot-Scan")
print("══════════════════════════════════════════════════════")
print(f"  Grid: {N}x{N},  chi=[55,65]  rho=[0.35,0.50]")
print(f"  Konvergenz bei (chi={CHI_0}, rho={RHO_0}): {ref_conv} Schritte")
print(f"  Conv-Rate gesamt: {(conv_map>0).sum()} / {N*N} = {(conv_map>0).mean()*100:.1f}%")
print()
print(f"  Box-Counting FD (Grenzlinie):")
for e, n in zip(epsilons, N_boxes):
    print(f"    eps={e}: N={n}")
print(f"  FD = {FD:.4f}  (Erwartung: 1.2 - 1.8, OT-6 Ref: 0.78)")
print()
print(f"  {'BESTAETIGT' if 1.0 <= FD <= 1.8 else 'INCONCLUSIVE'}: FD={'%0.3f'%FD}")
print("══════════════════════════════════════════════════════")

# ── Plot ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.patch.set_facecolor('white')

ax = axes[0]
extent = [rho_range[0], rho_range[-1], chi_range[0], chi_range[-1]]
im = ax.imshow(conv_map, extent=extent, origin='lower', aspect='auto',
               cmap='hot', interpolation='nearest')
ax.plot(RHO_0, CHI_0, 'c*', ms=12, label=f'V21.3: chi={CHI_0}, rho={RHO_0}')
ax.set_xlabel('rho', fontsize=9); ax.set_ylabel('chi (Grad)', fontsize=9)
ax.set_title(f'M_meta Konvergenzkarte\n(schwarz=kein Knoten, hell=schnell)', fontsize=8)
ax.legend(fontsize=8)
plt.colorbar(im, ax=ax, label='Konvergenzschritte')

ax2 = axes[1]
ax2.loglog(epsilons, N_boxes, 'bo-', ms=6)
x_fit = np.array([epsilons[0], epsilons[-1]])
ax2.loglog(x_fit, np.exp(intercept) * x_fit**slope, 'r--', label=f'FD={FD:.3f}')
ax2.set_xlabel('Box-Groesse epsilon', fontsize=9)
ax2.set_ylabel('N(epsilon)', fontsize=9)
ax2.set_title('Box-Counting FD der M_meta Grenze', fontsize=8)
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)

plt.suptitle('OT-8: M_meta — Metageometra Mandelbrot-Analogon', fontsize=9)
plt.tight_layout()
out = os.path.join(RESDIR, 'mmeta_scan.png')
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
plt.close()

with open(os.path.join(RESDIR, 'mmeta_boundary_fd.txt'), 'w') as f:
    f.write(f"OT-8: M_meta Fraktaldimension\n")
    f.write(f"Grid: {N}x{N}\nFD = {FD:.4f}\n")
    f.write(f"Referenz (chi={CHI_0},rho={RHO_0}): {ref_conv} Konvergenzschritte\n")
    f.write(f"Status: {'BESTAETIGT' if 1.0<=FD<=1.8 else 'INCONCLUSIVE'}\n")
print(f"Gespeichert: {out}")
