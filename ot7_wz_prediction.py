"""
OT-7: w(z) DESI DR2 Fit (V21.3)
================================
w(z) = -1 + (1-D_f/2) * delta_w * (1+z)^(3*delta_w), D_f=1.8
"""
import math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

BASE   = os.path.dirname(os.path.abspath(__file__))
RESDIR = os.path.join(BASE, "results")
os.makedirs(RESDIR, exist_ok=True)

D_F = 1.8
PREFACTOR = 1 - D_F/2  # = 0.1

# ── DESI DR2 Datenpunkte ─────────────────────────────────
desi_z  = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 1.49])
desi_w  = np.array([-0.87, -0.93, -1.07, -0.95, -0.98, -1.01])
desi_sw = np.array([0.15,  0.09,  0.08,  0.10,  0.12,  0.14])

# DESI DR2 Best-fit CPL: w0=-0.827, wa=-0.75
CPL_w0 = -0.827; CPL_wa = -0.75

def w_htm(z, dw):
    return -1.0 + PREFACTOR * dw * (1+z)**(3*dw)

def w_cpl(z, w0=CPL_w0, wa=CPL_wa):
    return w0 + wa * z/(1+z)

# ── Fit delta_w ──────────────────────────────────────────
def chi2_htm(dw):
    wmod = np.array([w_htm(z, dw) for z in desi_z])
    return np.sum(((wmod - desi_w)/desi_sw)**2)

res = minimize_scalar(chi2_htm, bounds=(-3, -1e-6), method='bounded')
dw_best = res.x
chi2_htm_best = res.fun
dof = len(desi_z) - 1

wmod_best = np.array([w_htm(z, dw_best) for z in desi_z])
rms_htm   = np.sqrt(np.mean((wmod_best - desi_w)**2))

wmod_cpl  = np.array([w_cpl(z) for z in desi_z])
rms_cpl   = np.sqrt(np.mean((wmod_cpl - desi_w)**2))
chi2_cpl  = np.sum(((wmod_cpl - desi_w)/desi_sw)**2)

z_plot = np.linspace(0, 3, 300)
w_htm_plot = np.array([w_htm(z, dw_best) for z in z_plot])
w_cpl_plot = np.array([w_cpl(z) for z in z_plot])

print("══════════════════════════════════════════════════════")
print("OT-7: w(z) DESI DR2 Fit (V21.3)")
print("══════════════════════════════════════════════════════")
print(f"  D_f = {D_F}, Prefactor (1-D_f/2) = {PREFACTOR}")
print(f"  HTM Best-fit: delta_w = {dw_best:.5f}")
print(f"  w(z=0)  = {w_htm(0, dw_best):.4f}")
print(f"  w(z=0.5)= {w_htm(0.5, dw_best):.4f}")
print(f"  w(z=1)  = {w_htm(1, dw_best):.4f}")
print(f"  Chi2/dof HTM = {chi2_htm_best/dof:.4f}")
print(f"  RMS HTM = {rms_htm:.4f}")
print(f"  Chi2/dof CPL = {chi2_cpl/dof:.4f}")
print(f"  RMS CPL = {rms_cpl:.4f}")
print(f"  DESI DR2 CPL: w0={CPL_w0}, wa={CPL_wa}")
print()
cpl_better = rms_cpl < rms_htm
print(f"  CPL besser: {cpl_better}")
status = "BESTAETIGT" if chi2_htm_best/dof < 3 else "INCONCLUSIVE"
print(f"  Status: {status} (Chi2/dof={chi2_htm_best/dof:.2f})")
print("══════════════════════════════════════════════════════")

# ── Plot ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
ax.fill_between(z_plot, w_cpl_plot-0.05, w_cpl_plot+0.05, alpha=0.1, color='orange')
ax.plot(z_plot, w_cpl_plot, 'orange', lw=1.5, ls='--', label=f'CPL DESI DR2 (w0={CPL_w0}, wa={CPL_wa})')
ax.plot(z_plot, w_htm_plot, 'blue',   lw=2.0, label=f'HTM D_f={D_F}: dw={dw_best:.4f}')
ax.axhline(-1, color='k', lw=0.8, ls=':', label='w = -1 (Lambda CDM)')
ax.errorbar(desi_z, desi_w, yerr=desi_sw, fmt='rs', ms=6, capsize=3, label='DESI DR2 Datenpunkte')
ax.set_xlabel('Rotverschiebung z', fontsize=9)
ax.set_ylabel('w(z)', fontsize=9)
ax.set_title(f'OT-7: w(z) HTM vs DESI DR2\nD_f={D_F}, dw={dw_best:.4f}, Chi2/dof={chi2_htm_best/dof:.2f}', fontsize=9)
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_xlim(0, 2.5); ax.set_ylim(-1.4, -0.6)
plt.tight_layout()
out = os.path.join(RESDIR, 'w_z_vs_desi.png')
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
plt.close()

lines = [
    "════════════════════════════════════════════════════════════",
    "OT-7: w(z) DESI DR2 Fit (V21.3)",
    "════════════════════════════════════════════════════════════",
    f"HTM: delta_w = {dw_best:.5f}",
    f"Chi2/dof HTM = {chi2_htm_best/dof:.4f}",
    f"RMS HTM = {rms_htm:.4f}",
    f"RMS CPL = {rms_cpl:.4f}",
    f"w(z=0) = {w_htm(0,dw_best):.4f}",
    f"Status: {status}",
    "════════════════════════════════════════════════════════════",
]
with open(os.path.join(RESDIR,'OT_7_result.txt'),'w',encoding='utf-8') as f:
    f.write('\n'.join(lines)+'\n')
print(f"Gespeichert: {out}")
