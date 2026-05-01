"""
OT-39: 97-SMBH Katalog — 4-Arm KS-Test + Box-Counting FD (V21.3)
==================================================================
Shell-Test + 4-Arm Azimut-KS-Test + Monte Carlo
"""
import math, os
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, kstest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE   = os.path.dirname(os.path.abspath(__file__))
CATS   = os.path.join(BASE, "results", "catalogs")
RESDIR = os.path.join(BASE, "results")
os.makedirs(RESDIR, exist_ok=True)

CHI   = 59.1; DELTA = 1.0; RHO = 0.406
DPOL  = {'l': 305.0, 'b': 25.0}
BASE_AZ = 88.0; TILT = CHI/3

def shell_theta(n):
    d = math.radians(DELTA); c = math.radians(CHI)
    return math.degrees(math.acos(max(-1., min(1., math.cos(n*d)*math.cos(n*c)))))

shells = np.array([shell_theta(n) for n in range(1, 7)])

def ang_dist_az(l, b):
    bP = math.radians(DPOL['b']); lP = math.radians(DPOL['l'])
    bO = math.radians(b);          lO = math.radians(l)
    cos_d = (math.sin(bP)*math.sin(bO) + math.cos(bP)*math.cos(bO)*math.cos(lP-lO))
    dist  = math.degrees(math.acos(max(-1, min(1, cos_d))))
    y = math.cos(bO)*math.sin(lO-lP)
    x = math.cos(bP)*math.sin(bO) - math.sin(bP)*math.cos(bO)*math.cos(lO-lP)
    az = math.degrees(math.atan2(y, x)) % 360
    return dist, az

# ── Katalog laden ─────────────────────────────────────────
df = pd.read_csv(os.path.join(CATS, "smbh_extended.csv"), comment='#')
n_tot = len(df)

thetas, azs = [], []
for _, row in df.iterrows():
    # Convert RA/Dec to galactic (approx via RA/Dec directly if no galactic cols)
    # Use RA_deg, Dec_deg via equatorial D-pol conversion
    ra  = float(row['RA_deg']); dec = float(row['Dec_deg'])
    # Equatorial D-Pol from galactic (305,25) -> already computed in viz script
    # Use same equatorial dist_pa logic
    RA_NGP=192.859508; Dec_NGP=27.128336; l_NCP=122.932
    l_dp=305.; b_dp=25.
    _sin_dec = (math.sin(math.radians(b_dp))*math.sin(math.radians(Dec_NGP))
                + math.cos(math.radians(b_dp))*math.cos(math.radians(Dec_NGP))
                * math.cos(math.radians(l_NCP-l_dp)))
    DP_Dec = math.degrees(math.asin(max(-1.,min(1.,_sin_dec))))
    _cos_dec= math.cos(math.radians(DP_Dec))
    _sn = -math.cos(math.radians(b_dp))*math.sin(math.radians(l_dp-l_NCP))/_cos_dec
    _co = (math.sin(math.radians(b_dp))-math.sin(math.radians(Dec_NGP))*_sin_dec)/(math.cos(math.radians(Dec_NGP))*_cos_dec)
    DP_RA = math.degrees(math.atan2(_sn,_co)) + RA_NGP

    cos_a = (math.sin(math.radians(DP_Dec))*math.sin(math.radians(dec))
             + math.cos(math.radians(DP_Dec))*math.cos(math.radians(dec))
             * math.cos(math.radians(DP_RA-ra)))
    th = math.degrees(math.acos(max(-1.,min(1.,cos_a))))
    dra = ra - DP_RA
    y   = math.sin(math.radians(dra))*math.cos(math.radians(dec))
    x   = math.cos(math.radians(DP_Dec))*math.sin(math.radians(dec)) - math.sin(math.radians(DP_Dec))*math.cos(math.radians(dec))*math.cos(math.radians(dra))
    az  = math.degrees(math.atan2(y,x)) % 360
    thetas.append(th); azs.append(az)

thetas = np.array(thetas); azs = np.array(azs)

# ── Shell-Test (±2°) ──────────────────────────────────────
tol2 = 2.0
min_dist = np.array([np.min(np.abs(th - shells)) for th in thetas])
hits2  = (min_dist < tol2).sum()
p_exp2 = sum(min(2*tol2, 180) for _ in shells) / 180  # coverage fraction
p_exp2 = min(1.0, p_exp2)
from scipy.stats import binom
p_binom2 = 1 - binom.cdf(hits2-1, n_tot, p_exp2)
z_shell  = (hits2 - n_tot*p_exp2) / math.sqrt(n_tot*p_exp2*(1-p_exp2))

# ── 4-Arm KS-Test: OT-42 ─────────────────────────────────
arm_azimuths = np.array([BASE_AZ, BASE_AZ+TILT, BASE_AZ+180, BASE_AZ+180+TILT]) % 360
bandwidth = TILT / 2  # ±9.85 Grad

# Min-Abstand jedes SMBH zum naechsten Arm-Azimut
def min_az_dist(az, arms):
    diffs = np.abs(az - arms)
    diffs = np.minimum(diffs, 360 - diffs)
    return diffs.min()

az_dists = np.array([min_az_dist(az, arm_azimuths) for az in azs])
in_arm   = (az_dists < bandwidth).sum()
p_arm_exp = 4 * 2*bandwidth / 360
p_arm_binom = 1 - binom.cdf(in_arm-1, n_tot, p_arm_exp)
z_arm = (in_arm - n_tot*p_arm_exp) / math.sqrt(n_tot*p_arm_exp*(1-p_arm_exp))

# KS-Test: Azimutverteilung vs. uniform
ks_stat, ks_p = kstest(azs/360, 'uniform')

# ── Monte Carlo: Zufaellige D-Pole (100 Permutationen) ────
np.random.seed(42)
n_mc = 1000
z_mc = []
for _ in range(n_mc):
    # Zufaelliger Pol auf der Sphaere
    l_rand = np.random.uniform(0, 360)
    b_rand = np.degrees(math.asin(np.random.uniform(-1, 1)))
    th_rand = []
    for _, row in df.iterrows():
        ra_r = float(row['RA_deg']); dec_r = float(row['Dec_deg'])
        cos_a = (math.sin(math.radians(b_rand))*math.sin(math.radians(dec_r))
                 + math.cos(math.radians(b_rand))*math.cos(math.radians(dec_r))
                 * math.cos(math.radians(l_rand-ra_r)))
        th_rand.append(math.degrees(math.acos(max(-1.,min(1.,cos_a)))))
    th_rand = np.array(th_rand)
    md = np.array([np.min(np.abs(t-shells)) for t in th_rand])
    k_rand = (md < tol2).sum()
    zr = (k_rand - n_tot*p_exp2) / math.sqrt(n_tot*p_exp2*(1-p_exp2))
    z_mc.append(zr)
z_mc = np.array(z_mc)
p_mc = (z_mc >= z_shell).mean()

# ── Box-Counting FD der Azimutverteilung ─────────────────
az_bin = np.zeros(360)
for az in azs:
    az_bin[int(az) % 360] += 1
epsilons = [2, 5, 10, 20, 45]
N_boxes  = []
for eps in epsilons:
    cnt = sum(1 for i in range(0, 360, eps) if az_bin[i:i+eps].sum() > 0)
    N_boxes.append(cnt)
log_e = np.log(epsilons); log_n = np.log(np.maximum(N_boxes, 1))
slope_fd, _ = np.polyfit(log_e, log_n, 1)
FD_az = -slope_fd

print("══════════════════════════════════════════════════════")
print("OT-39: 97-SMBH — Shell-Test + 4-Arm KS-Test (V21.3)")
print("══════════════════════════════════════════════════════")
print(f"  N_SMBHs = {n_tot}")
print(f"  Shells: {[f'{s:.1f}' for s in shells]} Grad")
print()
print(f"  ── Shell-Test (tol=±{tol2}°) ──")
print(f"  Treffer:    {hits2}/{n_tot} ({hits2/n_tot*100:.1f}%)")
print(f"  p_exp:      {p_exp2:.4f} ({p_exp2*100:.1f}%)")
print(f"  z-Wert:     {z_shell:.3f} sigma")
print(f"  p(Binom):   {p_binom2:.6f}")
print(f"  MC p-Wert:  {p_mc:.4f}  ({n_mc} Permutationen)")
print()
print(f"  ── OT-42 4-Arm KS-Test ──")
print(f"  Arm-Azimuthe: {[f'{a:.1f}' for a in arm_azimuths]}")
print(f"  Bandwidth:    ±{bandwidth:.2f}°")
print(f"  In-Arm:       {in_arm}/{n_tot}")
print(f"  p_arm_exp:    {p_arm_exp:.4f}")
print(f"  z_arm:        {z_arm:.3f} sigma")
print(f"  p_arm_binom:  {p_arm_binom:.6f}")
print()
print(f"  KS-Test Azimut vs Uniform: D={ks_stat:.4f}, p={ks_p:.4f}")
print()
print(f"  Box-Counting FD (Azimut): {FD_az:.3f}")
print()
status_shell = "BESTAETIGT" if z_shell > 3 else ("INCONCLUSIVE" if z_shell > 1.5 else "OFFEN")
status_arm   = "BESTAETIGT" if z_arm   > 2 else "INCONCLUSIVE"
print(f"  Shell-Status: {status_shell} (z={z_shell:.2f})")
print(f"  Arm-Status:   {status_arm} (z={z_arm:.2f})")
print("══════════════════════════════════════════════════════")

# ── Plot ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.patch.set_facecolor('white')

ax = axes[0]
ax.hist(thetas, bins=36, range=(0,180), density=True, color='steelblue', alpha=0.7, label='SMBHs')
for s in shells:
    ax.axvline(s, color='red', lw=1.5, ls='--', alpha=0.8)
ax.set_xlabel('theta_D (Grad)', fontsize=8); ax.set_ylabel('Dichte', fontsize=8)
ax.set_title(f'Shell-Clustering\nz={z_shell:.2f}σ, {hits2}/{n_tot} Treffer (±{tol2}°)', fontsize=8)
ax.legend(fontsize=7); ax.grid(alpha=0.3)

ax2 = axes[1]
ax2.hist(azs, bins=36, range=(0,360), density=True, color='green', alpha=0.7, label='SMBHs')
for arm_az in arm_azimuths:
    ax2.axvline(arm_az, color='red', lw=1.5, ls='--')
    ax2.axvspan(arm_az-bandwidth, arm_az+bandwidth, alpha=0.1, color='red')
ax2.set_xlabel('Azimut (Grad)', fontsize=8); ax2.set_ylabel('Dichte', fontsize=8)
ax2.set_title(f'4-Arm KS-Test (OT-42)\nz_arm={z_arm:.2f}σ, {in_arm}/{n_tot} in Arm', fontsize=8)
ax2.legend(fontsize=7); ax2.grid(alpha=0.3)

ax3 = axes[2]
ax3.hist(z_mc, bins=30, density=True, color='gray', alpha=0.7, label='Random poles MC')
ax3.axvline(z_shell, color='red', lw=2, label=f'D-Pol z={z_shell:.2f}')
ax3.set_xlabel('z-Wert', fontsize=8); ax3.set_ylabel('Dichte', fontsize=8)
ax3.set_title(f'Monte Carlo: D-Pol vs {n_mc} Zufallspole\np_mc={p_mc:.3f}', fontsize=8)
ax3.legend(fontsize=7); ax3.grid(alpha=0.3)

plt.suptitle('OT-39/42: SMBH Shell + 4-Arm Test', fontsize=9)
plt.tight_layout()
out = os.path.join(RESDIR, 'ot39_ks_test.png')
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
plt.close()

lines = [
    "════════════════════════════════════════════════════════════",
    "OT-39/42: SMBH Shell-Test + 4-Arm KS-Test (V21.3)",
    "════════════════════════════════════════════════════════════",
    f"N_SMBHs = {n_tot}",
    f"Shell-Test (tol=±{tol2}°): {hits2}/{n_tot}, z={z_shell:.3f}, p={p_binom2:.6f}",
    f"MC p-Wert ({n_mc} Perms): {p_mc:.4f}",
    f"4-Arm Test: {in_arm}/{n_tot}, z={z_arm:.3f}, p={p_arm_binom:.6f}",
    f"KS D={ks_stat:.4f}, p={ks_p:.4f}",
    f"FD Azimut: {FD_az:.3f}",
    f"Shell-Status: {status_shell}",
    f"Arm-Status: {status_arm}",
    "════════════════════════════════════════════════════════════",
]
with open(os.path.join(RESDIR,'OT_39_result.txt'),'w',encoding='utf-8') as f:
    f.write('\n'.join(lines)+'\n')
print(f"Gespeichert: {out}")
