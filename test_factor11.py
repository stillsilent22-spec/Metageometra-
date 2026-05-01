"""Faktor-11-Test stratifiziert nach Entfernungs-Bins.
Prueft ob der Effekt tautologisch (verschwindet) oder real (bleibt) ist.
"""
import csv, math
import numpy as np
import scipy.stats as st

rows = list(csv.DictReader(open('results/catalogs/smbh_extended.csv', encoding='utf-8')))

data = []
for r in rows:
    if not r.get('dist_mpc', '').strip():
        continue
    d     = float(r['dist_mpc'])
    delta = float(r['shell_delta'])
    phys  = d * math.sin(math.radians(delta))
    hit   = r['hit_5deg'] == 'YES'
    data.append({'d': d, 'delta': delta, 'phys': phys, 'hit': hit, 'name': r['Name']})

BINS = [
    ('<  15 Mpc',   lambda d: d <  15),
    ('15-50 Mpc',   lambda d: 15 <= d <  50),
    ('50-200 Mpc',  lambda d: 50 <= d < 200),
    ('>200 Mpc',    lambda d: d >= 200),
]

print()
print('DISTANCE-BIN STRATIFIED ANALYSIS')
print('Physikalischer Schalenabstand = d x sin(delta_theta)')
print()
print(f'  {"Bin":<12} {"n_hit":>6} {"n_miss":>6}  '
      f'{"hit_mean":>9}  {"miss_mean":>9}  {"Faktor":>7}  {"p(MW)":>8}  Signal')
print('  ' + '-'*72)

results = []
for label, cond in BINS:
    sub    = [x for x in data if cond(x['d'])]
    h_vals = np.array([x['phys'] for x in sub if     x['hit']])
    m_vals = np.array([x['phys'] for x in sub if not x['hit']])
    n_h, n_m = len(h_vals), len(m_vals)
    if n_h < 2 or n_m < 2:
        print(f'  {label:<12} {n_h:>6} {n_m:>6}   -- zu wenig Daten (n_hit={n_h}, n_miss={n_m})')
        results.append((label, n_h, n_m, None, None, None, None))
        continue
    fac  = m_vals.mean() / h_vals.mean() if h_vals.mean() > 0 else float('inf')
    _, p = st.mannwhitneyu(h_vals, m_vals, alternative='less')
    sig  = '*** ECHT' if p < 0.05 else 'n.s.'
    print(f'  {label:<12} {n_h:>6} {n_m:>6}  '
          f'{h_vals.mean():>9.2f}  {m_vals.mean():>9.2f}  {fac:>7.1f}x  {p:>8.4f}  {sig}')
    results.append((label, n_h, n_m, h_vals.mean(), m_vals.mean(), fac, p))
    # Detail-Liste
    for x in sorted(sub, key=lambda z: z['d']):
        tag = 'HIT' if x['hit'] else '   '
        print(f'      {tag}  {x["name"]:<22} '
              f'd={x["d"]:6.1f}  delta={x["delta"]:5.1f}deg  phys={x["phys"]:6.2f} Mpc')
    print()

# Gesamt zum Vergleich
all_h = np.array([x['phys'] for x in data if     x['hit']])
all_m = np.array([x['phys'] for x in data if not x['hit']])
_, p_all = st.mannwhitneyu(all_h, all_m, alternative='less')
fac_all  = all_m.mean() / all_h.mean()
print(f'  {"GESAMT":<12} {len(all_h):>6} {len(all_m):>6}  '
      f'{all_h.mean():>9.2f}  {all_m.mean():>9.2f}  {fac_all:>7.1f}x  {p_all:>8.5f}  '
      f'*** SIGN.' if p_all < 0.05 else 'n.s.')

print()
print('INTERPRETATION:')
print('  Tautologisch = Effekt verschwindet innerhalb der Bins (n.s. in allen Bins)')
print('  Echtes Signal = Effekt bleibt auch innerhalb der Bins (sign. in mind. einem Bin)')
sig_bins = [r for r in results if r[6] is not None and r[6] < 0.05]
print()
if sig_bins:
    print(f'  ERGEBNIS: Effekt BLEIBT in {len(sig_bins)} Bin(s) signifikant -> ECHTES SIGNAL')
    for r in sig_bins:
        print(f'    Bin "{r[0]}": p={r[6]:.4f}, Faktor={r[5]:.1f}x')
else:
    print('  ERGEBNIS: Effekt verschwindet in allen Bins -> TAUTOLOGISCH (reiner Distanz-Bias)')

