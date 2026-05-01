"""Milliquas theta_D Verteilung prüfen"""
import csv, math

def dpole_ra_dec():
    ra_ngp = math.radians(192.85948)
    dec_ngp = math.radians(27.12825)
    l_ncp  = math.radians(122.93192)
    l_d = math.radians(305.0); b_d = math.radians(25.0)
    sin_b = math.sin(b_d); cos_b = math.cos(b_d)
    sin_l = math.sin(l_d - l_ncp); cos_l = math.cos(l_d - l_ncp)
    sin_dec_ngp = math.sin(dec_ngp); cos_dec_ngp = math.cos(dec_ngp)
    dec_dp = math.asin(sin_dec_ngp * sin_b + cos_dec_ngp * cos_b * sin_l)
    ra_dp  = math.atan2(cos_b * cos_l, sin_b * cos_dec_ngp - cos_b * sin_dec_ngp * sin_l) + ra_ngp
    return math.degrees(ra_dp) % 360, math.degrees(dec_dp)

ra_d, dec_d = dpole_ra_dec()
print(f"D-Pol RA={ra_d:.2f}  Dec={dec_d:.2f}")

ra_d_r = math.radians(ra_d); dec_d_r = math.radians(dec_d)
xd = math.cos(dec_d_r)*math.cos(ra_d_r)
yd = math.cos(dec_d_r)*math.sin(ra_d_r)
zd = math.sin(dec_d_r)

def theta_from_radec(ra_deg, dec_deg):
    ra_r = math.radians(ra_deg); dec_r = math.radians(dec_deg)
    xp = math.cos(dec_r)*math.cos(ra_r)
    yp = math.cos(dec_r)*math.sin(ra_r)
    zp = math.sin(dec_r)
    return math.degrees(math.acos(max(-1., min(1., xd*xp + yd*yp + zd*zp))))

import os
base = os.path.dirname(os.path.abspath(__file__))
fp   = os.path.join(base, "results", "milliquas_sample.csv")

thetas = []
dec_vals = []
with open(fp, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            ra  = float(row['RAJ2000'])
            dec = float(row['DEJ2000'])
            thetas.append(theta_from_radec(ra, dec))
            dec_vals.append(dec)
        except (KeyError, ValueError):
            pass

print(f"Total: {len(thetas)}")
print(f"Dec range: {min(dec_vals):.1f} ... {max(dec_vals):.1f}")

bins = [0] * 9
for t in thetas:
    bins[min(8, int(t // 20))] += 1
print("Theta-Bins (je 20 Grad):")
for i, b in enumerate(bins):
    lo, hi = i*20, (i+1)*20
    # isotrope Erwartung: f = (cos(lo)-cos(hi))/2 * N
    import math as m
    f = (m.cos(m.radians(lo)) - m.cos(m.radians(hi))) / 2.0 * len(thetas)
    print(f"  {lo:3d}-{hi:3d}: obs={b:4d}  iso_exp={f:6.1f}  ratio={b/max(1,f):.2f}")
