#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
METAGEOMETRA V21.3 — 2-Arm Helix Kreuzversatz  (v3)
=====================================================
KORREKTE STRUKTUR: 2 Arme (nicht 4), jeder Arm ist ein vollstaendiger
Kreisbahn-Pfad der Kugeloberfläche der BEIDE Hemisphären verbindet.

Kreuzversatz:
  Arm A (prograd):  near-side Az=88°,  far-side Az=268°
  Arm B (retrograd): near-side Az=107.7°, far-side Az=287.7°
  Versatz = chi/3 = 19.7° zwischen Arm A und Arm B.

An jedem Kehrtwende-Punkt (theta≈8°, nahe D-Pol, alle 6 Schritte)
wechselt der Arm seitenweise: near-side → far-side → near-side → ...
Der Kreuzungspunkt liegt beim D-Pol (theta≈8°), NICHT im Kugelzentrum.

Raum-Expansion (OT-43):
  eps_exp = (1-cos(Δchi))·(1-2ρ) = 2.31e-5 pro Zyklus
  Schritt n liegt auf Sphäre R(n) = (1 + eps_exp)^(n//6)
  Nach 1200 Zyklen: R ≈ 1.028  (sichtbare Expansion)

Kurz-Kurz-Lang-Kurz Doublett-Rhythmus:
  Schalen n=1 (59°) und n=5 (64°) → Doublett KURZ (nahe beieinander!)
  Schalen n=2 (118°) und n=4 (123°) → Doublett KURZ
  Schale  n=3 (176°) → anti-D-Pol → LANG (weit entfernt!)
  Schale  n=6 (8°)   → D-Pol Kehrtwende → KURZ

Steuerung:
  LMB + Drag  Kugel drehen
  Mausrad     Zoom
  SPACE       Pause / Weiter
  Pfeil auf/ab Geschwindigkeit +/-
  R           Reset
  [ / ]       chi +/- 0.1 deg (live)
  G           Gitter + Schalenringe
  S           Sphärenflaeche
  D           D-Pol-Achse
  E           Expansion ein/aus
  1 / 2       Arm A / Arm B einzeln
  ESC         Beenden
"""

import sys
import math
import numpy as np
import pygame
from pygame.locals import DOUBLEBUF, OPENGL
from OpenGL.GL import *
from OpenGL.GLU import *

# ======================================================================
# Framework-Konstanten
# ======================================================================
CHI   = 59.1     # Tesserakt-Rotationswinkel [deg]  (live via [ / ])
DELTA = 1.0      # Shell-Inkrement [deg]

ETA     = 5.58e-10
EPS_EXP = 2.31e-5   # Expansion pro Zyklus (OT-43)

DP_L, DP_B = 305.0, 25.0   # D-Pol galaktisch
BASE_ROT   = 88.0           # Basisazimut Arm A near-side

ARM_COLORS = [
    (0.25, 0.65, 1.00),   # Arm A  Blau   (prograd)
    (1.00, 0.30, 0.08),   # Arm B  Orange (retrograd)
]
ARM_NAMES = [
    "Arm A  prograd    near=88°/far=268°",
    "Arm B  retrograd  near=107.7°/far=287.7°",
]

# Kurz-Kurz-Lang-Kurz Rhythmus: n=1/5 Doublett, n=2/4 Doublett, n=3 LANG, n=6 Kehrtwende
SHELL_RHYTHM = ["KURZ", "KURZ", "LANG", "KURZ", "KURZ", "KURZ"]
# Doublett-Paare: (ni, nj) sind nahe beieinander auf der Kugeloberfläche
DOUBLET_PAIRS = [(0, 4), (1, 3)]   # 0-indexed: n=1/n=5 und n=2/n=4

# (Name, theta_D [deg], az_D [deg], RGB-float, Referenz, Arm, Shell-n)
CONFIRMED_SMBHS = [
    ("Sgr A*",   58.65, 109.0, (1.00, 1.00, 0.15), "EHT 2022",    "A", 1),
    ("M31",      175.1, 270.0, (0.75, 0.75, 1.00), "Local Group", "A", 3),
    ("NGC 1052",  97.0, 270.0, (1.00, 0.55, 0.05), "Baczko 2016", "B", 2),
    ("NGC 0315",  80.0, 265.0, (0.15, 1.00, 0.50), "Daly 2023",   "B", 1),
    ("NGC 4751",   5.1, 105.0, (0.90, 0.40, 1.00), "OT-5 cand.",  "A", 6),
]

TOTAL_STEPS  = 7200   # 1200 Zyklen
SUBSTEPS     = 10     # Geodaetische Zwischenpunkte
TOTAL_SMOOTH = TOTAL_STEPS * SUBSTEPS   # 72 000 Punkte / Arm


# ======================================================================
# Koordinatensystem D-Pol-Frame -> kartesisch
# ======================================================================
_l = math.radians(DP_L)
_b = math.radians(DP_B)
DP_CART = np.array([
    math.cos(_b) * math.cos(_l),
    math.cos(_b) * math.sin(_l),
    math.sin(_b),
], dtype=np.float64)


def _build_rot(target):
    z  = np.array([0., 0., 1.])
    dp = target / np.linalg.norm(target)
    ax = np.cross(z, dp)
    n  = np.linalg.norm(ax)
    if n < 1e-10:
        return np.eye(3) * (1.0 if np.dot(z, dp) > 0 else -1.0)
    ax /= n
    ca = float(np.dot(z, dp))
    sa = float(n)
    K  = np.array([[0, -ax[2], ax[1]],
                   [ax[2], 0, -ax[0]],
                   [-ax[1], ax[0], 0]])
    return np.eye(3) + sa * K + (1 - ca) * (K @ K)


ROT = _build_rot(DP_CART)


def dpol_to_unit(theta_deg, phi_deg):
    """Gibt float64 Einheitsvektor (D-Pol-Frame -> kartesisch)."""
    t = math.radians(theta_deg)
    p = math.radians(phi_deg)
    v = np.array([math.sin(t) * math.cos(p),
                  math.sin(t) * math.sin(p),
                  math.cos(t)], dtype=np.float64)
    return ROT @ v


# ======================================================================
# Schalenwinkel theta_n fuer n=1..6
# ======================================================================
def shell_thetas(chi):
    return [
        math.degrees(math.acos(max(-1.0, min(1.0,
            math.cos(n * math.radians(DELTA)) * math.cos(n * math.radians(chi))
        ))))
        for n in range(1, 7)
    ]


# ======================================================================
# Slerp — Grosskreis-Interpolation mit linearer Radiusinterpolation
# ======================================================================
def _slerp_seg(p0, p1, substeps):
    """
    p0, p1: 3D-Vektoren (moeglicherweise mit Expansionsradius, nicht unit).
    Gibt `substeps` Punkte zurueck (t=0..S-1/S), Richtung per Slerp,
    Radius linear interpoliert.
    """
    r0 = np.linalg.norm(p0); r1 = np.linalg.norm(p1)
    u0 = p0 / max(r0, 1e-15)
    u1 = p1 / max(r1, 1e-15)

    dot = float(np.dot(u0, u1))
    dot = max(-1.0, min(1.0, dot))

    t_arr = np.arange(substeps, dtype=np.float64) / substeps
    r_arr = r0 + t_arr * (r1 - r0)   # Radius linear interpolieren

    if abs(dot) > 0.9999:
        dirs = u0[None, :] + t_arr[:, None] * (u1 - u0)
        nrm  = np.linalg.norm(dirs, axis=1, keepdims=True).clip(1e-15)
        return (dirs / nrm * r_arr[:, None]).astype(np.float32)

    omega = math.acos(dot)
    s     = math.sin(omega)
    w0    = (np.sin((1.0 - t_arr) * omega) / s)[:, None]
    w1    = (np.sin(t_arr  * omega) / s)[:, None]
    dirs  = w0 * u0 + w1 * u1
    return (dirs * r_arr[:, None]).astype(np.float32)


# ======================================================================
# 2-Arm Trajektorienberechnung (Kreuzversatz + Expansion)
# ======================================================================
def compute_arm(arm_idx, chi):
    """
    Erzeugt den vollen Kreisbahn-Pfad fuer einen der 2 Arme.

    arm_idx=0 (Arm A, prograd):
      near-side base = BASE_ROT = 88 deg
      far-side  base = BASE_ROT + 180 = 268 deg
    arm_idx=1 (Arm B, retrograd):
      near-side base = BASE_ROT + chi/3 = 107.7 deg
      far-side  base = BASE_ROT + chi/3 + 180 = 287.7 deg

    An jedem Kehrtwende-Punkt (end of cycle, theta=8.1°, theta nahe D-Pol)
    wechselt der Arm zur ANDEREN Seite der Kugel (Kreuzversatz).
    Der Kreuzungspunkt liegt bei theta≈8° (D-Pol-Naehe), NICHT im Zentrum!

    Expansion: Schritt n steht auf Sphere mit R = (1+EPS_EXP)^(n//6).
    """
    thetas  = shell_thetas(chi)
    tilt    = chi / 3.0
    drift_k = abs(chi - 60.0)

    near_base = (BASE_ROT + arm_idx * tilt) % 360.0
    far_base  = (near_base + 180.0)        % 360.0

    disc   = np.empty((TOTAL_STEPS, 3), dtype=np.float64)
    az_acc = 0.0

    for n in range(1, TOTAL_STEPS + 1):
        cycle = (n - 1) // 6
        sic   = ((n - 1) % 6) + 1
        dirn  = +1 if (cycle % 2 == 0) else -1
        az_acc += dirn * tilt

        # Kreuzversatz: near-side bei geraden Zyklen, far-side bei ungeraden
        base  = near_base if (cycle % 2 == 0) else far_base
        az    = base + az_acc + cycle * drift_k
        theta = thetas[sic - 1]

        # Expansion: Radius waechst pro Zyklus
        r = (1.0 + EPS_EXP) ** cycle
        disc[n - 1] = r * dpol_to_unit(theta, az)

    # Glatten Pfad per Slerp
    smooth = np.empty((TOTAL_SMOOTH, 3), dtype=np.float32)
    for i in range(TOTAL_STEPS - 1):
        smooth[i * SUBSTEPS: (i + 1) * SUBSTEPS] = \
            _slerp_seg(disc[i], disc[i + 1], SUBSTEPS)
    smooth[(TOTAL_STEPS - 1) * SUBSTEPS:] = disc[-1].astype(np.float32)

    return smooth


def compute_all(chi):
    tilt = chi / 3.0
    print(f"  2-Arm Kreuzversatz: chi={chi:.2f}  tilt={tilt:.3f}  "
          f"near=[{BASE_ROT:.1f}, {BASE_ROT+tilt:.1f}]  "
          f"far=[{BASE_ROT+180:.1f}, {BASE_ROT+180+tilt:.1f}]  "
          f"drift/Zyklus={abs(chi-60):.2f} deg ...", flush=True)
    result = [compute_arm(i, chi) for i in range(2)]
    print("  Fertig.", flush=True)
    return result


# ======================================================================
# Initialer Build
# ======================================================================
print("METAGEOMETRA V21.3 — 2-Arm Kreuzversatz + Expansion (v3)")
SMOOTH_PATHS = compute_all(CHI)

SMBH_POS = [
    (nm, dpol_to_unit(td, ad).astype(np.float32), col, ref, arm_id, sn)
    for nm, td, ad, col, ref, arm_id, sn in CONFIRMED_SMBHS
]


# ======================================================================
# OpenGL-Zeichenfunktionen
# ======================================================================

def draw_sphere_fill(r=1.0, alpha=0.08):
    glColor4f(0.04, 0.06, 0.18, alpha)
    glDepthMask(GL_FALSE)
    q = gluNewQuadric()
    gluQuadricDrawStyle(q, GLU_FILL)
    gluSphere(q, r, 40, 20)
    glDepthMask(GL_TRUE)


def draw_sphere_wire(r=1.0):
    glColor4f(0.12, 0.22, 0.50, 0.14)
    glLineWidth(0.6)
    q = gluNewQuadric()
    gluQuadricDrawStyle(q, GLU_LINE)
    gluSphere(q, r, 36, 18)


def draw_shell_rings(chi, r=1.0):
    """
    Schalenringe mit Kurz-Kurz-Lang-Kurz Farbkodierung.
    Doublett-Paare (n=1/5 und n=2/4) durch weisse Verbindungslinie hervorgehoben.
    n=3 (LANG, anti-D-Pol) leuchtet orange.
    """
    thetas = shell_thetas(chi)

    # Farben: n=1..6, gelb/orange für LANG (n=3), cyan für Kehrtwende (n=6)
    ring_styles = [
        (0.90, 0.40, 0.40, 0.20, 1.2),   # n=1 KURZ – rot
        (0.40, 0.90, 0.40, 0.20, 1.2),   # n=2 KURZ – gruen
        (1.00, 0.70, 0.10, 0.35, 2.5),   # n=3 LANG – orange (dicker!)
        (0.40, 0.90, 0.40, 0.12, 0.8),   # n=4 KURZ (Doublett n=2) – gedimmt
        (0.90, 0.40, 0.40, 0.12, 0.8),   # n=5 KURZ (Doublett n=1) – gedimmt
        (0.20, 0.90, 0.95, 0.30, 2.0),   # n=6 Kehrtwende – cyan
    ]

    pts = np.empty((80, 3), dtype=np.float32)
    glEnableClientState(GL_VERTEX_ARRAY)
    for i, (theta, (rc, gc, bc, alpha, lw)) in enumerate(zip(thetas, ring_styles)):
        for j in range(80):
            pts[j] = (r * dpol_to_unit(theta, 360.0 * j / 80)).astype(np.float32)
        glColor4f(rc, gc, bc, alpha)
        glLineWidth(lw)
        glVertexPointer(3, GL_FLOAT, 0, pts)
        glDrawArrays(GL_LINE_LOOP, 0, 80)
    glDisableClientState(GL_VERTEX_ARRAY)

    # Doublett-Brücken: verbinde n=1/n=5 und n=2/n=4 sichtbar
    glLineWidth(1.4)
    for (i0, i1), az_ref in [((0, 4), 109.0), ((1, 3), 270.0)]:
        glColor4f(1.0, 1.0, 1.0, 0.35)
        glBegin(GL_LINES)
        glVertex3fv((r * dpol_to_unit(thetas[i0], az_ref)).astype(np.float32))
        glVertex3fv((r * dpol_to_unit(thetas[i1], az_ref)).astype(np.float32))
        glEnd()


def draw_dpol_axis(r=1.0):
    dp = (DP_CART / np.linalg.norm(DP_CART)).astype(np.float32)
    tip = dp * r * 1.52
    bot = dp * r * -1.42

    glLineWidth(2.5)
    glBegin(GL_LINES)
    glColor4f(0.0, 1.0, 1.0, 0.95)
    glVertex3fv(tip)
    glColor4f(0.0, 0.4, 0.6, 0.18)
    glVertex3fv(bot)
    glEnd()

    glPointSize(16.0)
    glColor4f(0.0, 1.0, 1.0, 1.0)
    glBegin(GL_POINTS)
    glVertex3fv(tip)
    glEnd()

    # Ring am Pol
    prp = np.cross(dp, np.array([0., 1., 0.], dtype=np.float32))
    if np.linalg.norm(prp) < 0.05:
        prp = np.cross(dp, np.array([1., 0., 0.], dtype=np.float32))
    prp  = (prp / np.linalg.norm(prp)).astype(np.float32)
    prp2 = np.cross(dp, prp).astype(np.float32)
    prp2 = (prp2 / np.linalg.norm(prp2)).astype(np.float32)

    pts = np.empty((80, 3), dtype=np.float32)
    glEnableClientState(GL_VERTEX_ARRAY)
    for rr, alpha, lw in [(0.09, 0.85, 2.2), (0.15, 0.22, 5.0)]:
        rscl = rr * r
        for i in range(80):
            a  = 2 * math.pi * i / 80
            p  = tip + math.cos(a) * rscl * prp + math.sin(a) * rscl * prp2
            n_ = np.linalg.norm(p)
            pts[i] = (p / n_ * r).astype(np.float32) if n_ > 1e-10 else tip
        glLineWidth(lw)
        glColor4f(0.0, 1.0, 1.0, alpha)
        glVertexPointer(3, GL_FLOAT, 0, pts)
        glDrawArrays(GL_LINE_LOOP, 0, 80)
    glDisableClientState(GL_VERTEX_ARRAY)


def draw_arm(smooth, n_pts, color, arm_active):
    """
    Zeichnet den Arm als EINE durchgehende Linie (GL_LINE_STRIP).
    Koordinaten enthalten bereits die Expansionsskalierung.
    Aeltere Abschnitte: gedimmt. Aktuelle: hell.
    """
    if not arm_active or n_pts < 2:
        return

    r, g, b = color
    n = min(n_pts, TOTAL_SMOOTH)
    RECENT = 600
    split  = max(0, n - RECENT)

    glEnableClientState(GL_VERTEX_ARRAY)

    if split >= 2:
        glColor4f(r * 0.38, g * 0.38, b * 0.38, 0.52)
        glLineWidth(1.1)
        glVertexPointer(3, GL_FLOAT, 0, smooth)
        glDrawArrays(GL_LINE_STRIP, 0, split + 1)

    if n - split >= 2:
        glColor4f(r, g, b, 0.95)
        glLineWidth(2.8)
        glVertexPointer(3, GL_FLOAT, 0, smooth[split:])
        glDrawArrays(GL_LINE_STRIP, 0, n - split)

    glDisableClientState(GL_VERTEX_ARRAY)

    # Kopf-Glow
    p = smooth[n - 1]
    glPointSize(32.0); glColor4f(r, g, b, 0.25)
    glBegin(GL_POINTS); glVertex3fv(p); glEnd()
    glPointSize(12.0); glColor4f(r, g, b, 0.95)
    glBegin(GL_POINTS); glVertex3fv(p); glEnd()
    glPointSize(4.5);  glColor4f(1.0, 1.0, 1.0, 1.0)
    glBegin(GL_POINTS); glVertex3fv(p); glEnd()


def draw_smbhs(time_sec, r=1.0):
    """
    SMBHs immer sichtbar (Depth-Test deaktiviert), pulsierend.
    Skaliert mit aktuellem Expansionsradius r.
    """
    glDisable(GL_DEPTH_TEST)
    pulse = 0.60 + 0.40 * math.sin(time_sec * 2.8)

    for nm, unit_pos, col, ref, arm_id, sn in SMBH_POS:
        cr, cg, cb = col
        p  = (r * unit_pos).astype(np.float32)
        dp = (p / max(float(np.linalg.norm(p)), 1e-15)).astype(np.float32)

        prp = np.cross(dp, np.array([0., 1., 0.], dtype=np.float32))
        if np.linalg.norm(prp) < 0.05:
            prp = np.cross(dp, np.array([1., 0., 0.], dtype=np.float32))
        prp  = (prp / np.linalg.norm(prp)).astype(np.float32)
        prp2 = np.cross(dp, prp).astype(np.float32)
        prp2 = (prp2 / np.linalg.norm(prp2)).astype(np.float32)

        pts = np.empty((60, 3), dtype=np.float32)
        glEnableClientState(GL_VERTEX_ARRAY)
        for rr, alpha, lw in [(0.055 * r, 0.45 * pulse, 2.5),
                              (0.095 * r, 0.20 * pulse, 5.0)]:
            for i in range(60):
                a  = 2 * math.pi * i / 60
                pt = p + math.cos(a) * rr * prp + math.sin(a) * rr * prp2
                n_ = float(np.linalg.norm(pt))
                pts[i] = (pt / n_ * r).astype(np.float32) if n_ > 1e-10 else p
            glColor4f(cr, cg, cb, alpha)
            glLineWidth(lw)
            glVertexPointer(3, GL_FLOAT, 0, pts)
            glDrawArrays(GL_LINE_LOOP, 0, 60)
        glDisableClientState(GL_VERTEX_ARRAY)

        d = 0.072 * r
        glLineWidth(2.2)
        glColor4f(cr, cg, cb, 0.92)
        for ax_v in [np.array([1., 0., 0.]), np.array([0., 1., 0.]), np.array([0., 0., 1.])]:
            glBegin(GL_LINES)
            glVertex3fv((p - ax_v * d).astype(np.float32))
            glVertex3fv((p + ax_v * d).astype(np.float32))
            glEnd()

        glPointSize(22.0); glColor4f(cr, cg, cb, pulse)
        glBegin(GL_POINTS); glVertex3fv(p); glEnd()
        glPointSize(8.0);  glColor4f(1., 1., 1., 1.0)
        glBegin(GL_POINTS); glVertex3fv(p); glEnd()

    glEnable(GL_DEPTH_TEST)


# ======================================================================
# 2D HUD
# ======================================================================
def _init_fonts():
    pygame.font.init()
    try:
        fs = pygame.font.SysFont("Consolas,Courier New,Lucida Console", 13)
        fm = pygame.font.SysFont("Consolas,Courier New,Lucida Console", 14, bold=True)
        fl = pygame.font.SysFont("Consolas,Courier New,Lucida Console", 18, bold=True)
    except Exception:
        fs = pygame.font.Font(None, 17)
        fm = pygame.font.Font(None, 20)
        fl = pygame.font.Font(None, 25)
    return fs, fm, fl


def _blit(font, text, x, y, color, H):
    surf = font.render(text, True, color)
    w, h = surf.get_size()
    if x + w > 0 and y + h > 0:
        glRasterPos2i(x, H - y - h)
        glDrawPixels(w, h, GL_RGBA, GL_UNSIGNED_BYTE,
                     pygame.image.tostring(surf, "RGBA", True))


def hud_begin(W, H):
    glMatrixMode(GL_PROJECTION)
    glPushMatrix(); glLoadIdentity()
    gluOrtho2D(0, W, 0, H)
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix(); glLoadIdentity()
    glDisable(GL_DEPTH_TEST)


def hud_end():
    glEnable(GL_DEPTH_TEST)
    glMatrixMode(GL_MODELVIEW); glPopMatrix()
    glMatrixMode(GL_PROJECTION); glPopMatrix()
    glMatrixMode(GL_MODELVIEW)


def draw_hud(W, H, fs, fm, fl, chi, n_pts, paused, arm_vis,
             anim_speed, r_now, show_expand):
    tilt   = chi / 3.0
    drift  = abs(chi - 60.0)
    thetas = shell_thetas(chi)

    step_n   = min(n_pts // SUBSTEPS, TOTAL_STEPS - 1)
    cycle    = step_n // 6
    sic      = (step_n % 6) + 1
    side     = "near-side" if (cycle % 2 == 0) else "far-side (Kreuzversatz!)"
    theta_n  = thetas[sic - 1]
    rhythm   = SHELL_RHYTHM[sic - 1]

    hud_begin(W, H)

    def sm(t, x, y, c=(185, 215, 255)): _blit(fs, t, x, y, c, H)
    def md(t, x, y, c=(120, 200, 255)): _blit(fm, t, x, y, c, H)
    def lg(t, x, y, c=(90, 195, 255)):  _blit(fl, t, x, y, c, H)

    # Titel
    lg("METAGEOMETRA  V21.3", 14, 12, (90, 195, 255))
    md("2-Arm Kreuzversatz  |  Kehrtwende am D-Pol  |  Raum-Expansion  (v3)", 14, 34,
       (70, 160, 240))

    # Parameter
    sm(f"chi={chi:.2f}°   chi/3={tilt:.3f}°   Drift/Zyk={drift:.3f}°   "
       f"eta={ETA:.2e}   eps_exp={EPS_EXP:.2e}",
       14, 57)

    # Chi + Expansion Info rechts
    c_chi = (90, 255, 90) if abs(chi - 59.1) < 0.05 else (255, 210, 60)
    sm("[ ] = chi +/-0.1   E = Expansion", W - 240, 12, c_chi)
    r_col = (255, 180, 80) if show_expand else (100, 100, 100)
    sm(f"R(Zyklus {cycle}) = {r_now:.6f}   {'aktiv' if show_expand else 'AUS'}",
       W - 240, 27, r_col)

    # Schritt-Info mit Rhythmus
    rhy_col = (255, 100, 100) if rhythm == "LANG" else (100, 255, 100)
    md(f"Schritt {step_n:5d}/{TOTAL_STEPS}   Zyklus {cycle:4d}"
       f"   n={sic}   theta={theta_n:.1f}°   [{rhythm}]   {side}",
       14, 77, (255, 220, 70))

    sm(f"Az-Drift kumuliert: {cycle * drift:.1f}°   "
       f"Az near=[{BASE_ROT:.1f}°, {BASE_ROT+tilt:.1f}°]   "
       f"far=[{BASE_ROT+180:.1f}°, {BASE_ROT+180+tilt:.1f}°]",
       14, 96)

    # Schalen mit Rhythmus-Annotation
    tstr_parts = []
    for i, t in enumerate(thetas):
        rhy = SHELL_RHYTHM[i]
        marker = "*" if i in (0, 4) or i in (1, 3) else "~" if i == 5 else "!"
        tstr_parts.append(f"n{i+1}={t:.1f}°[{rhy}]{marker}")
    sm("Schalen:  " + "  ".join(tstr_parts), 14, 112, (150, 170, 210))
    sm("Doublett: n1/n5 (59°/64°) KURZ-Paar ↔ n2/n4 (118°/123°) KURZ-Paar   n3 (176°) LANG   n6 (8°) Kehrtwende",
       14, 127, (180, 180, 130))

    # Arm-Legende
    y = 146
    sm("Arme  (1/2 ein/aus):", 14, y, (200, 200, 200))
    y += 15
    for i, (nm, col) in enumerate(zip(ARM_NAMES, ARM_COLORS)):
        rc = (int(col[0]*255), int(col[1]*255), int(col[2]*255))
        c  = rc if arm_vis[i] else (70, 70, 70)
        ex = "" if arm_vis[i] else "  [AUS]"
        sm(f"  {i+1}. {nm}{ex}", 14, y, c)
        y += 14

    # SMBH-Legende
    y += 5
    sm("Bestaetigte SMBHs   (Tiefe: Depth-Test AUS, immer sichtbar):", 14, y, (200, 200, 200))
    y += 14
    for nm, unit_pos, col, ref, arm_id, sn in SMBH_POS:
        td_v = [s[1] for s in CONFIRMED_SMBHS if s[0] == nm][0]
        az_v = [s[2] for s in CONFIRMED_SMBHS if s[0] == nm][0]
        rhy  = SHELL_RHYTHM[sn - 1]
        rc   = (int(col[0]*255), int(col[1]*255), int(col[2]*255))
        sm(f"  o  {nm:10s}  n={sn}[{rhy}]  theta={td_v:.1f}°  Az={az_v}°  Arm {arm_id}  [{ref}]",
           14, y, rc)
        y += 13

    # D-Pol
    y += 3
    sm(f"D-Pol: l={DP_L}°  b={DP_B}°  |  KREUZUNGSPUNKT = theta≈8° (D-Pol-Naehe), NICHT Kugelzentrum!",
       14, y, (70, 215, 215))

    # Steuerung
    sm("SPACE=Pause  Arr=Speed  R=Reset  G=Gitter  S=Sph  D=DPol  "
       "E=Expansion  [ ]=chi  1/2=Arme  ESC=Ende",
       14, H - 20, (105, 105, 160))

    if paused:
        md("PAUSE", W // 2 - 40, H // 2 - 10, (255, 200, 0))

    # Fortschrittsbalken
    bx, by, bw, bh = 14, H - 36, W - 170, 5
    prog = n_pts / TOTAL_SMOOTH
    glColor4f(0.12, 0.22, 0.45, 0.55)
    glBegin(GL_QUADS)
    glVertex2f(bx,          by);      glVertex2f(bx + bw,         by)
    glVertex2f(bx + bw,     by + bh); glVertex2f(bx,              by + bh)
    glEnd()
    glColor4f(0.28, 0.68, 1.00, 0.90)
    glBegin(GL_QUADS)
    glVertex2f(bx,              by);      glVertex2f(bx + bw * prog,  by)
    glVertex2f(bx + bw * prog,  by + bh); glVertex2f(bx,              by + bh)
    glEnd()

    hud_end()


# ======================================================================
# Hauptschleife
# ======================================================================
def main():
    global SMOOTH_PATHS, CHI, SMBH_POS

    pygame.init()
    fs, fm, fl = _init_fonts()

    W, H = 1340, 820
    pygame.display.set_mode((W, H), DOUBLEBUF | OPENGL)
    pygame.display.set_caption(
        f"METAGEOMETRA V21.3  |  chi={CHI:.1f}°  |  2-Arm Kreuzversatz + Expansion")

    glViewport(0, 0, W, H)
    glMatrixMode(GL_PROJECTION); glLoadIdentity()
    gluPerspective(42.0, W / H, 0.01, 200.0)
    glMatrixMode(GL_MODELVIEW)

    glEnable(GL_DEPTH_TEST);  glDepthFunc(GL_LEQUAL)
    glEnable(GL_BLEND);       glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glEnable(GL_LINE_SMOOTH); glEnable(GL_POINT_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST)

    # Zustand
    rot_x, rot_y = 18.0, -22.0
    dragging     = False
    last_pos     = (0, 0)
    zoom         = 3.6

    smooth_n    = 0.0
    anim_speed  = 3.5
    paused      = False
    show_grid   = True
    show_fill   = True
    show_dpol   = True
    show_expand = True
    arm_vis     = [True, True]

    clock = pygame.time.Clock()
    t0    = pygame.time.get_ticks()

    while True:
        clock.tick(60)
        time_sec = (pygame.time.get_ticks() - t0) / 1000.0

        # Events
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit(); sys.exit()

            elif ev.type == pygame.KEYDOWN:
                k = ev.key
                if   k == pygame.K_ESCAPE:   pygame.quit(); sys.exit()
                elif k == pygame.K_SPACE:    paused = not paused
                elif k == pygame.K_UP:       anim_speed = min(60.0, anim_speed * 1.5)
                elif k == pygame.K_DOWN:     anim_speed = max(0.05, anim_speed / 1.5)
                elif k == pygame.K_r:        smooth_n = 0.0
                elif k == pygame.K_g:        show_grid = not show_grid
                elif k == pygame.K_s:        show_fill = not show_fill
                elif k == pygame.K_d:        show_dpol = not show_dpol
                elif k == pygame.K_e:        show_expand = not show_expand
                elif k == pygame.K_1:        arm_vis[0] = not arm_vis[0]
                elif k == pygame.K_2:        arm_vis[1] = not arm_vis[1]
                elif k == pygame.K_LEFTBRACKET:
                    CHI = round(max(10.0, CHI - 0.1), 2)
                    SMOOTH_PATHS = compute_all(CHI)
                    SMBH_POS = [
                        (nm, dpol_to_unit(td, ad).astype(np.float32), col, ref, aid, sn)
                        for nm, td, ad, col, ref, aid, sn in CONFIRMED_SMBHS
                    ]
                    smooth_n = 0.0
                    pygame.display.set_caption(
                        f"METAGEOMETRA V21.3  |  chi={CHI:.2f}°")
                elif k == pygame.K_RIGHTBRACKET:
                    CHI = round(min(90.0, CHI + 0.1), 2)
                    SMOOTH_PATHS = compute_all(CHI)
                    SMBH_POS = [
                        (nm, dpol_to_unit(td, ad).astype(np.float32), col, ref, aid, sn)
                        for nm, td, ad, col, ref, aid, sn in CONFIRMED_SMBHS
                    ]
                    smooth_n = 0.0
                    pygame.display.set_caption(
                        f"METAGEOMETRA V21.3  |  chi={CHI:.2f}°")

            elif ev.type == pygame.MOUSEBUTTONDOWN:
                if   ev.button == 1: dragging = True;  last_pos = ev.pos
                elif ev.button == 4: zoom = max(1.5, zoom - 0.18)
                elif ev.button == 5: zoom = min(8.0, zoom + 0.18)

            elif ev.type == pygame.MOUSEBUTTONUP and ev.button == 1:
                dragging = False

            elif ev.type == pygame.MOUSEMOTION and dragging:
                dx = ev.pos[0] - last_pos[0]
                dy = ev.pos[1] - last_pos[1]
                rot_y   += dx * 0.36
                rot_x   += dy * 0.36
                rot_x    = max(-90.0, min(90.0, rot_x))
                last_pos = ev.pos

        # Animation
        if not paused:
            smooth_n += anim_speed
            if smooth_n >= TOTAL_SMOOTH:
                smooth_n = 0.0

        n_pts = max(2, int(smooth_n))

        # Aktueller Expansionsradius
        step_now  = min(n_pts // SUBSTEPS, TOTAL_STEPS - 1)
        cycle_now = step_now // 6
        r_now = (1.0 + EPS_EXP) ** cycle_now if show_expand else 1.0

        # 3D Rendering
        glClearColor(0.007, 0.007, 0.045, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glMatrixMode(GL_MODELVIEW); glLoadIdentity()
        glTranslatef(0.0, 0.0, -zoom)
        glRotatef(rot_x, 1.0, 0.0, 0.0)
        glRotatef(rot_y, 0.0, 0.0, 1.0)

        if show_fill:
            draw_sphere_fill(r=r_now)
        if show_grid:
            draw_shell_rings(CHI, r=r_now)
            draw_sphere_wire(r=r_now)
        if show_dpol:
            draw_dpol_axis(r=r_now)

        for i in range(2):
            draw_arm(SMOOTH_PATHS[i], n_pts, ARM_COLORS[i], arm_vis[i])

        draw_smbhs(time_sec, r=r_now)

        # HUD
        draw_hud(W, H, fs, fm, fl, CHI, n_pts, paused, arm_vis,
                 anim_speed, r_now, show_expand)

        pygame.display.flip()


if __name__ == "__main__":
    main()
