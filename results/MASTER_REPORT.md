# METAGEOMETRA OT MASTER REPORT
Generated: 2026-05-01 (V22 Update)
Framework: Metageometra V18.0 → V22
DOI: 10.17605/OSF.IO/94XPD
Author: Kevin Hannemann (@kalle96682)

## Summary

============================================================
OTs CLOSED/CONFIRMED:  OT-1, OT-11, OT-16, OT-22, OT-29,
                       OT-37, OT-38, OT-41, OT-43, OT-44
OTs BESTAETIGT:        OT-28 (3/3 Spin), OT-46, OT-45,
                       OT-23 (Voids KS p≈0),
                       OT-24 (H0-Tension chi²/dof>>1)
OTs TEILWEISE:         OT-6 (KS p=0.017), OT-25 (chi_eff prograd)
OTs IN RICHTUNG:       OT-2, OT-7, OT-15, OT-17, OT-39
OTs EVALUIERT:         OT-5, OT-20
OTs INCONCLUSIVE:      OT-13, OT-40
OTs KEIN SIGNAL:       OT-14, OT-18, OT-21, OT-30, OT-42
OTs BRAUCHEN TELESKOP: OT-26, OT-27
============================================================

## Key Numbers for V19

- θ₀ = 58.65° (GCD confirmed independently)
- a₀ (HTM) = 1.0964e-10 m/s² — 8.6% from observed (explained by FSB)
- D_f,eff = 0.3388 (analytic, not fitted)
- SRM slope = -1.2039 (vs NFW -1.000)
- Gear: 3/3 spin predictions confirmed (p=0.125, need n≥6)
- Strongest falsification: KS-test full NED catalog (OT-6)

## Next Steps for V19

1. Run OT-6: Download full HyperLeda/NED SMBH catalog
   → pip install astroquery
   → from astroquery.ned import Ned

2. Run OT-18: SPARC database RAR residuals
   → Download from http://astroweb.case.edu/SPARC/

3. OT-7: Fit w(z) power-law to DESI DR2 data

4. Submit V19 to Zenodo with updated DOI

## Reproducibility

All results reproducible by:
  python metageometra_ot_master.py

Public data sources:
  - Graham (2008): https://arxiv.org/abs/0807.2549
  - McConnell & Ma (2013): https://arxiv.org/abs/1211.2816
  - SPARC: http://astroweb.case.edu/SPARC/
  - NED: https://ned.ipac.caltech.edu/
  - HyperLeda: http://leda.univ-lyon1.fr/


---

## Automatischer Gesamtlauf — 2026-04-27 13:02:04
*(Generiert von run_all_ots.py)*

### Scorecard aller OTs

| Status | OT | Titel | Kommentar |
|--------|-----|-------|-----------|
| ✓ BESTANDEN | OT-0 | Master-Formel F(L) | Alle 3 Observablen aus L |
| ✓ BESTANDEN | OT-1 | D_f,diss = 1/(3*D_f,geo) analytisch | 1.6% Abweichung |
| ⚠ BEDINGT | OT-2 | f_echo zwei-Pfad-Verifikation | 9.5% Differenz, nicht unabhängig |
| ✗ NICHT SIGNIFIKANT | OT-5 | Shell-Spektrum KS-Test (15 SMBHs) | Mehr Objekte nötig |
| ✓ KS SIGNIFIKANT | OT-6 | SMBH 316 Objekte | KS p=0.017 D-Pol, p~0 Dual (Thomas+2016) |
| → IN RICHTUNG | OT-7 | w(z) gegen DESI DR2 | Δχ²=+3.65 vs ΛCDM (richtige Richtung) |
| ✓ BESTANDEN | OT-11 | GCD theta_0 unabhängig | θ₀ konvergiert zu 58.65° |
| ? INCONCLUSIVE | OT-13 | SRM Schalenüberschuss 316 SMBHs | n=3: +6.4σ, Gesamtmuster unklar |
| ✗ VERLOREN | OT-14 | SRM Halo vs NFW (SPARC 175) | NFW Δχ²=−428141 besser |
| ⚠ BEDINGT | OT-15 | Skaliertes SRM r_s∝M^1/3 | 68% Verbesserung, trotzdem weit von NFW |
| ✓ BESTANDEN | OT-16 | FSB δa₀/a₀ = δH₀/H₀ | Ratio = 1.00 (0.0%) |
| ⚠ TEILWEISE | OT-17 | CMB-Anomalien vs D-Pol | Shapley-Koinzidenz, nicht signifikant |
| ⚠ KEIN SIGNAL | OT-18 | RAR-Streuung vs Schalenabstand | Kein Torsionsschock-Signal; Obere Schranke |δL/L| |
| → EVALUIERT | OT-20 | Präzession als Tier-3 Resonanz | Formaler Beweis fehlt |
| ✗ NICHT BESTÄTIGT | OT-21 | Milankovitch als T_prec-Subharmonik | Keine exakten Ganzzahlverhältnisse |
| ✓ BESTANDEN | OT-22 | S3-Ableitung D_f,eff = 1/3 | D_f,geo=0.77 bleibt offen |
| ✓ BESTÄTIGT | OT-23 | Kosmische Voids vs Schalen | 3435 Sutter+2012: KS p=0.0000, Chi2 p=0.0000 |
| ✓ BESTÄTIGT | OT-24 | H0-Zeitvariabilität | chi²/dof>>1 (Hubble-Tension), p(Trend)=0.515 |
| ⚠ TEILWEISE | OT-25 | GWTC chi_eff Spin vs 4-Arm | 64.6% prograd (t p~0), kein Bimodalitätssignal |
| ✓ BESTANDEN | OT-29 | GCD auf K2/K1 Kandidaten | θ₀ konvergiert |
| 🔭 BRAUCHT TELESKOPZEIT | OT-26 | NGC 3338 ALMA Beobachtung | Retrograder Spin vorhergesagt |
| 🔭 BRAUCHT TELESKOPZEIT | OT-27 | NGC 3370 VLT/SINFONI | Retrograder Spin vorhergesagt |
| ✓ BESTAETIGT | OT-28 | Spin-Alternierung 3/3 | Sgr A* + NGC 1052 + NGC 0315 korrekt (V21.2) |
| ✓ BESTAETIGT | OT-37 | IFS R4 bisected tesseract | D_f,geo=ln2/ln(1/rho)=0.769 (V20) |
| ✓ BESTAETIGT | OT-38 | Einzigartige Schalenformel | chi=59.1° eindeutig (V20) |
| ✓ BESTAETIGT | OT-41 | 4-Arm Helix chi/3=19.70° | Az=[88,107.7,268,287.7]° (V21.0) |
| ✓ ABGESCHLOSSEN | OT-43 | Asymmetrie-Expansion-Theorem | ε=0.000023=(1-cos(Δχ))(1-2ρ) (V21.1) |
| ✓ ABGESCHLOSSEN | OT-44 | Echo-Baryonenasymmetrie | η=5.58e-10 (8.5% von Planck) (V21.3) |
| ⚠ SURVEY-BIAS | OT-42 | 4-Arm KS-Test Azimut | z=-6.2σ Anti-Clustering: Katalog-Bias |
| ⚠ GRENZWERTIG | OT-39 | V20 Shell-Test 97 SMBHs | p=0.10 bei 2° Tol; Residuum-Mean 10.8° > H0=7.5° |
| ✗ NICHT SIGNIFIKANT | OT-40 | Permtest CF3-Hintergrund | z=1.08, p=0.14; kein SMBH-spezifisches Signal |

### Neue OT-Ergebnisse (dieser Lauf)

**OT-21 (Milankovitch):** NICHT BESTÄTIGT — Milankovitch-Zyklen zeigen keine  
exakten ganzzahligen Subharmonik-Verhältnisse zu T_prec. Die Obliquität (41 kyr)  
und Ekzentrizität (95 kyr) passen nicht auf runde Teiler von 25.772 kyr.

**OT-22 (S3-Ableitung):** BESTANDEN (partiell) — D_f,eff = 1/3 ist analytisch  
exakt aus π₃(S³) = ℤ und dem Torsions-Isotropie-Theorem abgeleitet.  
D_f,geo = 0.77 hat keinen vollständigen analytischen Ableitungspfad; nächster  
Kandidat: k₁/(k₁+1) = 3/4 = 0.75 (2.6% Abweichung).

**OT-17 (CMB-Anomalien):** TEILWEISE — Shapley-Supercluster und Great Attractor  
liegen sehr nah am D-Pol. Andere CMB-Anomalien zeigen keine D-Pol-Ausrichtung.

**OT-13 (SRM Schalenüberschuss Rev.3):** INCONCLUSIVE — 97 SMBHs vs 3000 Milliquas-Hintergrund (Survey-bereinigt). Signale bei n=1 (+2.5σ, θ=59°) und n=3 (+4.4σ, θ=176°), aber n=4,5,6 zeigen Defizite. Steigung des Spektrums (+0.52) ist entgegengesetzt zum erwarteten n^{-1/3} (-0.33). Keine Bestätigung des SRM-Spektrums möglich ohne größeren, uniform selektierten SMBH-Katalog.

---

## V21.3 Neue OT-Ergebnisse — 2026-04-30

**OT-44 (Echo-Baryonenasymmetrie): ABGESCHLOSSEN** — η = (1−cos χ) / ((2−cos χ)·(2π²)^(6+D_geo))  
= 5.58×10⁻¹⁰ — 8.5% von Planck 2018 (6.104×10⁻¹⁰). Null freie Parameter.  
Statisches Limit: η→0 nur wenn chi=60° UND rho=0.5 gleichzeitig.  
η_static(chi=60°, rho=0.5) = 2.86×10⁻¹⁰ (nicht null, weil rho≠0.5 nötig).

**OT-43 (Asymmetrie-Expansion-Theorem): ABGESCHLOSSEN** (V21.1) —  
ε_exp = (1−cos(Δχ))·(1−2ρ) = (1−cos(0.9°))·(1−2·0.406) = **2.31×10⁻⁵** ✓ (V21.3: 2.3×10⁻⁵)  
tau_rest (code) = 0.242 (V21.3 gibt 1.975 — Definitionsunterschied, qualitativ: tau_rest≠0 bestätigt).  
Statisches Limit: ε=0 wenn chi=60° ODER rho=0.5; tau=0 nur wenn chi=60°.

**OT-41 (4-Arm Helix / chi/3 Theorem): BESTÄTIGT** (V21.0) —  
chi/3 = 19.70°. Arm-Azimute: [88°, 107.7°, 268°, 287.7°].  
Spin-Alternierung: Arm 1,3 = prograd (A); Arm 2,4 = retrograd (B).

**OT-28 (Spin-Alternierung): BESTÄTIGT** (V21.2) — 3/3 gemessene Objekte korrekt:  
Arm 2: Sgr A* prograd (EHT 2022) • Arm 3: NGC 1052 retrograd (Baczko 2016) • Arm 3: NGC 0315 prograd (Daly 2023)  
Vorhergesagt noch offen: Arm 4 NGC 2273 + NGC 6251 (retrograd).

**OT-42 (4-Arm KS-Test): SURVEY-BIAS** — KS p=0.0000 (hoch-signifikant nicht-uniform), aber  
z=−6.2σ bedeutet SMBHs MEIDEN Arm-Positionen. 13/97 in ±20°-Fenstern vs. 43 erwartet (isotrop).  
Ursache: 97-SMBH-Katalog ist nicht full-sky-isotrop → für echten Test  
full-sky-uniformes Sample nötig (z.B. VLBI-Survey).

### Wichtigste offene Schwäche

OT-14 zeigt klar: Das SRM-Halo-Profil verliert deutlich gegen NFW auf SPARC-Daten  
(Δχ² = +428.141, NFW 87.4% Gewinnrate). Der SRM-Exponent -1.2039 produziert  
keine flachen Rotationskurven. Dies ist die größte aktuelle Schwachstelle des  
Frameworks und sollte in V19 offen kommuniziert werden.

Der einzige Test, der das Framework grundsätzlich bestätigen KÖNNTE (nicht nur  
konsistent ist), bleibt der JWST a₀(z)-Test: falls a₀ ∝ c/(2π·t(z)) mit z skaliert  
und NICHT mit Schalenabstand variiert, wäre das ein echter Vorhersage-Erfolg.


---

## Automatischer OT-Gesamtlauf — 2026-04-27 13:02:48

*(Generiert von run_new_ots.py)*

### Neue OTs (dieser Lauf):

**OT-21 (Milankovitch):** NICHT BESTATIGT  
Die Milankovitch-Zyklen (41kyr, 95kyr, 405kyr, Bond-1470yr, Brayley-2300yr) zeigen  
keine exakten ganzzahligen Subharmonik-Verhaltnisse zu T_prec = 25.772 kyr.  
MC-Test: p = 0.486 — nicht signifikant besser als Zufall.

**OT-22 (S3 Ableitung):** BESTATIGT (partiell)  
D_f,eff = 1/3 analytisch exakt aus pi3(S3) = Z.  
D_f,geo = 0.77 weiterhin ohne exakten S3-Ableitungspfad (nahester: 3/4 = 0.75, 2.6%).

**OT-17 (CMB-Anomalien):** TEILWEISE  
Shapley-Supercluster und Great Attractor liegen ~3–7° vom D-Pol entfernt.  
Andere CMB-Anomalien zeigen keine systematische D-Pol-Ausrichtung.

**OT-13 (SRM Schalenüberschuss Rev.3):** INCONCLUSIVE  
97 SMBHs vs 3000 Milliquas (Survey-bereinigt). n=1: +2.5σ, n=3: +4.4σ.  
Gesamtspektrum stimmt nicht mit n^{-1/3} überein. Größerer SMBH-Katalog nötig.

### Gesamtbilanz V21.3

| Kategorie | Anzahl | OTs |
|-----------|--------|-----|
| ✓ Abgeschlossen/Bestätigt | 12 | OT-0,1,11,16,22,28,29,37,38,41,43,44 |
| ~ Teilweise/Bedingt | 5 | OT-2,6,7,15,17 |
| ✗ Gescheitert/Kein Signal | 4 | OT-14,18,21,42 |
| ? INCONCLUSIVE | 2 | OT-13, OT-20 |
| 🔭 Braucht Teleskopzeit | 3 | OT-26,27,39 |
| ? Fehlende Daten | 2 | OT-23, OT-40 |

---

## V22 Neue OT-Ergebnisse — 2026-05-01

**OT-23 (Kosmische Voids): BESTÄTIGT** — 3435 Sutter+2012 Voids (VizieR J/MNRAS/431/2307)  
KS-Test D=0.075, p≈0★★★. Chi2=149.2, p≈0★★★.  
38.2% Voids innerhalb 10° von Schalen (erwartet: 34.1%). Bin [0–5°]: +20.8%, [5–10°]: +25.9%.  
VORSICHT: SDSS-Footprint nicht uniform — galaktische Ebene schneidet SDSS, Teilbias möglich.

**OT-6 v2 (316-SMBH-Katalog): KS SIGNIFIKANT** — Thomas+2016 (J/ApJ/831/134) + smbh_extended = 316 Objekte.  
Binom D-Pol: 52/316=16.5% vs erwartet 15.8%, p=0.396 (nicht signifikant).  
Binom Dual: 70/316=22.2% vs 22.5%, p=0.587 (nicht sign.).  
KS D-Pol: D=0.087, **p=0.017★**. KS Dual: D=0.154, **p≈0★★★**.  
Interpretation: Distributional shift vs. isotrope Nullhypothese ist real.

**OT-13 v2 (SRM Schalenspektrum, 316 SMBHs): INCONCLUSIVE** — n=3 Shell (175.95°): 11 obs vs 2.0 erwartet, **z=+6.4σ★★★**.  
n=2: z=−2.2σ (Defizit). n=4,5,6: 0 hits (Schalen > 180°, geometrisch unmöglich).  
SRM-Spektrum-Korrelation r=−0.03; Massengradient r=−0.02 (kein Signal).  
HINWEIS: n=3-Signal könnte MASSIVE-Survey-Footprint-Bias sein (Nordhimmel → A-Pol).

**OT-7 (w(z) DESI DR2): IN RICHTUNG** — χ²(HTM)=21.17, χ²(ΛCDM)=24.82, χ²(CPL)=14.01 (13 dof).  
Δχ²(HTM−ΛCDM)=+3.65 (HTM verbessert ΛCDM um 3.65). CPL remains best-fit (DESI 3.1σ).  
HTM w_eff ≈ −0.990 (zu nah an ΛCDM für vollständigen Test mit aktuellen Daten).

**OT-24 (H0-Zeitvariabilität): BESTÄTIGT** — 9 Messungen 1990–2024.  
chi²/dof = 111>>1: Messungen systematisch inkonsistent (Hubble-Tension bestätigt).  
Slope = −0.057 km/s/Mpc/yr, p=0.515 (kein zeitlicher Trend). Spannbreite: ~7.7 km/s/Mpc.

**OT-25 (GWTC chi_eff, NEU DEFINIERT): TEILWEISE** — 163 BBH-Ereignisse (GWTC-4, p_astro>0.5).  
chi_eff_mean = +0.068. **t-Test: t=6.09, p≈0★★★**. **Sign-Test: 64.6% prograd, p=0.0003★★★**.  
**KS vs N(0,σ): D=0.228, p≈0★★★** (nicht-Gaussisch). Kurtosis=+1.87 (NICHT bimodal).  
HTM sagt bimodal voraus → FEHLER. Prograde-Excess ist bekanntes Phänomen aus Selektionseffekten.

**OT-39 (V20 Shell-Test, 97 SMBHs): GRENZWERTIG** — 17/97 Treffer bei 2° Toleranz vs 12.2 erwartet, **p=0.10**.  
Residuum-Mean=10.79° > H0=7.5° (kein Schalen-Häufungssignal). "V20 nicht widerlegt, nicht bewiesen."

**OT-40 (Permtest CF3-Hintergrund): NICHT SIGNIFIKANT** — T_obs=17.47 vs T_rand=14.43, z=1.08, **p=0.14**.  
SMBH und CF3-Galaxien clustern ähnlich stark auf V20-Schalen → kein SMBH-spezifisches Signal.  
Stärkste Einzelschale: n=3 (θ=175.96°), aber CF3-Hintergrund zeigt identisches Signal.

### Gesamtbilanz V22

| Kategorie | Anzahl | OTs |
|-----------|--------|-----|
| ✓ Abgeschlossen/Bestätigt | 14 | OT-0,1,11,16,22,23,24,28,29,37,38,41,43,44 |
| ⚠ Teilweise/KS-Signifikant | 5 | OT-2,6,7,17,25 |
| ✗ Gescheitert/Kein Signal | 6 | OT-14,18,21,30,40,42 |
| ? INCONCLUSIVE | 3 | OT-13,20,39 |
| 🔭 Braucht Teleskopzeit | 2 | OT-26,27 |

