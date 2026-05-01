# METAGEOMETRA OT MASTER — COPILOT INSTRUCTIONS

## Deine einzige Aufgabe

Der User wirft Dateien rein. Du führst `metageometra_ot_master.py` aus.
Das war's. Keine Erklärungen. Kein Umschreiben. Keine Vorschläge.

---

## Was der User tun wird

1. Eine oder mehrere dieser Dateien einfügen:
   - Ein SMBH-Katalog (CSV, TXT, FITS) — NED / HyperLeda / Graham / McConnell
   - SPARC Rotationskurven (CSV)
   - DESI w(z) Daten (CSV, JSON)
   - Milliquas Katalog (TXT, FITS)
   - Irgendeine andere astronomische Datei

2. Sagen: "Führ es aus" / "Run it" / nichts sagen und nur die Datei droppen

---

## Was du tust

### Schritt 1 — Datei erkennen

| Dateiname enthält       | Datei ist                  | Zuständiger OT |
|-------------------------|----------------------------|----------------|
| graham, mcconnell, smbh | SMBH dynamischer Katalog   | OT-5, OT-6    |
| sparc, rotation         | Rotationskurven            | OT-14, OT-18  |
| milliquas, qso, agn     | Quasar/AGN Katalog         | OT-6          |
| desi, w_z, dark_energy  | Dark Energy Daten          | OT-7          |
| void, sdss              | Void Katalog               | OT-23         |
| leda, hyperleda         | HyperLeda Katalog          | OT-5, OT-6    |

### Schritt 2 — Katalog einlesen

```python
import pandas as pd
import numpy as np

# Versuche automatisch:
try:
    df = pd.read_csv(filepath, comment='#', sep=None, engine='python')
except:
    try:
        from astropy.io import fits
        hdul = fits.open(filepath)
        df = pd.DataFrame(hdul[1].data)
    except:
        df = pd.read_csv(filepath, delim_whitespace=True, comment='#')
```

### Schritt 3 — Koordinaten finden

Suche automatisch nach Spalten die so heißen könnten:
- RA: `ra`, `RA`, `RAJ2000`, `ra_deg`, `_RA`
- Dec: `dec`, `Dec`, `DEJ2000`, `dec_deg`, `_DE`
- Name: `name`, `Name`, `ID`, `Galaxy`, `Object`

### Schritt 4 — In metageometra_ot_master.py einspeisen

Füge die Daten als `external_catalog` Variable ein:

```python
external_catalog = [
    (row['name'], row['ra'], row['dec'])  # name, RA_deg, Dec_deg
    for _, row in df.iterrows()
]
```

Dann rufe auf:
```python
run_OT05(catalog_override=external_catalog)
run_OT06(catalog=external_catalog)
```

### Schritt 5 — Script ausführen

```bash
python metageometra_ot_master.py
```

---

## Koordinaten-Transformation (immer verwenden)

```python
import math

def equatorial_to_galactic(ra_deg, dec_deg):
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    ra_ngp  = math.radians(192.859508)
    dec_ngp = math.radians(27.128336)
    l_ncp   = math.radians(122.932)
    sin_b = (math.sin(dec) * math.sin(dec_ngp) +
             math.cos(dec) * math.cos(dec_ngp) * math.cos(ra - ra_ngp))
    b = math.asin(max(-1.0, min(1.0, sin_b)))
    x = math.cos(dec) * math.sin(ra - ra_ngp)
    y = (math.sin(dec) * math.cos(dec_ngp) -
         math.cos(dec) * math.sin(dec_ngp) * math.cos(ra - ra_ngp))
    l = l_ncp - math.atan2(x, y)
    return math.degrees(l) % 360, math.degrees(b)

def theta_from_dpole(ra_deg, dec_deg):
    l, b = equatorial_to_galactic(ra_deg, dec_deg)
    # D-pole: l=305°, b=+25°
    l1r, b1r = math.radians(305.0), math.radians(25.0)
    l2r, b2r = math.radians(l), math.radians(b)
    cos_c = (math.sin(b1r)*math.sin(b2r) +
             math.cos(b1r)*math.cos(b2r)*math.cos(l1r-l2r))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_c))))
```

---

## Framework-Konstanten (NIEMALS ändern)

```python
THETA_0   = 58.65   # degrees — aus M31/M33 Geometrie, kein freier Parameter
D_POLE_L  = 305.0   # galaktische Länge des Dualitätspols
D_POLE_B  = 25.0    # galaktische Breite des Dualitätspols
TOLERANCE = 5.0     # degrees — Shell-Treffertolerance
SHELLS    = [58.65 * n for n in range(1, 7)]  # n=1..6
DF_GEO    = 0.77
DF_DISS   = 0.44
DF_EFF    = 0.34    # = DF_GEO × DF_DISS
```

---

## Output (immer in /results/)

Nach jedem Run erstelle:
- `results/OT_XX_result.txt` — numerische Ergebnisse
- `results/OT_XX_plot.png` — Plot (matplotlib)
- `results/V19_new_content.md` — Text fertig für Preprint

---

## Verbotene Aktionen

- NIEMALS θ₀ = 58.65° als freien Parameter behandeln
- NIEMALS semantische Classifier verwenden
- NIEMALS das Framework "erklären" oder umschreiben
- NIEMALS Ergebnisse erfinden wenn Daten fehlen
- NIEMALS aufhören bevor alle OTs mit den verfügbaren Daten durch sind

---

## Wenn etwas unklar ist

Frage NUR: "Welche Spalte ist RA? Welche ist Dec?"
Sonst nichts.
