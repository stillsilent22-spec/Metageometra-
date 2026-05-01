"""Pre-registration check: Extract L/a0 definitions from Metageometra version docs"""
import os, sys, re

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

DOWNLOADS = r"c:\Users\kalle\Downloads"

def search_text(text, label):
    patterns = [
        (r'seepage', 'seepage'),
        (r'L_seepage|L\s+seepage', 'L_seepage'),
        (r'global\s+(?:leak|seepage|rate|ooze)', 'global-rate'),
        (r'shell.local|local\s+seepage|per.shell', 'shell-local'),
        (r'a_?0\s*=\s*c\s*/|c\s*/\s*2.*pi.*T|Grundbeschleunigung', 'a0-def'),
        (r'Versickerungsrate|Sickerstrom|Leckrate', 'German-seepage'),
        (r'was defined|defined as|defined in|first appear|introduced in', 'first-def'),
        (r'L\s*=\s*\\?rho|L\s*=\s*\d|L\s*wird', 'L-assignment'),
        (r'v0\.2|v4\.0|version\s+0\.2|version\s+4', 'version-ref'),
    ]
    found = []
    for pat, tag in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE | re.MULTILINE):
            ctx = text[max(0,m.start()-120):m.start()+250]
            found.append((tag, m.start(), ctx.strip()))
    if found:
        print(f"\n=== {label}: {len(found)} hits ===")
        for tag, pos, ctx in found[:25]:
            print(f"  [{tag}]: {ctx[:220]}")
            print()
    else:
        print(f"\n=== {label}: no hits ===")
    return found

def try_pdf(path):
    try:
        import fitz
        doc = fitz.open(path)
        text = ''.join(page.get_text() for page in doc)
        print(f"  fitz OK: {len(text)} chars")
        return text
    except Exception:
        pass
    try:
        import pdfplumber
        with pdfplumber.open(path) as pdf:
            text = ''.join(p.extract_text() or '' for p in pdf.pages)
        print(f"  pdfplumber OK: {len(text)} chars")
        return text
    except Exception as e:
        print(f"  PDF read error: {e}")
    return ""

def try_docx(path):
    try:
        import docx
        doc = docx.Document(path)
        text = '\n'.join(p.text for p in doc.paragraphs)
        print(f"  docx OK: {len(text)} chars")
        return text
    except ImportError:
        pass
    try:
        import zipfile, xml.etree.ElementTree as ET
        with zipfile.ZipFile(path) as z:
            with z.open('word/document.xml') as f:
                tree = ET.parse(f)
        ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        text = ' '.join(t.text for t in tree.findall('.//w:t', ns) if t.text)
        print(f"  docx-zip OK: {len(text)} chars")
        return text
    except Exception as e:
        print(f"  DOCX read error: {e}")
        return ""

docs = {
    "V18_PDF":   os.path.join(DOWNLOADS, "Metageometra_V18_full (2).pdf"),
    "V11_PDF":   os.path.join(DOWNLOADS, "Metageometra_V11_Master_Synthesis.pdf"),
    "V10_DOCX":  os.path.join(DOWNLOADS, "Metageometra_V10_Master_Synthesis.docx"),
    "Pos_Paper": os.path.join(DOWNLOADS, "Hannemann_Torsion_Position_Paper.docx"),
    "Echo_PDF":  os.path.join(DOWNLOADS, "-Kopie-Echo-Inversion_im_Hannemann-Torsionsmodell_(HTM).pdf"),
    "Fraktal_DOCX": os.path.join(DOWNLOADS, "Fraktale_Noether_Dissipation.docx"),
}

all_found = {}
for label, path in docs.items():
    print(f"\nProcessing {label}: {os.path.basename(path)}")
    if not os.path.exists(path):
        print("  FILE NOT FOUND")
        continue
    text = try_pdf(path) if path.endswith('.pdf') else try_docx(path)
    if text:
        all_found[label] = search_text(text, label)

for fname in ["COPILOT_PROMPT.md", "V19_new_content.md"]:
    fpath = os.path.join(DOWNLOADS, fname)
    if os.path.exists(fpath):
        with open(fpath, encoding='utf-8', errors='replace') as f:
            text = f.read()
        all_found[fname] = search_text(text, fname)

print("\n" + "="*70)
print("SUMMARY:")
for label, hits in all_found.items():
    print(f"  {label}: {[h[0] for h in hits]}")

