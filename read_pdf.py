import pypdf
pdf_path = r'C:\Users\kalle\Neuer Ordner (4)\v6_1_framework.pdf'
reader = pypdf.PdfReader(pdf_path)
print(f"SEITEN: {len(reader.pages)}")
for i, page in enumerate(reader.pages):
    txt = page.extract_text() or "(keine text-layer)"
    print(f"\n{'='*60}")
    print(f"SEITE {i+1}")
    print('='*60)
    print(txt)
