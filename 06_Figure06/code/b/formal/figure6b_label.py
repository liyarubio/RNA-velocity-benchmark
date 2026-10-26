"""Clarify the formal Figure 6b axis while preserving all saved curves."""
from pathlib import Path
import argparse, json
import pymupdf as fitz
import numpy as np

HERE = Path(__file__).resolve().parent
LABEL = 'Mean velocity cosine similarity to full-depth prediction'
LINES = ['Mean velocity cosine similarity', 'to full-depth prediction']

def relabel(page):
    spans = [s for b in page.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]
    hits = [s for s in spans if s['text'] == 'Mean cosine to 1.0× reads']
    if not hits and any(s['text'].replace('\u00a0', ' ') == LINES[0] for s in spans):
        return None
    assert len(hits) == 1, 'Expected exactly one original read-depth axis label'
    old = hits[0]; box = fitz.Rect(old['bbox']); size = old['size']
    center = (box.y0 + box.y1) / 2
    page.add_redact_annot(box + (-.3, -.3, .3, .3), fill=(1, 1, 1))
    page.apply_redactions(images=0, graphics=0)
    # Liberation Sans matches the existing figure; Helvetica is a portable fallback.
    fontpath = Path('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf')
    font = fitz.Font(fontfile=str(fontpath)) if fontpath.exists() else fitz.Font('helv')
    page.insert_font(fontname='ReadDepthAxis', fontbuffer=font.buffer)
    bounds = fitz.Rect(box)
    for i, line in enumerate(LINES):
        x = old['origin'][0] + (i - 1) * size * 1.15
        width = font.text_length(line, fontsize=size)
        y = center + width / 2
        page.insert_text((x, y), line, fontsize=size, fontname='ReadDepthAxis', rotate=90)
        bounds |= fitz.Rect(x - size * 1.1, center - width / 2 - 1, x + size * .3, y + 1)
    return bounds + (-1, -1, 1, 1)

def pixels(page):
    pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
    return np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3)

def render(source=None):
    source = source or fitz.open(HERE / 'assets/Figure06_before.pdf')
    doc = fitz.open(); doc.insert_pdf(source, from_page=0, to_page=0)
    bounds = relabel(doc[0])
    assert bounds is not None
    doc = fitz.open(stream=doc.tobytes(garbage=4, deflate=True), filetype='pdf')
    before, after = pixels(source[0]), pixels(doc[0])
    mask = np.ones(before.shape[:2], bool)
    x0, y0, x1, y1 = bounds
    mask[max(0, int(y0*2)-2):int(y1*2)+3, max(0, int(x0*2)-2):int(x1*2)+3] = False
    assert np.array_equal(before[mask], after[mask]), 'Pixels changed outside the axis label'
    text = doc[0].get_text().replace('\u00a0', ' ').replace('\u00ad', '-')
    assert all(line in text for line in LINES) and 'Mean cosine to 1.0× reads' not in text
    return doc, dict(status='pass', label=LABEL, changed_bounds=list(bounds),
                     other_artwork_pixels_identical=True, numerical_values_changed=False)

def export(doc, stem, crop=None):
    stem = Path(stem); stem.parent.mkdir(parents=True, exist_ok=True)
    preview = doc
    if crop:
        preview = fitz.open(); page = preview.new_page(width=crop.width, height=crop.height)
        page.show_pdf_page(page.rect, doc, 0, clip=crop)
    preview.save(stem.with_suffix('.pdf'), garbage=4, deflate=True)
    preview[0].get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False).save(stem.with_suffix('.png'))
    stem.with_suffix('.svg').write_text(preview[0].get_svg_image(text_as_path=True))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); doc, audit = render()
    export(doc, args.out / 'Figure06')
    export(doc, args.out / 'Figure06b_cosine', fitz.Rect(300, 115, 590, 403))
    (args.out / 'axis_validation.json').write_text(json.dumps(audit, indent=2)+'\n')

if __name__ == '__main__':
    main()
