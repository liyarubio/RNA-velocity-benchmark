"""Publish numbered Extended Data assets without changing numerical source IDs."""
from pathlib import Path
import csv, hashlib, html, json, shutil, sys
from bs4 import BeautifulSoup
from figure_order_1019 import NUMBER, LEGACY_ORDER, apply_figure_order
from figure2_legend_style import draw_caption
sys.path.insert(0, '/data1/liyaru/RNA_velocity_materials/Manuscript_1010/Code/vendor')
import pymupdf as fitz

B = Path(__file__).resolve().parents[1]
R = B / 'Reports/figure_order_20261019'
BEFORE = B.parent / 'Manuscript_1019_revisions/figure_order_20261019_before'
FONT = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, x): p.write_text(json.dumps(x, ensure_ascii=False, indent=2) + '\n')
def normalized(t): return ' '.join(t.replace('\xad', '-').split())

def update_main_cross_reference(doc):
    """Main Figure 2 also cites the renumbered simulation comparison."""
    outputs = {}; checked = []
    caption = doc.find(id='ppt-figure-2').figcaption.get_text(' ', strip=True)
    for rel in ['submission_figures/Figure02.pdf', 'figures/main/Figure02.pdf', 'submission_figures/Main_Figures_1-7_combined.pdf']:
        backup=BEFORE/rel
        if not backup.exists():
            backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(B/rel,backup)
    for rel in ['submission_figures/Figure02.pdf', 'figures/main/Figure02.pdf']:
        old=fitz.open(BEFORE/rel);out=fitz.open();out.insert_pdf(old)
        page=out[0]
        block=next(b for b in page.get_text('blocks') if normalized(b[4]).startswith('Fig. 2 |'))
        box=fitz.Rect(block[:4]);clip=fitz.Rect(0,0,page.rect.width,box.y0-2)
        page.add_redact_annot(box+(-1,-1,1,1),fill=(1,1,1));page.apply_redactions(images=0,graphics=0)
        style=draw_caption(page,fitz.Rect(32,box.y0,page.rect.width-32,page.rect.height-32),caption,fitz)
        assert old[0].get_pixmap(clip=clip,alpha=False).samples==page.get_pixmap(clip=clip,alpha=False).samples
        out.set_metadata({'title':'Figure 2 | Schematic overview of the benchmarking workflow','subject':'Cross-reference updated to Extended Data Fig. 4c; scientific artwork preserved'})
        tmp=(B/rel).with_suffix('.tmp.pdf');out.save(tmp,garbage=3,deflate=True);tmp.replace(B/rel)
        checked.append({'path':rel,'artwork_pixels_identical':True,'clip':list(clip),'caption_style':style})
        outputs[rel]=sha(B/rel)
    rel='submission_figures/Main_Figures_1-7_combined.pdf'
    old=fitz.open(BEFORE/rel);out=fitz.open();updated=fitz.open(B/'submission_figures/Figure02.pdf')
    out.insert_pdf(old,from_page=0,to_page=0);out.insert_pdf(updated);out.insert_pdf(old,from_page=2)
    out.set_toc(old.get_toc());out.set_metadata(old.metadata)
    for i in [0]+list(range(2,len(old))):
        assert old[i].get_pixmap(alpha=False).samples==out[i].get_pixmap(alpha=False).samples
    tmp=(B/rel).with_suffix('.tmp.pdf');out.save(tmp,garbage=3,deflate=True);tmp.replace(B/rel)
    outputs[rel]=sha(B/rel)
    return {'artwork_checks':checked,'combined_unchanged_pages':[1]+list(range(3,len(old)+1)),
        'updated_cross_reference':'Extended Data Fig. 5c -> Extended Data Fig. 4c'},outputs

def main():
    R.mkdir(exist_ok=True)
    artwork = B / 'figures/extended_numbered'
    final = B / 'submission_figures/numbered_extended'
    legends = final / 'legends'
    for path in [artwork, final, legends]: path.mkdir(parents=True, exist_ok=True)
    doc = BeautifulSoup((BEFORE/'Manuscript_1019_NM_bilingual.html').read_text(), 'html.parser')
    apply_figure_order(doc)
    rows = json.loads((BEFORE/'figure_data_code_map.json').read_text())
    merged = fitz.open(); toc = []; validation = []; mapping = []
    for old in LEGACY_ORDER:
        new = NUMBER[old]
        key = f'Extended_Data_Fig{new:02d}'
        oldkey = f'Extended_Data_Fig{old:02d}'
        section = doc.find(id=f'ed{new}')
        title = section.h3.get_text(' ', strip=True)
        caption = section.figcaption.find(lang='en').get_text(' ', strip=True)
        for ext in ['png', 'svg']:
            source = B / f'figures/extended/{oldkey}.{ext}'
            target = artwork / f'{key}.{ext}'
            shutil.copy2(source, target)
            assert sha(source) == sha(target)
        original = fitz.open(B/f'submission_figures/{oldkey}.pdf')
        assert len(original) == 1
        output = fitz.open(); output.insert_pdf(original)
        page = output[0]
        blocks = [b for b in page.get_text('blocks') if normalized(b[4]).startswith('Extended Data Fig.')]
        assert len(blocks) <= 1
        clip = None
        if blocks:
            # Two inherited figures place their caption on the artwork page.
            rect = fitz.Rect(blocks[0][:4])
            clip = fitz.Rect(0, 0, page.rect.width, rect.y0 - 2)
            page.add_redact_annot(rect + (-1, -1, 1, 1), fill=(1, 1, 1))
            page.apply_redactions(images=0, graphics=0)
            box = fitz.Rect(36, rect.y0, page.rect.width-36, page.rect.height-32)
            size = 9.5
            result = page.insert_textbox(box, title + '. ' + caption, fontname='CaptionArial', fontfile=FONT, fontsize=size, lineheight=1.25)
            assert result >= 0, (old, result)
        else:
            legend = output.new_page(width=page.rect.width, height=page.rect.height)
            result = legend.insert_textbox(fitz.Rect(36, 42, legend.rect.width-36, legend.rect.height-36),
                title + '\n\n' + caption, fontname='CaptionArial', fontfile=FONT, fontsize=10, lineheight=1.35)
            assert result >= 0, (old, result)
        # This compares every pixel of the scientific artwork at original placement.
        pixel_identical = original[0].get_pixmap(clip=clip, alpha=False).samples == output[0].get_pixmap(clip=clip, alpha=False).samples
        assert pixel_identical, old
        output.set_metadata({'title': title, 'subject': f'Manuscript 1019; current ED{new:02d}; legacy source ED{old:02d}; numbering only'})
        output.set_toc([[1, title, 1]])
        output.save(final/f'{key}.pdf', garbage=3, deflate=True)
        shutil.copy2(final/f'{key}.pdf', artwork/f'{key}.pdf')
        legend_doc = fitz.open(); legend_doc.insert_pdf(output, from_page=len(output)-1, to_page=len(output)-1)
        if not blocks: legend_doc.save(legends/f'{key}_legend.pdf', garbage=3, deflate=True)
        toc.append([1, title, len(merged)+1]); merged.insert_pdf(output)
        validation.append({'old':old,'new':new,'artwork_pixels_identical':pixel_identical,
            'artwork_comparison_clip':list(clip) if clip else 'entire original artwork page',
            'png_identical':True,'svg_identical':True,'pages':len(output),'caption':caption})
        row = next(r for r in rows if r['figure'] == f'ED{old:02d}')
        row.update(figure=f'ED{new:02d}', legacy_figure=f'ED{old:02d}',
            artwork=f'figures/extended_numbered/{key}', final_pdf=f'submission_figures/numbered_extended/{key}.pdf')
        mapping.append({'old':old,'new':new,'title':title,'source_data':row['source_data'],
            'artwork':row['artwork'],'final_pdf':row['final_pdf']})
    merged.set_toc(toc)
    merged.set_metadata({'title':'Manuscript 1019 — Extended Data Figures 1–12 in first-citation order'})
    target = B/'submission_figures/Extended_Data_Figures_combined.pdf'
    temporary = target.with_suffix('.tmp.pdf'); merged.save(temporary, garbage=3, deflate=True); temporary.replace(target)
    for row in rows:
        if row['figure'].startswith('Figure'):
            row.update(legacy_figure=row['figure'],artwork=f'figures/main/{row["figure"]}',final_pdf=f'submission_figures/{row["figure"]}.pdf')
    dump(B/'figure_data_code_map.json', rows)
    with (B/'figure_data_code_map.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    dump(R/'mapping.json', mapping)
    main_check,main_outputs=update_main_cross_reference(doc)
    dump(R/'asset_validation.json', {'status':'pass','figures_checked':12,'combined_pages':len(merged),'figures':validation,
        'main_figure2':main_check,
        'outputs_sha256':{'submission_figures/Extended_Data_Figures_combined.pdf':sha(target),**main_outputs},
        'source_id_policy':'Legacy artwork and numerical workbook paths retained; current exports are explicitly numbered.',
        'models_rerun':False,'numeric_inputs_changed':False})
    print(json.dumps({'figures':12,'pages':len(merged),'artwork_pixels_identical':True}))

if __name__ == '__main__': main()
