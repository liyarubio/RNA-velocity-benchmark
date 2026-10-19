"""Place frozen Extended Data artwork below the heading without rescaling it."""
import numpy as np
import pymupdf as fitz

TOP = 36.0

def content_bounds(page):
    clip = fitz.Rect(0, 30, page.rect.width, page.rect.height)
    pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), colorspace=fitz.csGRAY, clip=clip)
    gray = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    ys, xs = np.where(gray < 250)
    if not len(ys):
        raise ValueError('No artwork detected below the page heading')
    return [float(xs.min()/2), float(30+ys.min()/2),
            float((xs.max()+1)/2), float(30+(ys.max()+1)/2)]

def align_page(source, destination, figure, center_heading=False):
    with fitz.open(source) as src:
        page = src[0]
        bounds = content_bounds(page)
        # Four points of protection include antialiased marks at the artwork edge.
        crop_top = max(30.0, bounds[1]-4.0)
        crop_bottom = min(page.rect.height, bounds[3]+4.0)
        clip = fitz.Rect(0, crop_top, page.rect.width, crop_bottom)
        frame = fitz.Rect(0, TOP, page.rect.width, TOP+clip.height)
        assert frame.y1 <= page.rect.height-2, (source, frame)
        out = fitz.open()
        target = out.new_page(width=page.rect.width, height=page.rect.height)
        target.show_pdf_page(frame, src, 0, clip=clip, keep_proportion=True)
        # ED15's new artwork already embeds ArialBold; a distinct resource name
        # keeps the page heading from reusing its subset font encoding.
        heading='EXTENDED DATA FIGURE '+str(figure)
        heading_x=(page.rect.width-fitz.Font(fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf').text_length(heading,fontsize=8))/2 if center_heading else 24
        target.insert_text((heading_x,20), heading, fontname='ED15Heading' if figure==15 else 'ArialBold', fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf', fontsize=8)
        out.subset_fonts()
        tmp = destination.with_suffix('.tmp.pdf')
        out.save(tmp, garbage=4, deflate=True)
        tmp.replace(destination)
        target.get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(destination.with_suffix('.png'))
        destination.with_suffix('.svg').write_text(target.get_svg_image(text_as_path=False))
        current = content_bounds(target)
        record = {'figure':figure,'source_bounds_pt':bounds,'current_bounds_pt':current,
                  'source_clip_pt':list(clip),'target_frame_pt':list(frame),
                  'translation_y_pt':TOP-crop_top,'scale':1.0}
        out.close()
        return record

def verify_alignment(source, target, record):
    with fitz.open(source) as src, fitz.open(target) as dst:
        a, b = src[0], dst[0]
        assert abs(a.rect.width-b.rect.width)<.01 and abs(a.rect.height-b.rect.height)<.01
        bounds = content_bounds(b)
        assert abs(bounds[1]-40)<1.1, (target, bounds)
        pa = a.get_pixmap(matrix=fitz.Matrix(2,2), clip=fitz.Rect(record['source_clip_pt']))
        pb = b.get_pixmap(matrix=fitz.Matrix(2,2), clip=fitz.Rect(record['target_frame_pt']))
        assert (pa.width,pa.height)==(pb.width,pb.height), (source, target)
        av=np.frombuffer(pa.samples,dtype=np.uint8).astype(float)
        bv=np.frombuffer(pb.samples,dtype=np.uint8).astype(float)
        error=float(np.abs(av-bv).mean())
        assert error<.5, (target,error)
        # Existing image content may be historical raster; vector-only sources stay vector.
        assert len(a.get_images())==len(b.get_images()), target
        return {'figure':record['figure'],'file':str(target.name),
                'artwork_top_pt':bounds[1],'scale':1.0,
                'mean_pixel_difference_after_translation':error}
