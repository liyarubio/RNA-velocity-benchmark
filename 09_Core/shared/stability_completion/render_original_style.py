"""Refresh only original-style Figures 5 and 6, preserving native PDF artwork."""
from pathlib import Path
import sys,json,hashlib,xml.etree.ElementTree as ET
R=Path(__file__).resolve().parent.parent
D=R/'original_style_20261008'
p=D/'build_original_style.py'
s=p.read_text()
prefix=s[:s.index('# Restore the workflow artwork')]
start=s.index('# Figure 5:')
end=s.index("subprocess.run([sys.executable,str(D/'preserve_native_pdf.py')]")
exec(compile(prefix+s[start:end],str(p),'exec'),{'__file__':str(p),'__name__':'stability_original_style'})
# Retain PDF image masks, transparency and framework artwork. Do not regenerate Figure2.
p=D/'preserve_native_pdf.py';s=p.read_text()
prefix=s[:s.index('# Figure 2:')]
tail=s[s.index('# Figure 6b:'):]
exec(compile(prefix+tail,str(p),'exec'),{'__file__':str(p),'__name__':'stability_native_pdf'})
path=D/'figure_manifest.json'
manifest=json.loads(path.read_text())
for n in [5,6]:
    p=R/f'Figures_original_style/Figure{n}/Figure{n}_original_style.svg'
    root=ET.parse(p).getroot()
    manifest[str(n)].update(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        svg_images=len(root.findall('.//{http://www.w3.org/2000/svg}image')),
        stability_completion='20261008')
path.write_text(json.dumps(manifest,indent=2)+'\n')
