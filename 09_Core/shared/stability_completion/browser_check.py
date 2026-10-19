"""Check all four local readers using the temporary Firefox BiDi session."""
import base64
import hashlib
import json
import os
from pathlib import Path
import websocket

LOG = Path(__file__).resolve().parent
ROOT = LOG.parent
AUDIT = ROOT.parent / 'audit_manuscript_20260930'
for key in list(os.environ):
    if key.lower().endswith('_proxy'):
        os.environ.pop(key)
ws = websocket.create_connection('ws://127.0.0.1:9268/session', timeout=45, suppress_origin=True)
seq = 0


def call(method, params):
    global seq
    seq += 1
    ws.send(json.dumps({'id': seq, 'method': method, 'params': params}))
    while True:
        response = json.loads(ws.recv())
        if response.get('id') == seq:
            assert response.get('type') != 'error', response
            return response['result']


call('session.new', {'capabilities': {}})
ctx = call('browsingContext.getTree', {})['contexts'][0]['context']


def js(expression):
    result = call('script.evaluate', {'expression': expression, 'target': {'context': ctx}, 'awaitPromise': True})
    assert result.get('type') != 'exception', result
    return result['result'].get('value')


checks = {}
try:
    call('browsingContext.setViewport', {'context': ctx, 'viewport': {'width': 1400, 'height': 1050}})
    paths = [ROOT / 'Manuscript_0528_minimal_bilingual.html', ROOT / 'Manuscript_0528_minimal_bilingual_standalone.html',
             AUDIT / 'manuscript_nature_methods_full.html', AUDIT / 'manuscript_nature_methods_standalone.html']
    for path in paths:
        call('browsingContext.navigate', {'context': ctx, 'url': path.as_uri(), 'wait': 'complete'})
        js('Promise.all([...document.images].map(i=>{i.loading="eager";return i.decode()})).then(()=>true)')
        row = json.loads(js('JSON.stringify({pairs:document.querySelectorAll(".pair").length,english:document.querySelectorAll(".en").length,chinese:document.querySelectorAll(".zh").length,loaded:[...document.images].filter(i=>i.complete&&i.naturalWidth>0).length,overflow:document.documentElement.scrollWidth>innerWidth,additions:document.querySelectorAll(".algorithm-added").length,brokenAnchors:[...document.querySelectorAll("a[href]")].filter(a=>a.getAttribute("href").startsWith("#")&&!document.getElementById(a.hash.slice(1))).length})'))
        assert row == {'pairs': 203, 'english': 203, 'chinese': 195, 'loaded': 7, 'overflow': False, 'additions': 0, 'brokenAnchors': 0}, row
        row['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        checks[str(path)] = row
    for n in range(3,7):
        js(f'document.getElementById("ppt-figure-{n}").scrollIntoView();true')
        capture=call('browsingContext.captureScreenshot', {'context':ctx})
        (LOG/f'html_figure{n}.png').write_bytes(base64.b64decode(capture['data']))
    assert js('document.documentElement.dataset.figureDataRevision') == 'figures-and-data-20261008'
    assert js('document.querySelectorAll(".figure-data-updated").length') == 6
    assert 'noise_' not in js('document.querySelector("#ppt-figure-6 figcaption").textContent')
    positions = json.loads(js('JSON.stringify([...document.querySelectorAll("#all-figures-gallery > figure")].map(f=>f.id))'))
    assert positions == ['ppt-figure-'+str(n) for n in range(1,8)]
    checks['gallery_figure_order'] = positions
    assert js('document.querySelectorAll("figure[data-placement-after]").length') == 0
    js('document.querySelector("#P041 .en .figure-ref").click();true')
    assert js('location.hash') == '#ppt-figure-3'
    js('document.querySelector("[data-lang=en]").click();document.getElementById("P036").scrollIntoView();true')
    assert js('getComputedStyle(document.querySelector("#P028 .zh")).display') == 'none'
    screenshot = call('browsingContext.captureScreenshot', {'context': ctx})
    (LOG / 'reader_figure2.png').write_bytes(base64.b64decode(screenshot['data']))
    red = js('getComputedStyle(document.querySelector("#P051 .en .figure-data-updated")).color')
    js('document.getElementById("red").click();true')
    black = js('getComputedStyle(document.querySelector("#P051 .en .figure-data-updated")).color')
    assert red != black
    js('document.getElementById("red").click();[...document.querySelectorAll("a")].find(a=>a.hash==="#reference-72").click();true')
    assert js('location.hash') == '#reference-72'
    js('document.querySelector("[data-lang=both]").click();document.getElementById("P028").scrollIntoView();true')
    call('browsingContext.setViewport', {'context': ctx, 'viewport': {'width': 390, 'height': 844}})
    assert js('document.documentElement.scrollWidth<=innerWidth')
    checks['interaction'] = {'language_toggle': True, 'redline_toggle': True, 'reference_click': True, 'mobile_no_overflow': True}
    call('browsingContext.navigate', {'context':ctx,'url':(ROOT/'Figures_and_data_20261008.html').as_uri(),'wait':'complete'})
    js('Promise.all([...document.images].map(i=>{i.loading="eager";return i.decode()})).then(()=>true)')
    assert js('document.querySelectorAll("figure").length') == 7
    assert js('document.documentElement.scrollWidth<=innerWidth')
    checks['separate_gallery']={'figures':7,'loaded':7,'mobile_no_overflow':True}
finally:
    call('session.end', {})
    ws.close()
checks['status'] = 'pass'
checks['current_version'] = 'Manuscript_0528_original_style_gallery_20261008'
checks['index_sha256'] = hashlib.sha256((AUDIT / 'index.html').read_bytes()).hexdigest()
payload = json.dumps(checks, ensure_ascii=False, indent=2)
(LOG / 'browser_validation.json').write_text(payload)
(AUDIT / 'full_html_browser_validation.json').write_text(payload)
print('PASS: four readers plus gallery; all seven figures together; references and mobile layout.')
