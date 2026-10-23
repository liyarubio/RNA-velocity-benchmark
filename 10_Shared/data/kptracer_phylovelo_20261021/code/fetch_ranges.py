"""Resume the existing prefix and fetch the remainder using bounded HTTP ranges."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests,json,time,shutil,subprocess,sys
R=Path(__file__).resolve().parents[1]
meta=json.loads((R/'records/zenodo.json').read_text())['files'][0]
dest=R/'inputs'/meta['key'];part=dest.with_suffix(dest.suffix+'.part')
if dest.exists():
 subprocess.run([sys.executable,str(R/'code/fetch_inputs.py')],check=True)
 raise SystemExit(0)
prefix=R/'inputs/download_prefix.part';parts=R/'inputs/download_parts';parts.mkdir(exist_ok=True)
if not prefix.exists() and part.exists():part.rename(prefix)
start=prefix.stat().st_size if prefix.exists() else 0
chunk=32*1024*1024
ranges=[(n,a,min(a+chunk,meta['size'])-1) for n,a in enumerate(range(start,meta['size'],chunk))]
def fetch(job):
 n,a,b=job;p=parts/f'{n:04d}_{a}_{b}.part'
 if p.exists() and p.stat().st_size==b-a+1:return p
 for attempt in range(4):
  try:
   response=requests.get(meta['links']['self'],headers={'Range':f'bytes={a}-{b}'},stream=True,timeout=(30,90))
   response.raise_for_status()
   assert response.status_code==206 and response.headers.get('Content-Range')==f'bytes {a}-{b}/{meta["size"]}'
   with p.open('wb') as out:
    for data in response.iter_content(1024*1024):out.write(data)
   assert p.stat().st_size==b-a+1
   return p
  except Exception:
   if attempt==3:raise
   time.sleep(2+attempt*2)
with ThreadPoolExecutor(max_workers=8) as pool:
 done=start
 for fut in as_completed([pool.submit(fetch,j) for j in ranges]):
  p=fut.result();done+=p.stat().st_size
  print(f'Complete ranges: {done/1e6:.1f} / {meta["size"]/1e6:.1f} MB',flush=True)
with part.open('wb') as out:
 if prefix.exists():
  with prefix.open('rb') as f:shutil.copyfileobj(f,out,8*1024*1024)
 for n,a,b in ranges:
  with (parts/f'{n:04d}_{a}_{b}.part').open('rb') as f:shutil.copyfileobj(f,out,8*1024*1024)
assert part.stat().st_size==meta['size']
part.rename(dest)
subprocess.run([sys.executable,str(R/'code/fetch_inputs.py')],check=True)
