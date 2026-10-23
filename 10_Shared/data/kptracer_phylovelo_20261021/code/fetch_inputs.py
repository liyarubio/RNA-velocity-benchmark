from pathlib import Path
import requests, json, hashlib, tarfile, time
R=Path(__file__).resolve().parents[1]
meta=json.loads((R/'records/zenodo.json').read_text())['files'][0]
dest=R/'inputs'/meta['key'];part=dest.with_suffix(dest.suffix+'.part')
if not dest.exists():
 start=part.stat().st_size if part.exists() else 0
 response=requests.get(meta['links']['self'],headers={'Range':f'bytes={start}-'} if start else {},stream=True,timeout=(30,90))
 response.raise_for_status()
 mode='ab' if start and response.status_code==206 else 'wb'
 if mode=='wb':start=0
 total=start;last=time.monotonic()
 with part.open(mode) as out:
  for chunk in response.iter_content(4*1024*1024):
   if not chunk:continue
   out.write(chunk);total+=len(chunk)
   if time.monotonic()-last>15:
    print(f'Downloaded {total/1e6:.1f} / {meta["size"]/1e6:.1f} MB',flush=True);last=time.monotonic()
 assert part.stat().st_size==meta['size']
 part.rename(dest)
h=hashlib.md5()
with dest.open('rb') as f:
 for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
assert 'md5:'+h.hexdigest()==meta['checksum']
print('Archive MD5 verified',h.hexdigest(),flush=True)
entries=[];extracted=[]
with tarfile.open(dest,'r:gz') as tar:
 for member in tar:
  if not member.isfile():continue
  entries.append({'name':member.name,'size':member.size})
  if Path(member.name).name in ['3726_NT_T1_tree.nwk','adata_processed.nt.h5ad']:
   target=R/'inputs'/Path(member.name).name
   source=tar.extractfile(member)
   sha=hashlib.sha256()
   with target.open('wb') as out:
    for chunk in iter(lambda:source.read(8*1024*1024),b''):out.write(chunk);sha.update(chunk)
   extracted.append({'member':member.name,'file':str(target.relative_to(R)),'sha256':sha.hexdigest(),'bytes':target.stat().st_size})
   print('Extracted',member.name,member.size,flush=True)
(R/'records/archive_members.json').write_text(json.dumps(entries,indent=2))
(R/'records/download.json').write_text(json.dumps({'source_record':'https://zenodo.org/records/5847462','url':meta['links']['self'],'archive_md5':h.hexdigest(),'archive_bytes':dest.stat().st_size,'extracted':extracted},indent=2))
assert len(extracted)==2,extracted
