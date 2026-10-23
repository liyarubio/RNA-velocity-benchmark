"""Run independent saved-data panel files; contains no plotting implementation."""
from pathlib import Path
import argparse, json, subprocess, sys
ROOT=Path(__file__).resolve().parents[2]
REGISTRY=Path(__file__).with_name('panel_registry.json')

def main():
    entries=json.loads(REGISTRY.read_text())
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'panel_checks')
    parser.add_argument('--panel',action='append',help='Select Figure06/a, ED05/c, etc.; repeat to select several.')
    parser.add_argument('--list',action='store_true')
    args=parser.parse_args()
    available={e['figure']+'/'+e['panel']:e for e in entries}
    if args.list:
        for key,e in available.items():print(key+'\t'+e['path']+'\t'+e['scope'])
        return
    unknown=set(args.panel or [])-set(available)
    if unknown:parser.error('Unknown panels: '+', '.join(sorted(unknown)))
    selected=[available[k] for k in args.panel] if args.panel else entries
    if len({e['output'] for e in selected})!=len(selected):
        parser.error('Selected panel output names must be unique')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True);log=[]
    for e in selected:
        print('Rendering '+e['figure']+'/'+e['panel'],flush=True)
        result=subprocess.run([sys.executable,str(ROOT/e['path']),'--out',str(out)],text=True,capture_output=True)
        if result.returncode:
            sys.stderr.write(result.stderr);raise RuntimeError(e['path']+' failed: '+result.stdout)
        records=json.loads(result.stdout)
        if len(records)!=1 or records[0]['output']!=e['output']:
            raise RuntimeError(e['path']+' must render exactly its registered plot')
        log.extend(dict(row,figure=e['figure'],panel=e['panel'],entrypoint=e['path'],scope=e['scope']) for row in records)
    (out/'plot_manifest.json').write_text(json.dumps(dict(status='complete',plots=log,models_rerun=False,publication_artwork_modified=False),indent=2)+'\n')
    print(str(len(log))+' independent panel redraws written to '+str(out))

if __name__=='__main__':main()
