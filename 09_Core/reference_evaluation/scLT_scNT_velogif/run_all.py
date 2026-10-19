"""Durable sequential run, with per-stage logs and status updates."""
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/task')


def main():
    states = []
    for dataset in ('scLT', 'scNT'):
        config = str(ROOT / 'configs' / (dataset + '.yaml'))
        stages = [('check', ['-m', 'velogif', 'check', '--config', config]),
                  ('run', ['-m', 'velogif', 'run', '--config', config, '--no-plot']),
                  ('plot', ['-m', 'velogif', 'plot', '--config', config]),
                  ('reference', [str(ROOT / 'scripts/evaluate_reference.py'), '--dataset', dataset])]
        for stage, args in stages:
            path = ROOT / 'records' / (dataset + '_' + stage + '.log')
            if path.exists():
                path.rename(path.with_name(path.name + '.' + datetime.now().strftime('%Y%m%dT%H%M%S')))
            environment = 'scvelo' if stage == 'reference' else 'cuda_env'
            command = ['conda', 'run', '--no-capture-output', '-n', environment, 'python'] + args
            item = {'dataset': dataset, 'stage': stage, 'started': datetime.now(timezone.utc).isoformat(),
                    'status': 'running', 'log': str(path), 'command': command}
            states.append(item)
            (ROOT / 'records/pipeline.json').write_text(json.dumps(states, indent=2))
            print('START', dataset, stage, flush=True)
            with path.open('w') as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                        cwd='/task/code', env=dict(os.environ))
            item.update(status='completed' if result.returncode == 0 else 'failed',
                        exit_code=result.returncode, ended=datetime.now(timezone.utc).isoformat())
            (ROOT / 'records/pipeline.json').write_text(json.dumps(states, indent=2))
            print('END', dataset, stage, result.returncode, flush=True)
            if stage == 'check' and result.returncode:
                break
    return int(any(s['status'] == 'failed' for s in states))


if __name__ == '__main__':
    raise SystemExit(main())
