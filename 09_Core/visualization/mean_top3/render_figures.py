"""Render both complete figure artwork pages with Mean-only summaries."""
import argparse
import json
from pathlib import Path
from mean_top3 import render_figure, save_artwork


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('demo_output/mean_top3'))
    args = parser.parse_args()
    for name,panels in [('Figure03',['Figure03a']),('Extended_Data_Fig04',['ED04a','ED04b'])]:
        doc,audit = render_figure(panels)
        save_artwork(doc,args.out/name)
        (args.out/(name+'_validation.json')).write_text(json.dumps(audit,indent=2)+'\n')
        print(name+': passed',flush=True)


if __name__ == '__main__':
    main()
